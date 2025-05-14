#include "band_energy.h"

#include "evolve_elec.h"
#include "module_base/lapack_connector.h"
#include "module_base/module_container/ATen/kernels/blas.h"
#include "module_base/scalapack_connector.h"

#include <complex>
#include <iostream>

namespace module_tddft
{
#ifdef __MPI

inline int globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

void compute_ekb(const Parallel_Orbitals* pv,
                 const int nband,
                 const int nlocal,
                 const std::complex<double>* Htmp,
                 const std::complex<double>* psi_k,
                 double* ekb,
                 std::ofstream& ofs_running)
{
    assert(pv->nloc_wfc > 0 && pv->nloc > 0);

    std::complex<double>* tmp1 = new std::complex<double>[pv->nloc_wfc];
    ModuleBase::GlobalFunc::ZEROS(tmp1, pv->nloc_wfc);

    std::complex<double>* eij = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(eij, pv->nloc);

    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nlocal,
                             1.0,
                             Htmp,
                             1,
                             1,
                             pv->desc,
                             psi_k,
                             1,
                             1,
                             pv->desc_wfc,
                             0.0,
                             tmp1,
                             1,
                             1,
                             pv->desc_wfc);

    ScalapackConnector::gemm('C',
                             'N',
                             nband,
                             nband,
                             nlocal,
                             1.0,
                             psi_k,
                             1,
                             1,
                             pv->desc_wfc,
                             tmp1,
                             1,
                             1,
                             pv->desc_wfc,
                             0.0,
                             eij,
                             1,
                             1,
                             pv->desc_Eij);

    if (PARAM.inp.td_print_eij > 0.0)
    {
        ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
        ofs_running << " Eij:" << std::endl;
        for (int i = 0; i < pv->nrow_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol_bands; j++)
            {
                double aa = eij[in + j].real();
                double bb = eij[in + j].imag();
                if (std::abs(aa) < PARAM.inp.td_print_eij)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < PARAM.inp.td_print_eij)
                {
                    bb = 0.0;
                }
                if (aa > 0.0 || bb > 0.0)
                {
                    ofs_running << i << " " << j << " " << aa << "+" << bb << "i " << std::endl;
                }
            }
        }
        ofs_running << std::endl;
        ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
    }

    int info = 0;
    int naroc[2] = {0, 0};

    assert(nband > 0);
    double* eii = new double[nband];
    ModuleBase::GlobalFunc::ZEROS(eii, nband);

    for (int iprow = 0; iprow < pv->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < pv->dim1; ++ipcol)
        {
            if (iprow == pv->coord[0] && ipcol == pv->coord[1])
            {
                naroc[0] = pv->nrow;
                naroc[1] = pv->ncol;
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, pv->nb, pv->dim1, ipcol);
                    if (igcol >= nband)
                    {
                        continue;
                    }
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, pv->nb, pv->dim0, iprow);
                        if (igrow >= nband)
                        {
                            continue;
                        }
                        if (igcol == igrow)
                        {
                            eii[igcol] = eij[j * naroc[0] + i].real();
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow
    info = MPI_Allreduce(eii, ekb, nband, MPI_DOUBLE, MPI_SUM, pv->comm());

    delete[] tmp1;
    delete[] eij;
    delete[] eii;
}

void compute_ekb_tensor(const Parallel_Orbitals* pv,
                        const int nband,
                        const int nlocal,
                        const ct::Tensor& Htmp,
                        const ct::Tensor& psi_k,
                        ct::Tensor& ekb,
                        std::ofstream& ofs_running)
{
    assert(pv->nloc_wfc > 0 && pv->nloc > 0);

    // Create Tensor objects for temporary data
    ct::Tensor tmp1(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({pv->nloc_wfc}));
    tmp1.zero();

    ct::Tensor eij(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({pv->nloc}));
    eij.zero();

    // Perform matrix multiplication: tmp1 = Htmp * psi_k
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nlocal,
                             1.0,
                             Htmp.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc,
                             psi_k.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_wfc,
                             0.0,
                             tmp1.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_wfc);

    // Perform matrix multiplication: eij = psi_k^dagger * tmp1
    ScalapackConnector::gemm('C',
                             'N',
                             nband,
                             nband,
                             nlocal,
                             1.0,
                             psi_k.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_wfc,
                             tmp1.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_wfc,
                             0.0,
                             eij.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_Eij);

    if (PARAM.inp.td_print_eij >= 0.0)
    {
        ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
        ofs_running << " Eij:" << std::endl;
        for (int i = 0; i < pv->nrow_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol_bands; j++)
            {
                double aa = eij.data<std::complex<double>>()[in + j].real();
                double bb = eij.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < PARAM.inp.td_print_eij)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < PARAM.inp.td_print_eij)
                {
                    bb = 0.0;
                }
                if (aa > 0.0 || bb > 0.0)
                {
                    ofs_running << i << " " << j << " " << aa << "+" << bb << "i " << std::endl;
                }
            }
        }
        ofs_running << std::endl;
        ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
    }

    int info = 0;
    int naroc[2] = {0, 0};

    // Create a Tensor for eii
    assert(nband > 0);
    ct::Tensor eii(ct::DataType::DT_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({nband}));
    eii.zero();

    for (int iprow = 0; iprow < pv->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < pv->dim1; ++ipcol)
        {
            if (iprow == pv->coord[0] && ipcol == pv->coord[1])
            {
                naroc[0] = pv->nrow;
                naroc[1] = pv->ncol;
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, pv->nb, pv->dim1, ipcol);
                    if (igcol >= nband)
                    {
                        continue;
                    }
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, pv->nb, pv->dim0, iprow);
                        if (igrow >= nband)
                        {
                            continue;
                        }
                        if (igcol == igrow)
                        {
                            eii.data<double>()[igcol] = eij.data<std::complex<double>>()[j * naroc[0] + i].real();
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow

    // Perform MPI reduction to compute ekb
    info = MPI_Allreduce(eii.data<double>(), ekb.data<double>(), nband, MPI_DOUBLE, MPI_SUM, pv->comm());
}

template <typename Device>
void compute_ekb_tensor_lapack(const Parallel_Orbitals* pv,
                               const int nband,
                               const int nlocal,
                               const ct::Tensor& Htmp,
                               const ct::Tensor& psi_k,
                               ct::Tensor& ekb,
                               std::ofstream& ofs_running)
{
    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    // Create Tensor objects for temporary data
    ct::Tensor tmp1(ct::DataType::DT_COMPLEX_DOUBLE,
                    ct_device_type,
                    ct::TensorShape({nlocal * nband})); // tmp1 shape: nlocal * nband
    tmp1.zero();

    ct::Tensor eij(ct::DataType::DT_COMPLEX_DOUBLE,
                   ct_device_type,
                   ct::TensorShape({nlocal * nlocal})); // eij shape: nlocal * nlocal
    // Why not use nband * nband ?????
    eij.zero();

    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta = {0.0, 0.0};

    // Perform matrix multiplication: tmp1 = Htmp * psi_k
    ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                              'N',
                                                              nlocal,
                                                              nband,
                                                              nlocal,
                                                              &alpha,
                                                              Htmp.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of Htmp
                                                              psi_k.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of psi_k
                                                              &beta,
                                                              tmp1.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of tmp1

    // Perform matrix multiplication: eij = psi_k^dagger * tmp1
    ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('C',
                                                              'N',
                                                              nband,
                                                              nband,
                                                              nlocal,
                                                              &alpha,
                                                              psi_k.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of psi_k
                                                              tmp1.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of tmp1
                                                              &beta,
                                                              eij.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of eij

    if (PARAM.inp.td_print_eij >= 0.0)
    {
        ct::Tensor eij_cpu = eij.to_device<ct::DEVICE_CPU>();

        ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
        ofs_running << " Eij:" << std::endl;
        for (int i = 0; i < nband; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nband; j++)
            {
                double aa = eij_cpu.data<std::complex<double>>()[in + j].real();
                double bb = eij_cpu.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < PARAM.inp.td_print_eij)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < PARAM.inp.td_print_eij)
                {
                    bb = 0.0;
                }
                if (aa > 0.0 || bb > 0.0)
                {
                    ofs_running << i << " " << j << " " << aa << "+" << bb << "i " << std::endl;
                }
            }
        }
        ofs_running << std::endl;
        ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
    }

    // Extract diagonal elements of eij into ekb
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        // GPU implementation
        for (int i = 0; i < nband; ++i)
        {
            base_device::memory::synchronize_memory_op<double, Device, Device>()(
                ekb.data<double>() + i,
                reinterpret_cast<const double*>(eij.data<std::complex<double>>() + i * nlocal + i),
                1);
        }
    }
    else
    {
        // CPU implementation
        for (int i = 0; i < nband; ++i)
        {
            ekb.data<double>()[i] = eij.data<std::complex<double>>()[i * nlocal + i].real();
        }
    }
}

// Explicit instantiation of template functions
template void compute_ekb_tensor_lapack<base_device::DEVICE_CPU>(const Parallel_Orbitals* pv,
                                                                 const int nband,
                                                                 const int nlocal,
                                                                 const ct::Tensor& Htmp,
                                                                 const ct::Tensor& psi_k,
                                                                 ct::Tensor& ekb,
                                                                 std::ofstream& ofs_running);

#if ((defined __CUDA) /* || (defined __ROCM) */)
template void compute_ekb_tensor_lapack<base_device::DEVICE_GPU>(const Parallel_Orbitals* pv,
                                                                 const int nband,
                                                                 const int nlocal,
                                                                 const ct::Tensor& Htmp,
                                                                 const ct::Tensor& psi_k,
                                                                 ct::Tensor& ekb,
                                                                 std::ofstream& ofs_running);
#endif // __CUDA
#endif // __MPI

} // namespace module_tddft
