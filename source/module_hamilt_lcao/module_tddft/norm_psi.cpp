#include "norm_psi.h"

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

void norm_psi(const Parallel_Orbitals* pv,
              const int nband,
              const int nlocal,
              const std::complex<double>* Stmp,
              std::complex<double>* psi_k,
              std::ofstream& ofs_running,
              const int print_matrix)
{
    assert(pv->nloc_wfc > 0 && pv->nloc > 0);

    std::complex<double>* tmp1 = new std::complex<double>[pv->nloc_wfc];
    ModuleBase::GlobalFunc::ZEROS(tmp1, pv->nloc_wfc);

    std::complex<double>* Cij = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(Cij, pv->nloc);

    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nlocal,
                             1.0,
                             Stmp,
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
                             Cij,
                             1,
                             1,
                             pv->desc_Eij);

    if (print_matrix)
    {
        ofs_running << "original Cij :" << std::endl;
        for (int i = 0; i < pv->ncol; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->nrow; j++)
            {
                double aa = Cij[in + j].real();
                double bb = Cij[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    int naroc[2] = {0, 0}; // maximum number of row or column

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
                            Cij[j * naroc[0] + i] = {1.0 / sqrt(Cij[j * naroc[0] + i].real()), 0.0};
                        }
                        else
                        {
                            Cij[j * naroc[0] + i] = {0.0, 0.0};
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow

    BlasConnector::copy(pv->nloc_wfc, psi_k, 1, tmp1, 1);

    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nband,
                             1.0,
                             tmp1,
                             1,
                             1,
                             pv->desc_wfc,
                             Cij,
                             1,
                             1,
                             pv->desc_Eij,
                             0.0,
                             psi_k,
                             1,
                             1,
                             pv->desc_wfc);

    if (print_matrix)
    {
        ofs_running << " Cij:" << std::endl;
        for (int i = 0; i < pv->ncol; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->nrow; j++)
            {
                ofs_running << Cij[in + j].real() << "+" << Cij[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa = psi_k[in + j].real();
                double bb = psi_k[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " psi_k before normalization:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa = tmp1[in + j].real();
                double bb = tmp1[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
    }

    delete[] tmp1;
    delete[] Cij;
}

void norm_psi_tensor(const Parallel_Orbitals* pv,
                     const int nband,
                     const int nlocal,
                     const ct::Tensor& Stmp,
                     ct::Tensor& psi_k,
                     std::ofstream& ofs_running,
                     const int print_matrix)
{
    assert(pv->nloc_wfc > 0 && pv->nloc > 0);

    // Create Tensor objects for temporary data
    ct::Tensor tmp1(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({pv->nloc_wfc}));
    tmp1.zero();

    ct::Tensor Cij(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({pv->nloc}));
    Cij.zero();

    // Perform matrix multiplication: tmp1 = Stmp * psi_k
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nlocal,
                             1.0,
                             Stmp.data<std::complex<double>>(),
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

    // Perform matrix multiplication: Cij = psi_k^dagger * tmp1
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
                             Cij.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_Eij);

    if (print_matrix)
    {
        ofs_running << "original Cij :" << std::endl;
        for (int i = 0; i < pv->ncol; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->nrow; j++)
            {
                double aa = Cij.data<std::complex<double>>()[in + j].real();
                double bb = Cij.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    int naroc[2] = {0, 0}; // maximum number of row or column

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
                            Cij.data<std::complex<double>>()[j * naroc[0] + i]
                                = {1.0 / sqrt(Cij.data<std::complex<double>>()[j * naroc[0] + i].real()), 0.0};
                        }
                        else
                        {
                            Cij.data<std::complex<double>>()[j * naroc[0] + i] = {0.0, 0.0};
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow

    // Copy psi_k to tmp1 (using deep copy)
    // tmp1.CopyFrom(psi_k); // Does not work because this will cause tmp1 and psi_k to share the same data
    tmp1 = psi_k; // operator= overload for Tensor class

    // Perform matrix multiplication: psi_k = tmp1 * Cij
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nband,
                             1.0,
                             tmp1.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_wfc,
                             Cij.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_Eij,
                             0.0,
                             psi_k.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_wfc);

    if (print_matrix)
    {
        ofs_running << " Cij:" << std::endl;
        for (int i = 0; i < pv->ncol; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->nrow; j++)
            {
                ofs_running << Cij.data<std::complex<double>>()[in + j].real() << "+"
                            << Cij.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa = psi_k.data<std::complex<double>>()[in + j].real();
                double bb = psi_k.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " psi_k before normalization:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa = tmp1.data<std::complex<double>>()[in + j].real();
                double bb = tmp1.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
    }
}

template <typename Device>
void norm_psi_tensor_lapack(const Parallel_Orbitals* pv,
                            const int nband,
                            const int nlocal,
                            const ct::Tensor& Stmp,
                            ct::Tensor& psi_k,
                            std::ofstream& ofs_running,
                            const int print_matrix)
{
    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    // Create Tensor objects for temporary data
    ct::Tensor tmp1(
        ct::DataType::DT_COMPLEX_DOUBLE,
        ct_device_type,
        ct::TensorShape({nlocal * nband})); // tmp1 shape: nlocal * nband (under 2D block cyclic is pv->nloc_wfc)
    tmp1.zero();

    ct::Tensor Cij(ct::DataType::DT_COMPLEX_DOUBLE,
                   ct_device_type,
                   ct::TensorShape({nlocal * nlocal})); // Cij shape: nlocal * nlocal
    Cij.zero();

    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta = {0.0, 0.0};

    // Perform matrix multiplication: tmp1 = Stmp * psi_k
    ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                              'N',
                                                              nlocal,
                                                              nband,
                                                              nlocal,
                                                              &alpha,
                                                              Stmp.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of Stmp
                                                              psi_k.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of psi_k
                                                              &beta,
                                                              tmp1.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of tmp1

    // Perform matrix multiplication: Cij = psi_k^dagger * tmp1
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
                                                              Cij.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of Cij

    if (print_matrix)
    {
        ct::Tensor Cij_print_cpu = Cij.to_device<ct::DEVICE_CPU>();

        ofs_running << "original Cij :" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                double aa = Cij_print_cpu.data<std::complex<double>>()[in + j].real();
                double bb = Cij_print_cpu.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    // Normalize Cij: set diagonal elements to 1/sqrt(Cij[i][i]), off-diagonal elements to 0
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        // Step 1: Copy Cij from GPU to CPU
        ct::Tensor Cij_cpu = Cij.to_device<ct::DEVICE_CPU>();

        // Step 2: Perform normalization on CPU
        for (int i = 0; i < nband; ++i)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nband; ++j)
            {
                if (i == j)
                {
                    Cij_cpu.data<std::complex<double>>()[in + j]
                        = {1.0 / sqrt(Cij_cpu.data<std::complex<double>>()[in + j].real()), 0.0};
                }
                else
                {
                    Cij_cpu.data<std::complex<double>>()[in + j] = {0.0, 0.0};
                }
            }
        }

        // Step 3: Copy normalized Cij back to GPU
        Cij = Cij_cpu.to_device<ct_Device>();
    }
    else
    {
        // CPU implementation
        for (int i = 0; i < nband; ++i)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nband; ++j)
            {
                if (i == j)
                {
                    Cij.data<std::complex<double>>()[in + j]
                        = {1.0 / sqrt(Cij.data<std::complex<double>>()[in + j].real()), 0.0};
                }
                else
                {
                    Cij.data<std::complex<double>>()[in + j] = {0.0, 0.0};
                }
            }
        }
    }

    // Copy psi_k to tmp1 (using deep copy)
    // tmp1.CopyFrom(psi_k); // Does not work because this will cause tmp1 and psi_k to share the same data
    tmp1 = psi_k; // operator= overload for Tensor class

    // Perform matrix multiplication: psi_k = tmp1 * Cij
    ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                              'N',
                                                              nlocal,
                                                              nband,
                                                              nband,
                                                              &alpha,
                                                              tmp1.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of tmp1
                                                              Cij.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of Cij
                                                              &beta,
                                                              psi_k.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of psi_k

    if (print_matrix)
    {
        ct::Tensor Cij_print_cpu = Cij.to_device<ct::DEVICE_CPU>();
        ct::Tensor psi_k_cpu = psi_k.to_device<ct::DEVICE_CPU>();
        ct::Tensor tmp1_cpu = tmp1.to_device<ct::DEVICE_CPU>();

        ofs_running << " Cij:" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                ofs_running << Cij_print_cpu.data<std::complex<double>>()[in + j].real() << "+"
                            << Cij_print_cpu.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < nband; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                double aa = psi_k_cpu.data<std::complex<double>>()[in + j].real();
                double bb = psi_k_cpu.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " psi_k before normalization:" << std::endl;
        for (int i = 0; i < nband; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                double aa = tmp1_cpu.data<std::complex<double>>()[in + j].real();
                double bb = tmp1_cpu.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
    }
}

// Explicit instantiation of template functions
template void norm_psi_tensor_lapack<base_device::DEVICE_CPU>(const Parallel_Orbitals* pv,
                                                              const int nband,
                                                              const int nlocal,
                                                              const ct::Tensor& Stmp,
                                                              ct::Tensor& psi_k,
                                                              std::ofstream& ofs_running,
                                                              const int print_matrix);
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void norm_psi_tensor_lapack<base_device::DEVICE_GPU>(const Parallel_Orbitals* pv,
                                                              const int nband,
                                                              const int nlocal,
                                                              const ct::Tensor& Stmp,
                                                              ct::Tensor& psi_k,
                                                              std::ofstream& ofs_running,
                                                              const int print_matrix);
#endif // __CUDA
#endif // __MPI
} // namespace module_tddft
