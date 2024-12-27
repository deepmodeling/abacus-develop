#include "bandenergy.h"

#include "evolve_elec.h"
#include "module_base/lapack_connector.h"
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
                 double* ekb)
{

    std::complex<double>* tmp1 = new std::complex<double>[pv->nloc_wfc];
    ModuleBase::GlobalFunc::ZEROS(tmp1, pv->nloc_wfc);

    std::complex<double>* Eij = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(Eij, pv->nloc);

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
                             Eij,
                             1,
                             1,
                             pv->desc_Eij);

    if (Evolve_elec::td_print_eij > 0.0)
    {
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
        GlobalV::ofs_running << " Eij:" << std::endl;
        for (int i = 0; i < pv->nrow_bands; i++)
        {
            for (int j = 0; j < pv->ncol_bands; j++)
            {
                double aa = 0.0, bb = 0.0;
                aa = Eij[i * pv->ncol + j].real();
                bb = Eij[i * pv->ncol + j].imag();
                if (std::abs(aa) < Evolve_elec::td_print_eij)
                    aa = 0.0;
                if (std::abs(bb) < Evolve_elec::td_print_eij)
                    bb = 0.0;
                if (aa > 0.0 || bb > 0.0)
                {
                    GlobalV::ofs_running << i << " " << j << " " << aa << "+" << bb << "i " << std::endl;
                }
            }
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
    }

    int info = 0;
    int naroc[2];

    double* Eii = new double[nband];
    ModuleBase::GlobalFunc::ZEROS(Eii, nband);
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
                        continue;
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, pv->nb, pv->dim0, iprow);
                        if (igrow >= nband)
                            continue;
                        if (igcol == igrow)
                        {
                            Eii[igcol] = Eij[j * naroc[0] + i].real();
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow
    info = MPI_Allreduce(Eii, ekb, nband, MPI_DOUBLE, MPI_SUM, pv->comm());

    delete[] tmp1;
    delete[] Eij;
    delete[] Eii;
}

void compute_ekb_tensor(const Parallel_Orbitals* pv,
                        const int nband,
                        const int nlocal,
                        const container::Tensor& Htmp,
                        const container::Tensor& psi_k,
                        container::Tensor& ekb)
{
    // Create Tensor objects for temporary data
    container::Tensor tmp1(container::DataType::DT_COMPLEX_DOUBLE,
                           container::DeviceType::CpuDevice,
                           container::TensorShape({pv->nloc_wfc}));
    tmp1.zero();

    container::Tensor Eij(container::DataType::DT_COMPLEX_DOUBLE,
                          container::DeviceType::CpuDevice,
                          container::TensorShape({pv->nloc}));
    Eij.zero();

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

    // Perform matrix multiplication: Eij = psi_k^dagger * tmp1
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
                             Eij.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_Eij);

    if (Evolve_elec::td_print_eij >= 0.0)
    {
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
        GlobalV::ofs_running << " Eij:" << std::endl;
        for (int i = 0; i < pv->nrow_bands; i++)
        {
            for (int j = 0; j < pv->ncol_bands; j++)
            {
                double aa = 0.0, bb = 0.0;
                aa = Eij.data<std::complex<double>>()[i * pv->ncol + j].real();
                bb = Eij.data<std::complex<double>>()[i * pv->ncol + j].imag();
                if (std::abs(aa) < Evolve_elec::td_print_eij)
                    aa = 0.0;
                if (std::abs(bb) < Evolve_elec::td_print_eij)
                    bb = 0.0;
                if (aa > 0.0 || bb > 0.0)
                {
                    GlobalV::ofs_running << i << " " << j << " " << aa << "+" << bb << "i " << std::endl;
                }
            }
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
    }

    int info = 0;
    int naroc[2];

    // Create a Tensor for Eii
    container::Tensor Eii(container::DataType::DT_DOUBLE,
                          container::DeviceType::CpuDevice,
                          container::TensorShape({nband}));
    Eii.zero();

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
                        continue;
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, pv->nb, pv->dim0, iprow);
                        if (igrow >= nband)
                            continue;
                        if (igcol == igrow)
                        {
                            Eii.data<double>()[igcol] = Eij.data<std::complex<double>>()[j * naroc[0] + i].real();
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow

    // Perform MPI reduction to compute ekb
    info = MPI_Allreduce(Eii.data<double>(), ekb.data<double>(), nband, MPI_DOUBLE, MPI_SUM, pv->comm());
}

#endif

} // namespace module_tddft
