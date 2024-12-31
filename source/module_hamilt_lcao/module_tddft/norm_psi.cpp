#include "norm_psi.h"

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

void norm_psi(const Parallel_Orbitals* pv,
              const int nband,
              const int nlocal,
              const std::complex<double>* Stmp,
              std::complex<double>* psi_k,
              const int print_matrix)
{
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
        GlobalV::ofs_running << "original Cij :" << std::endl;
        for (int i = 0; i < pv->ncol; i++)
        {
            for (int j = 0; j < pv->nrow; j++)
            {
                double aa, bb;
                aa = Cij[i * pv->ncol + j].real();
                bb = Cij[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    int naroc[2]; // maximum number of row or column

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
        GlobalV::ofs_running << " Cij:" << std::endl;
        for (int i = 0; i < pv->ncol; i++)
        {
            for (int j = 0; j < pv->nrow; j++)
            {
                GlobalV::ofs_running << Cij[i * pv->ncol + j].real() << "+" << Cij[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa, bb;
                aa = psi_k[i * pv->ncol + j].real();
                bb = psi_k[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " psi_k before normalization:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa, bb;
                aa = tmp1[i * pv->ncol + j].real();
                bb = tmp1[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
    }

    delete[] tmp1;
    delete[] Cij;
}

void norm_psi_tensor(const Parallel_Orbitals* pv,
                     const int nband,
                     const int nlocal,
                     const container::Tensor& Stmp,
                     container::Tensor& psi_k,
                     const int print_matrix)
{
    // Create Tensor objects for temporary data
    container::Tensor tmp1(container::DataType::DT_COMPLEX_DOUBLE,
                           container::DeviceType::CpuDevice,
                           container::TensorShape({pv->nloc_wfc}));
    tmp1.zero();

    container::Tensor Cij(container::DataType::DT_COMPLEX_DOUBLE,
                          container::DeviceType::CpuDevice,
                          container::TensorShape({pv->nloc}));
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
        GlobalV::ofs_running << "original Cij :" << std::endl;
        for (int i = 0; i < pv->ncol; i++)
        {
            for (int j = 0; j < pv->nrow; j++)
            {
                double aa, bb;
                aa = Cij.data<std::complex<double>>()[i * pv->ncol + j].real();
                bb = Cij.data<std::complex<double>>()[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    int naroc[2]; // maximum number of row or column

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
        GlobalV::ofs_running << " Cij:" << std::endl;
        for (int i = 0; i < pv->ncol; i++)
        {
            for (int j = 0; j < pv->nrow; j++)
            {
                GlobalV::ofs_running << Cij.data<std::complex<double>>()[i * pv->ncol + j].real() << "+"
                                     << Cij.data<std::complex<double>>()[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa, bb;
                aa = psi_k.data<std::complex<double>>()[i * pv->ncol + j].real();
                bb = psi_k.data<std::complex<double>>()[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " psi_k before normalization:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa, bb;
                aa = tmp1.data<std::complex<double>>()[i * pv->ncol + j].real();
                bb = tmp1.data<std::complex<double>>()[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
    }
}

#endif
} // namespace module_tddft
