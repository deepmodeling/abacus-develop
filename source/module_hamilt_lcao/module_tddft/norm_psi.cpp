#include "norm_psi.h"

#include <complex>
#include <iostream>

#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"

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
                    aa = 0.0;
                if (std::abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    int info;
    int myid;
    MPI_Comm_rank(pv->comm_2D, &myid);
    int naroc[2]; // maximum number of row or column

    for (int iprow = 0; iprow < pv->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < pv->dim1; ++ipcol)
        {
            const int coord[2] = {iprow, ipcol};
            int src_rank;
            info = MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
            if (myid == src_rank)
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
    }     // loop iprow

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
                    aa = 0.0;
                if (std::abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " psi_k before normalization:" << std::endl;
        for (int i = 0; i < pv->ncol; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa, bb;
                aa = tmp1[i * pv->ncol + j].real();
                bb = tmp1[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                    aa = 0.0;
                if (std::abs(bb) < 1e-8)
                    bb = 0.0;
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
#endif
} // namespace module_tddft
