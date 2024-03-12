#include "bandenergy.h"

#include <complex>
#include <iostream>

#include "evolve_elec.h"
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
                double aa, bb;
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

    int info;
    int myid;
    int naroc[2];
    MPI_Comm_rank(pv->comm_2D, &myid);

    double* Eii = new double[nband];
    ModuleBase::GlobalFunc::ZEROS(Eii, nband);
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
                            Eii[igcol] = Eij[j * naroc[0] + i].real();
                        }
                    }
                }
            }
        } // loop ipcol
    }     // loop iprow
    info = MPI_Allreduce(Eii, ekb, nband, MPI_DOUBLE, MPI_SUM, pv->comm_2D);

    delete[] tmp1;
    delete[] Eij;
    delete[] Eii;
}

#endif

} // namespace module_tddft