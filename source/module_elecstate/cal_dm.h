#ifndef CAL_DM_H
#define CAL_DM_H

// #include "math_tools.h"
#include "module_elecstate/math_tools.h"  // modyfied by JG for testing
#include "module_base/timer.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"

namespace elecstate
{

// for gamma_only(double case) and multi-k(complex<double> case)
inline void cal_dm(const Parallel_Orbitals* ParaV, const ModuleBase::matrix& wg, const psi::Psi<double>& wfc, std::vector<ModuleBase::matrix>& dm)
{
    ModuleBase::TITLE("elecstate", "cal_dm");
    ModuleBase::timer::tick("elecstate","cal_dm");

    //dm.resize(wfc.get_nk(), ParaV->ncol, ParaV->nrow);
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    // dm = wfc.T * wg * wfc.conj()
    // dm[is](iw1,iw2) = \sum_{ib} wfc[is](ib,iw1).T * wg(is,ib) * wfc[is](ib,iw2).conj()
    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        wfc.fix_k(ik);
        //dm.fix_k(ik);
        dm[ik].create(ParaV->ncol, ParaV->nrow);
        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        psi::Psi<double> wg_wfc(wfc, 1);

        int ib_global = 0;
        for (int ib_local = 0; ib_local < nbands_local; ++ib_local)
        {
            while (ib_local != ParaV->global2local_col(ib_global))
            {
                ++ib_global;
                if (ib_global >= wg.nc)
                {
                    break;
                    ModuleBase::WARNING_QUIT("ElecStateLCAO::cal_dm", "please check global2local_col!");
                }
            }
            if (ib_global >= wg.nc) continue;
            const double wg_local = wg(ik, ib_global);
            double* wg_wfc_pointer = &(wg_wfc(0, ib_local, 0));
            BlasConnector::scal(nbasis_local, wg_local, wg_wfc_pointer, 1);
        }

        // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
#ifdef __MPI
        psiMulPsiMpi(wg_wfc, wfc, dm[ik], ParaV->desc_wfc, ParaV->desc);
#else
        psiMulPsi(wg_wfc, wfc, dm[ik]);
#endif
    }
    ModuleBase::timer::tick("elecstate","cal_dm");

    return;
}

inline void cal_dm(const Parallel_Orbitals* ParaV, const ModuleBase::matrix& wg, const psi::Psi<std::complex<double>>& wfc, std::vector<ModuleBase::ComplexMatrix>& dm)
{
    ModuleBase::TITLE("elecstate", "cal_dm");
    ModuleBase::timer::tick("elecstate","cal_dm");

    //dm.resize(wfc.get_nk(), ParaV->ncol, ParaV->nrow);
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    // dm = wfc.T * wg * wfc.conj()
    // dm[is](iw1,iw2) = \sum_{ib} wfc[is](ib,iw1).T * wg(is,ib) * wfc[is](ib,iw2).conj()
    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        wfc.fix_k(ik);
        //dm.fix_k(ik);
        dm[ik].create(ParaV->ncol, ParaV->nrow);
        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        psi::Psi<std::complex<double>> wg_wfc(1, wfc.get_nbands(), wfc.get_nbasis(), nullptr);
        const std::complex<double>* pwfc = wfc.get_pointer();
        std::complex<double>* pwg_wfc = wg_wfc.get_pointer();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int i = 0;i<wg_wfc.size();++i)
        {
            pwg_wfc[i] = conj(pwfc[i]);
        }

        int ib_global = 0;
        for (int ib_local = 0; ib_local < nbands_local; ++ib_local)
        {
            while (ib_local != ParaV->global2local_col(ib_global))
            {
                ++ib_global;
                if (ib_global >= wg.nc)
                {
                    break;
                    ModuleBase::WARNING_QUIT("ElecStateLCAO::cal_dm", "please check global2local_col!");
                }
            }
            if (ib_global >= wg.nc) continue;
            const double wg_local = wg(ik, ib_global);
            std::complex<double>* wg_wfc_pointer = &(wg_wfc(0, ib_local, 0));
            BlasConnector::scal(nbasis_local, wg_local, wg_wfc_pointer, 1);
        }

        // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
#ifdef __MPI
        psiMulPsiMpi(wg_wfc, wfc, dm[ik], ParaV->desc_wfc, ParaV->desc);
#else
        psiMulPsi(wg_wfc, wfc, dm[ik]);
#endif
    }

    ModuleBase::timer::tick("elecstate","cal_dm");
    return;
}

}//namespace elecstate

#endif
