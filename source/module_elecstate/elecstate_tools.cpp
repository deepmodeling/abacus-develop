#include "elecstate_tools.h"
#include "occupy.h"
namespace elecstate
{
    void calEBand(const ModuleBase::matrix& ekb,const ModuleBase::matrix& wg,fenergy& f_en)
    {
        ModuleBase::TITLE("ElecState", "calEBand");
        // calculate ebands using wg and ekb
        double eband = 0.0;
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2) reduction(+ : eband)
    #endif
        for (int ik = 0; ik < ekb.nr; ++ik)
        {
            for (int ibnd = 0; ibnd < ekb.nc; ibnd++)
            {
                eband += ekb(ik, ibnd) * wg(ik, ibnd);
            }
        }
        f_en.eband = eband;

    #ifdef __MPI
        const int npool = GlobalV::KPAR * PARAM.inp.bndpar;
        Parallel_Reduce::reduce_double_allpool(npool, GlobalV::NPROC_IN_POOL, f_en.eband);
    #endif
        return;
    }

    void calculate_weights(const ModuleBase::matrix& ekb,
                           ModuleBase::matrix& wg,
                           const K_Vectors* klist,
                           efermi& eferm,
                           fenergy& f_en,
                           std::vector<double>& nelec_spin,
                           const bool skip_weights=false)
    {
        ModuleBase::TITLE("ElecState", "calculate_weights");
        if (skip_weights==true) return;
        
        const int nbands = ekb.nc;
        const int nks = ekb.nr;
        if (!(Occupy::use_gaussian_broadening || Occupy::fixed_occupations))
        {
            if (PARAM.globalv.two_fermi)
            {
                Occupy::iweights(nks,
                                klist->wk,
                                nbands,
                                nelec_spin[0],
                                ekb,
                                eferm.ef_up,
                                wg,
                                0,
                                klist->isk);
                Occupy::iweights(nks,
                                klist->wk,
                                nbands,
                                nelec_spin[1],
                                ekb,
                                eferm.ef_dw,
                                wg,
                                1,
                                klist->isk);
                // ef = ( ef_up + ef_dw ) / 2.0_dp need??? mohan add 2012-04-16
            }
            else
            {
                // -1 means don't need to consider spin.
                Occupy::iweights(nks,
                                klist->wk,
                                nbands,
                                PARAM.inp.nelec,
                                ekb,
                                eferm.ef,
                                wg,
                                -1,
                                klist->isk);
            }
        }
        else if (Occupy::use_gaussian_broadening)
        {
            if (PARAM.globalv.two_fermi)
            {
                double demet_up = 0.0;
                double demet_dw = 0.0;
                Occupy::gweights(nks,
                                klist->wk,
                                nbands,
                                nelec_spin[0],
                                Occupy::gaussian_parameter,
                                Occupy::gaussian_type,
                                ekb,
                                eferm.ef_up,
                                demet_up,
                                wg,
                                0,
                                klist->isk);
                Occupy::gweights(nks,
                                klist->wk,
                                nbands,
                                nelec_spin[1],
                                Occupy::gaussian_parameter,
                                Occupy::gaussian_type,
                                ekb,
                                eferm.ef_dw,
                                demet_dw,
                                wg,
                                1,
                                klist->isk);
                f_en.demet = demet_up + demet_dw;
            }
            else
            {
                // -1 means is no related to spin.
                Occupy::gweights(nks,
                                klist->wk,
                                nbands,
                                PARAM.inp.nelec,
                                Occupy::gaussian_parameter,
                                Occupy::gaussian_type,
                                ekb,
                                eferm.ef,
                                f_en.demet,
                                wg,
                                -1,
                                klist->isk);
            }
    #ifdef __MPI
            const int npool = GlobalV::KPAR * PARAM.inp.bndpar;
            Parallel_Reduce::reduce_double_allpool(npool, GlobalV::NPROC_IN_POOL, f_en.demet);
    #endif
        }
        else if (Occupy::fixed_occupations)
        {
            ModuleBase::WARNING_QUIT("calculate_weights", "other occupations, not implemented");
        }

        return;
    }

    void fixed_weights(const std::vector<double>& ocp_kb,
                       const int& nbands,
                       const double& nelec,
                       const K_Vectors* klist,
                       ModuleBase::matrix& wg,
                       bool& skip_weights)
    {
        assert(nbands > 0);
        assert(nelec > 0.0);

        const double ne_thr = 1.0e-5;

        const int num = klist->get_nks() * nbands;
        if (num != ocp_kb.size())
        {
            ModuleBase::WARNING_QUIT("ElecState::fixed_weights",
                                    "size of occupation array is wrong , please check ocp_set");
        }

        double num_elec = 0.0;
        for (int i = 0; i < ocp_kb.size(); ++i)
        {
            num_elec += ocp_kb[i];
        }

        if (std::abs(num_elec - nelec) > ne_thr)
        {
            ModuleBase::WARNING_QUIT("ElecState::fixed_weights",
                                    "total number of occupations is wrong , please check ocp_set");
        }

        for (int ik = 0; ik < wg.nr; ++ik)
        {
            for (int ib = 0; ib < wg.nc; ++ib)
            {
                wg(ik, ib) = ocp_kb[ik * wg.nc + ib];
            }
        }
        skip_weights = true;

    }
}