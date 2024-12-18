#include "elecstate.h"
#include "module_parameter/parameter.h"
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/tool_title.h"
#include "occupy.h"

namespace elecstate
{

const double* ElecState::getRho(int spin) const
{
    return &(this->charge->rho[spin][0]);
}

void ElecState::fixed_weights(const std::vector<double>& ocp_kb, const int& nbands, const double& nelec)
{
    assert(nbands > 0);
    assert(nelec > 0.0);

    const double ne_thr = 1.0e-5;

    const int num = this->klist->get_nks() * nbands;
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

    for (int ik = 0; ik < this->wg.nr; ++ik)
    {
        for (int ib = 0; ib < this->wg.nc; ++ib)
        {
            this->wg(ik, ib) = ocp_kb[ik * this->wg.nc + ib];
        }
    }
    this->skip_weights = true;

    return;
}


void ElecState::init_nelec_spin()
{
    this->nelec_spin.resize(PARAM.inp.nspin);
    if (PARAM.inp.nspin == 2)
    {
        this->nelec_spin[0] = (PARAM.inp.nelec + PARAM.inp.nupdown) / 2.0;
        this->nelec_spin[1] = (PARAM.inp.nelec - PARAM.inp.nupdown) / 2.0;
    }
}


void ElecState::calculate_weights()
{
    ModuleBase::TITLE("ElecState", "calculate_weights");
    if (this->skip_weights)
    {
        return;
    }

    const int nbands = this->ekb.nc;
    const int nks = this->ekb.nr;

    if (!Occupy::use_gaussian_broadening && !Occupy::fixed_occupations)
    {
        if (PARAM.globalv.two_fermi)
        {
            Occupy::iweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[0],
                             this->ekb,
                             this->eferm.ef_up,
                             this->wg,
                             0,
                             this->klist->isk);
            Occupy::iweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[1],
                             this->ekb,
                             this->eferm.ef_dw,
                             this->wg,
                             1,
                             this->klist->isk);
            // ef = ( ef_up + ef_dw ) / 2.0_dp need??? mohan add 2012-04-16
        }
        else
        {
            // -1 means don't need to consider spin.
            Occupy::iweights(nks,
                             this->klist->wk,
                             nbands,
                             PARAM.inp.nelec,
                             this->ekb,
                             this->eferm.ef,
                             this->wg,
                             -1,
                             this->klist->isk);
        }
    }
    else if (Occupy::use_gaussian_broadening)
    {
        if (PARAM.globalv.two_fermi)
        {
            double demet_up = 0.0;
            double demet_dw = 0.0;
            Occupy::gweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[0],
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             this->ekb,
                             this->eferm.ef_up,
                             demet_up,
                             this->wg,
                             0,
                             this->klist->isk);
            Occupy::gweights(nks,
                             this->klist->wk,
                             nbands,
                             this->nelec_spin[1],
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             this->ekb,
                             this->eferm.ef_dw,
                             demet_dw,
                             this->wg,
                             1,
                             this->klist->isk);
            this->f_en.demet = demet_up + demet_dw;
        }
        else
        {
            // -1 means is no related to spin.
            Occupy::gweights(nks,
                             this->klist->wk,
                             nbands,
                             PARAM.inp.nelec,
                             Occupy::gaussian_parameter,
                             Occupy::gaussian_type,
                             this->ekb,
                             this->eferm.ef,
                             this->f_en.demet,
                             this->wg,
                             -1,
                             this->klist->isk);
        }
#ifdef __MPI
        // qianrui fix a bug on 2021-7-21
        Parallel_Reduce::reduce_double_allpool(GlobalV::KPAR, GlobalV::NPROC_IN_POOL, this->f_en.demet);
#endif
    }
    else if (Occupy::fixed_occupations)
    {
        ModuleBase::WARNING_QUIT("calculate_weights", "other occupations, not implemented");
    }

    return;
}


void ElecState::calEBand()
{
    ModuleBase::TITLE("ElecState", "calEBand");
    // calculate ebands using wg and ekb
    double eband = 0.0;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(+ : eband)
#endif
    for (int ik = 0; ik < this->ekb.nr; ++ik)
    {
        for (int ibnd = 0; ibnd < this->ekb.nc; ibnd++)
        {
            eband += this->ekb(ik, ibnd) * this->wg(ik, ibnd);
        }
    }
    this->f_en.eband = eband;
    if (GlobalV::KPAR != 1 && PARAM.inp.esolver_type != "sdft")
    {
        //==================================
        // Reduce all the Energy in each cpu
        //==================================
        this->f_en.eband /= GlobalV::NPROC_IN_POOL;
#ifdef __MPI
        Parallel_Reduce::reduce_all(this->f_en.eband);
#endif
    }
    return;
}


void ElecState::init_scf(const int istep, 
                         const UnitCell& ucell,
                         const Parallel_Grid& pgrid,
                         const ModuleBase::ComplexMatrix& strucfac, 
                         const bool* numeric,
                         ModuleSymmetry::Symmetry& symm, 
                         const void* wfcpw)
{
    //! core correction potential.
    if (!PARAM.inp.use_paw)
    {
        this->charge->set_rho_core(ucell,strucfac, numeric);
    }
    else
    {
        this->charge->set_rho_core_paw();
    }

    //! other effective potentials need charge density,
    // choose charge density from ionic step 0.
    if (istep == 0)
    {
        this->charge->init_rho(this->eferm,ucell, pgrid, strucfac, symm, (const void*)this->klist, wfcpw);
        this->charge->check_rho(); // check the rho
    }

    //! renormalize the charge density
    this->charge->renormalize_rho();

    //! initialize the potential
    this->pot->init_pot(istep, this->charge);
}


void ElecState::init_ks(Charge* chg_in, // pointer for class Charge
                        const K_Vectors* klist_in,
                        int nk_in,
                        ModulePW::PW_Basis* rhopw_in,
                        const ModulePW::PW_Basis_Big* bigpw_in)
{
    this->charge = chg_in;
    this->klist = klist_in;
    this->charge->set_rhopw(rhopw_in);
    this->bigpw = bigpw_in;
    // init nelec_spin with nelec and nupdown
    this->init_nelec_spin();
    // initialize ekb and wg
    this->ekb.create(nk_in, PARAM.inp.nbands);
    this->wg.create(nk_in, PARAM.inp.nbands);
}

} // namespace elecstate
