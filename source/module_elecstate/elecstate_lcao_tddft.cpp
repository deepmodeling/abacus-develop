#include "elecstate_lcao_tddft.h"

#include "cal_dm.h"
#include "module_base/timer.h"
#include "module_gint/grid_technique.h"
#include "src_pw/global.h"

namespace elecstate
{

// multi-k case
void ElecStateLCAO_TDDFT::psiToRho_td(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    this->calculate_weights_td();
    this->calEBand();

    ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");

    // this part for calculating dm_k in 2d-block format, not used for charge now
    //    psi::Psi<std::complex<double>> dm_k_2d();

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx"
        || GlobalV::KS_SOLVER == "lapack") // Peize Lin test 2019-05-15
    {
        cal_dm(this->loc->ParaV, this->wg, psi, this->loc->dm_k);
    }

    // this part for steps:
    // 1. psi_k transform from 2d-block to grid format
    // 2. psi_k_grid -> DM_R
    // 3. DM_R -> rho(r)
    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx" || GlobalV::KS_SOLVER == "lapack")
    {
        for (int ik = 0; ik < psi.get_nk(); ik++)
        {
            psi.fix_k(ik);
            this->print_psi(psi);
        }
    }

    this->loc->cal_dk_k(GlobalC::GridT, this->wg);
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    Gint_inout inout(this->loc->DM_R, this->charge, Gint_Tools::job_type::rho);
    this->uhm->GK.cal_gint(&inout);

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");
    return;
}

void ElecStateLCAO_TDDFT::calculate_weights_td()
{
    ModuleBase::TITLE("ElecState", "calculate_weights");

    if (GlobalV::ocp == 1)
    {
        for (int ik = 0; ik < GlobalC::kv.nks; ik++)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                this->wg(ik, ib) = GlobalV::ocp_kb[ik * GlobalV::NBANDS + ib];
            }
        }
    }
    return;
}

} // namespace elecstate