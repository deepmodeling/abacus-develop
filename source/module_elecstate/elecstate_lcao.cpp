#include "elecstate_lcao.h"
#include <vector>

#include "cal_dm.h"
#include "module_base/timer.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace elecstate
{

// multi-k case
template <>
void ElecStateLCAO<std::complex<double>>::psiToRho(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    this->calculate_weights();

// the calculations of dm, and dm -> rho are, technically, two separate functionalities, as we cannot
// rule out the possibility that we may have a dm from other sources, such as read from file.
// However, since we are not separating them now, I opt to add a flag to control how dm is obtained as of now
if(!GlobalV::dm_to_rho)
{
    this->calEBand();

    ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");

    // this part for calculating DMK in 2d-block format, not used for charge now
    //    psi::Psi<std::complex<double>> dm_k_2d();

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx" || GlobalV::KS_SOLVER == "lapack"
        || GlobalV::KS_SOLVER == "cusolver" || GlobalV::KS_SOLVER == "cg_in_lcao") // Peize Lin test 2019-05-15
    {
        //cal_dm(this->loc->ParaV, this->wg, psi, this->loc->dm_k);
        elecstate::cal_dm_psi(this->DM->get_paraV_pointer(), this->wg, psi, *(this->DM));
        this->DM->cal_DMR();

// interface for RI-related calculation, which needs loc.dm_k
#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            const K_Vectors* kv = this->DM->get_kv_pointer();
            this->loc->dm_k.resize(kv->get_nks());
            for (int ik = 0; ik < kv->get_nks(); ++ik)
            {
                this->loc->set_dm_k(ik, this->DM->get_DMK_pointer(ik));         
            }
        }
#endif

    }
}
    // old 2D-to-Grid conversion has been replaced by new Gint Refactor 2023/09/25
    //this->loc->cal_dk_k(*this->lowf->gridt, this->wg, (*this->klist));
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    this->gint_k->transfer_DM2DtoGrid(this->DM->get_DMR_vector()); // transfer DM2D to DM_grid in gint
    Gint_inout inout(this->loc->DM_R, this->charge->rho, Gint_Tools::job_type::rho);
    this->gint_k->cal_gint(&inout);

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[is], this->charge->nrxx);
        }
        Gint_inout inout1(this->loc->DM_R, this->charge->kin_r, Gint_Tools::job_type::tau);
        this->gint_k->cal_gint(&inout1);
    }

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");
    return;
}

// Gamma_only case
template <>
void ElecStateLCAO<double>::psiToRho(const psi::Psi<double>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    this->calculate_weights();
    this->calEBand();

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx" || GlobalV::KS_SOLVER == "lapack"
        || GlobalV::KS_SOLVER == "cusolver" || GlobalV::KS_SOLVER == "cg_in_lcao")
    {
        ModuleBase::timer::tick("ElecStateLCAO", "cal_dm_2d");

        // get DMK in 2d-block format
        //cal_dm(this->loc->ParaV, this->wg, psi, this->loc->dm_gamma);
        elecstate::cal_dm_psi(this->DM->get_paraV_pointer(), this->wg, psi, *(this->DM));
        this->DM->cal_DMR();
        if (this->loc->out_dm) // keep interface for old Output_DM until new one is ready
        {
            this->loc->dm_gamma.resize(GlobalV::NSPIN);
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                this->loc->set_dm_gamma(is, this->DM->get_DMK_pointer(is));
            }
        }
        ModuleBase::timer::tick("ElecStateLCAO", "cal_dm_2d");
        for (int ik = 0; ik < psi.get_nk(); ++ik)
        {
            // for gamma_only case, no convertion occured, just for print.
            // old 2D-to-Grid conversion has been replaced by new Gint Refactor 2023/09/25
            if (this->loc->out_dm) // keep interface for old Output_DM until new one is ready
            {
                this->loc->cal_dk_gamma_from_2D_pub();
            }
        }
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------
    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");

    this->gint_gamma->transfer_DM2DtoGrid(this->DM->get_DMR_vector()); // transfer DM2D to DM_grid in gint

    Gint_inout inout(this->loc->DM, this->charge->rho, Gint_Tools::job_type::rho);

    this->gint_gamma->cal_gint(&inout);

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[is], this->charge->nrxx);
        }
        Gint_inout inout1(this->loc->DM, this->charge->kin_r, Gint_Tools::job_type::tau);
        this->gint_gamma->cal_gint(&inout1);
    }

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");
    return;
}

template <typename TK>
void ElecStateLCAO<TK>::init_DM(const K_Vectors* kv, const Parallel_Orbitals* paraV, const int nspin)
{
    this->DM = new DensityMatrix<TK,double>(kv, paraV, nspin);
}

template<>
double ElecStateLCAO<double>::get_spin_constrain_energy()
{
    SpinConstrain<double, base_device::DEVICE_CPU>& sc = SpinConstrain<double>::getScInstance();
    return sc.cal_escon();
}

template<>
double ElecStateLCAO<std::complex<double>>::get_spin_constrain_energy()
{
    SpinConstrain<std::complex<double>, base_device::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>>::getScInstance();
    return sc.cal_escon();
}

#ifdef __PEXSI
template<>
void ElecStateLCAO<double>::dmToRho(std::vector<double*> pexsi_DM, std::vector<double*> pexsi_EDM)
{
    ModuleBase::timer::tick("ElecStateLCAO", "dmToRho");

    int nspin = GlobalV::NSPIN;
    if (GlobalV::NSPIN == 4)
    {
        nspin = 1;
    }

    // old 2D-to-Grid conversion has been replaced by new Gint Refactor 2023/09/25
    if (this->loc->out_dm) // keep interface for old Output_DM until new one is ready
    {
        for (int is = 0; is < nspin; ++is)
        {
            this->loc->set_dm_gamma(is, pexsi_DM[is]);
        }
        this->loc->cal_dk_gamma_from_2D_pub();
    }

    this->get_DM()->pexsi_EDM = pexsi_EDM;
    
    for (int is = 0; is < nspin; is++)
    {
        this->DM->set_DMK_pointer(is, pexsi_DM[is]);
    }
    DM->cal_DMR();
    
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx); // mohan 2009-11-10
    }

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    this->gint_gamma->transfer_DM2DtoGrid(this->DM->get_DMR_vector()); // transfer DM2D to DM_grid in gint
    Gint_inout inout(this->loc->DM, this->charge->rho, Gint_Tools::job_type::rho);
    this->gint_gamma->cal_gint(&inout);
    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[0], this->charge->nrxx);
        }
        Gint_inout inout1(this->loc->DM, this->charge->kin_r, Gint_Tools::job_type::tau);
        this->gint_gamma->cal_gint(&inout1);
    }

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "dmToRho");
    return;
}

template<>
void ElecStateLCAO<std::complex<double>>::dmToRho(std::vector<std::complex<double>*> pexsi_DM, std::vector<std::complex<double>*> pexsi_EDM)
{
    ModuleBase::WARNING_QUIT("ElecStateLCAO", "pexsi is not completed for multi-k case");
}

#endif

template class ElecStateLCAO<double>; // Gamma_only case
template class ElecStateLCAO<std::complex<double>>; // multi-k case

} // namespace elecstate
