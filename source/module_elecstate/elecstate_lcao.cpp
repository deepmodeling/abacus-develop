#include "elecstate_lcao.h"

#include "cal_dm.h"
#include "module_base/timer.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"

#include <vector>

namespace elecstate
{

// multi-k case
template <>
void ElecStateLCAO<std::complex<double>>::psiToRho(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    // // the calculations of dm, and dm -> rho are, technically, two separate
    // // functionalities, as we cannot rule out the possibility that we may have a
    // // dm from other sources, such as read from file. However, since we are not
    // // separating them now, I opt to add a flag to control how dm is obtained as
    // // of now
    // if (!PARAM.inp.dm_to_rho)
    // {
    //     ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");

    //     // this part for calculating DMK in 2d-block format, not used for charge
    //     // now
    //     //    psi::Psi<std::complex<double>> dm_k_2d();

    //     if (PARAM.inp.ks_solver == "genelpa" || PARAM.inp.ks_solver == "elpa" || PARAM.inp.ks_solver ==
    //     "scalapack_gvx" || PARAM.inp.ks_solver == "lapack"
    //         || PARAM.inp.ks_solver == "cusolver" || PARAM.inp.ks_solver == "cusolvermp"
    //         || PARAM.inp.ks_solver == "cg_in_lcao") // Peize Lin test 2019-05-15
    //     {
    //         elecstate::cal_dm_psi(this->DM->get_paraV_pointer(),
    //                               this->wg,
    //                               psi,
    //                               *(this->DM));
    //         this->DM->cal_DMR();
    //     }
    // }

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is],
                                      this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    this->gint_k->transfer_DM2DtoGrid(this->DM->get_DMR_vector()); // transfer DM2D to DM_grid in gint
    Gint_inout inout(this->charge->rho, Gint_Tools::job_type::rho, PARAM.inp.nspin);
    this->gint_k->cal_gint(&inout);

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        this->cal_tau(psi);
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

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is],
                                      this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------
    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");

    this->gint_gamma->transfer_DM2DtoGrid(this->DM->get_DMR_vector()); // transfer DM2D to DM_grid in gint

    Gint_inout inout(this->charge->rho, Gint_Tools::job_type::rho, PARAM.inp.nspin);

    this->gint_gamma->cal_gint(&inout);

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        this->cal_tau(psi);
    }

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");
    return;
}

template <typename TK>
void ElecStateLCAO<TK>::init_DM(const K_Vectors* kv, const Parallel_Orbitals* paraV, const int nspin)
{
    const int nspin_dm = nspin == 2 ? 2 : 1;
    this->DM = new DensityMatrix<TK, double>(paraV, nspin_dm, kv->kvec_d, kv->get_nks() / nspin_dm);
}

template <>
double ElecStateLCAO<double>::get_spin_constrain_energy()
{
    spinconstrain::SpinConstrain<double>& sc = spinconstrain::SpinConstrain<double>::getScInstance();
    return sc.cal_escon();
}

template <>
double ElecStateLCAO<std::complex<double>>::get_spin_constrain_energy()
{
    spinconstrain::SpinConstrain<std::complex<double>>& sc
        = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
    return sc.cal_escon();
}

#ifdef __PEXSI
template <>
void ElecStateLCAO<double>::dmToRho(std::vector<double*> pexsi_DM, std::vector<double*> pexsi_EDM)
{
    ModuleBase::timer::tick("ElecStateLCAO", "dmToRho");

    int nspin = PARAM.inp.nspin;
    if (PARAM.inp.nspin == 4)
    {
        nspin = 1;
    }

    this->get_DM()->pexsi_EDM = pexsi_EDM;

    for (int is = 0; is < nspin; is++)
    {
        this->DM->set_DMK_pointer(is, pexsi_DM[is]);
    }
    DM->cal_DMR();

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is],
                                      this->charge->nrxx); // mohan 2009-11-10
    }

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    this->gint_gamma->transfer_DM2DtoGrid(this->DM->get_DMR_vector()); // transfer DM2D to DM_grid in gint
    Gint_inout inout(this->charge->rho, Gint_Tools::job_type::rho, PARAM.inp.nspin);
    this->gint_gamma->cal_gint(&inout);
    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[0], this->charge->nrxx);
        }
        Gint_inout inout1(this->charge->kin_r, Gint_Tools::job_type::tau);
        this->gint_gamma->cal_gint(&inout1);
    }

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "dmToRho");
    return;
}

template <>
void ElecStateLCAO<std::complex<double>>::dmToRho(std::vector<std::complex<double>*> pexsi_DM,
                                                  std::vector<std::complex<double>*> pexsi_EDM)
{
    ModuleBase::WARNING_QUIT("ElecStateLCAO", "pexsi is not completed for multi-k case");
}

#endif

template class ElecStateLCAO<double>;               // Gamma_only case
template class ElecStateLCAO<std::complex<double>>; // multi-k case

} // namespace elecstate
