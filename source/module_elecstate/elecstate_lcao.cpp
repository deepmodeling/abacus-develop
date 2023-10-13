#include "elecstate_lcao.h"

#include "cal_dm.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace elecstate
{
template <typename TK>
int ElecStateLCAO<TK>::out_wfc_lcao = 0;

template <typename TK>
int ElecStateLCAO<TK>::out_wfc_flag = 0;

template <typename TK>
bool ElecStateLCAO<TK>::need_psi_grid = 1;

template <>
void ElecStateLCAO<double>::print_psi(const psi::Psi<double>& psi_in, const int istep)
{
    if (!ElecStateLCAO<double>::out_wfc_lcao)
        return;

    // output but not do  "2d-to-grid" conversion
    double** wfc_grid = nullptr;
#ifdef __MPI
    this->lowf->wfc_2d_to_grid(istep, out_wfc_flag, psi_in.get_pointer(), wfc_grid, this->ekb, this->wg);
#endif
    return;
}

template <>
void ElecStateLCAO<std::complex<double>>::print_psi(const psi::Psi<std::complex<double>>& psi_in, const int istep)
{
    if (!ElecStateLCAO<std::complex<double>>::out_wfc_lcao && !ElecStateLCAO<std::complex<double>>::need_psi_grid)
        return;

    // output but not do "2d-to-grid" conversion
    std::complex<double>** wfc_grid = nullptr;
    int ik = psi_in.get_current_k();
    if (ElecStateLCAO<std::complex<double>>::need_psi_grid)
    {
        wfc_grid = this->lowf->wfc_k_grid[ik];
    }
#ifdef __MPI
    this->lowf->wfc_2d_to_grid(istep,
                               ElecStateLCAO<std::complex<double>>::out_wfc_flag,
                               psi_in.get_pointer(),
                               wfc_grid,
                               ik,
                               this->ekb,
                               this->wg,
                               this->klist->kvec_c);
#else
    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
        {
            this->lowf->wfc_k_grid[ik][ib][iw] = psi_in(ib, iw);
        }
    }
#endif

    // added by zhengdy-soc, rearrange the wfc_k_grid from [up,down,up,down...] to [up,up...down,down...],
    if (ElecStateLCAO<std::complex<double>>::need_psi_grid && GlobalV::NSPIN == 4)
    {
        int row = this->lowf->gridt->lgd;
        std::vector<std::complex<double>> tmp(row);
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            for (int iw = 0; iw < row / GlobalV::NPOL; iw++)
            {
                tmp[iw] = this->lowf->wfc_k_grid[ik][ib][iw * GlobalV::NPOL];
                tmp[iw + row / GlobalV::NPOL] = this->lowf->wfc_k_grid[ik][ib][iw * GlobalV::NPOL + 1];
            }
            for (int iw = 0; iw < row; iw++)
            {
                this->lowf->wfc_k_grid[ik][ib][iw] = tmp[iw];
            }
        }
    }

    return;
}

// multi-k case
template <>
void ElecStateLCAO<std::complex<double>>::psiToRho(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    this->calculate_weights();
    this->calEBand();

    ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");

    // this part for calculating DMK in 2d-block format, not used for charge now
    //    psi::Psi<std::complex<double>> dm_k_2d();

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx"
        || GlobalV::KS_SOLVER == "lapack") // Peize Lin test 2019-05-15
    {
        //cal_dm(this->loc->ParaV, this->wg, psi, this->loc->dm_k);
        elecstate::cal_dm_psi(this->DM->get_paraV_pointer(), this->wg, psi, *(this->DM));
        this->DM->cal_DMR();

// interface for RI-related calculation, which needs loc.dm_k
#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            const K_Vectors* kv = this->DM->get_kv_pointer();
            this->loc->dm_k.resize(kv->nks);
            for (int ik = 0; ik < kv->nks; ++ik)
            {
                this->loc->set_dm_k(ik, this->DM->get_DMK_pointer(ik));         
            }
        }
#endif

    }
    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx" || GlobalV::KS_SOLVER == "lapack")
    {
        for (int ik = 0; ik < psi.get_nk(); ik++)
        {
            psi.fix_k(ik);
            this->print_psi(psi);
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
    this->uhm->GK.transfer_DM2DtoGrid(this->DM->get_DMR_vector()); // transfer DM2D to DM_grid in gint
    Gint_inout inout(this->loc->DM_R, this->charge->rho, Gint_Tools::job_type::rho);
    this->uhm->GK.cal_gint(&inout);

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[0], this->charge->nrxx);
        Gint_inout inout1(this->loc->DM_R, this->charge->kin_r, Gint_Tools::job_type::tau);
        this->uhm->GK.cal_gint(&inout1);
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

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx" || GlobalV::KS_SOLVER == "lapack")
    {
        ModuleBase::timer::tick("ElecStateLCAO", "cal_dm_2d");
        // get DMK in 2d-block format
        //cal_dm(this->loc->ParaV, this->wg, psi, this->loc->dm_gamma);
        elecstate::cal_dm_psi(this->DM->get_paraV_pointer(), this->wg, psi, *(this->DM));
        this->DM->cal_DMR();

// interface for RI-related calculation, which needs loc.dm_gamma    
#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx || this->loc->out_dm)
        {
            this->loc->dm_gamma.resize(GlobalV::NSPIN);
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                this->loc->set_dm_gamma(is, this->DM->get_DMK_pointer(is));    
            }
        }
#else
        if (this->loc->out_dm) // keep interface for old Output_DM until new one is ready
        {
            this->loc->dm_gamma.resize(GlobalV::NSPIN);
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                this->loc->set_dm_gamma(is, this->DM->get_DMK_pointer(is));    
            }
        }
#endif

        ModuleBase::timer::tick("ElecStateLCAO", "cal_dm_2d");

        for (int ik = 0; ik < psi.get_nk(); ++ik)
        {
            // for gamma_only case, no convertion occured, just for print.
            if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx")
            {
                psi.fix_k(ik);
                this->print_psi(psi);
            }
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
    this->uhm->GG.transfer_DM2DtoGrid(this->DM->get_DMR_vector()); // transfer DM2D to DM_grid in gint
    Gint_inout inout(this->loc->DM, this->charge->rho, Gint_Tools::job_type::rho);
    this->uhm->GG.cal_gint(&inout);
    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[0], this->charge->nrxx);
        }
        Gint_inout inout1(this->loc->DM, this->charge->kin_r, Gint_Tools::job_type::tau);
        this->uhm->GG.cal_gint(&inout1);
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
    SpinConstrain<double, psi::DEVICE_CPU>& sc = SpinConstrain<double>::getScInstance();
    return sc.cal_escon();
}

template<>
double ElecStateLCAO<std::complex<double>>::get_spin_constrain_energy()
{
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc = SpinConstrain<std::complex<double>>::getScInstance();
    return sc.cal_escon();
}

template class ElecStateLCAO<double>; // Gamma_only case
template class ElecStateLCAO<std::complex<double>>; // multi-k case

} // namespace elecstate