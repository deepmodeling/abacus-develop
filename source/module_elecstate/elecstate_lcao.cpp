#include "elecstate_lcao.h"

#include "cal_dm.h"
#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"

namespace elecstate
{
int ElecStateLCAO::out_wfc_lcao = 0;
int ElecStateLCAO::out_wfc_flag = 0;
bool ElecStateLCAO::need_psi_grid = 1;

// multi-k case
void ElecStateLCAO::psiToRho(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    this->calculate_weights();
    this->calEBand();

    ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");

    // this part for calculating dm_k in 2d-block format, not used for charge now
    //    psi::Psi<std::complex<double>> dm_k_2d();

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx"
        || GlobalV::KS_SOLVER == "lapack") // Peize Lin test 2019-05-15
    {
        cal_dm(this->loc->ParaV, this->wg, psi, this->loc->dm_k);
    }

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx" || GlobalV::KS_SOLVER == "lapack")
    {
        for (int ik = 0; ik < psi.get_nk(); ik++)
        {
            psi.fix_k(ik);
            this->print_psi(psi);
        }
    }

    this->loc->cal_dk_k(*this->lowf->gridt, this->wg, (*this->klist));
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
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
void ElecStateLCAO::psiToRho(const psi::Psi<double>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    this->calculate_weights();
    this->calEBand();

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx" || GlobalV::KS_SOLVER == "lapack")
    {
        ModuleBase::timer::tick("ElecStateLCAO", "cal_dm_2d");

        // psi::Psi<double> dm_gamma_2d;
        //  caution:wfc and dm
        cal_dm(this->loc->ParaV, this->wg, psi, this->loc->dm_gamma);

        ModuleBase::timer::tick("ElecStateLCAO", "cal_dm_2d");

        for (int ik = 0; ik < psi.get_nk(); ++ik)
        {
            // for gamma_only case, no convertion occured, just for print.
            if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx")
            {
                psi.fix_k(ik);
                this->print_psi(psi);
            }
            this->loc->cal_dk_gamma_from_2D_pub();
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

void ElecStateLCAO::print_psi(const psi::Psi<double>& psi_in, const int istep)
{
    if (!ElecStateLCAO::out_wfc_lcao)
        return;

    // output but not do "2d-to-grid" conversion
    double** wfc_grid = nullptr;
#ifdef __MPI
    this->lowf->wfc_2d_to_grid(istep, out_wfc_flag, psi_in.get_pointer(), wfc_grid, this->ekb, this->wg);
#endif
    return;
}
void ElecStateLCAO::print_psi(const psi::Psi<std::complex<double>>& psi_in, const int istep)
{
    if (!ElecStateLCAO::out_wfc_lcao && !ElecStateLCAO::need_psi_grid)
        return;

    // output but not do "2d-to-grid" conversion
    std::complex<double>** wfc_grid = nullptr;
    int ik = psi_in.get_current_k();
    if (ElecStateLCAO::need_psi_grid)
    {
        wfc_grid = this->lowf->wfc_k_grid[ik];
    }
#ifdef __MPI
    this->lowf->wfc_2d_to_grid(istep,
                               ElecStateLCAO::out_wfc_flag,
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
    if (ElecStateLCAO::need_psi_grid && GlobalV::NSPIN == 4)
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

} // namespace elecstate