//==========================================================
// Author: Jingang Han
// DATE : 2024-10-30
//==========================================================

#include "rdmft.h"
#include "module_rdmft/rdmft_tools.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/elecstate_lcao_cal_tau.h"

namespace rdmft
{



template <typename TK, typename TR>
void RDMFT<TK, TR>::update_ion(UnitCell& ucell_in, 
                               ModulePW::PW_Basis& rho_basis_in,
                               ModuleBase::matrix& vloc_in, 
                               ModuleBase::ComplexMatrix& sf_in)
{
    ucell = &ucell_in;
    rho_basis = &rho_basis_in;
    vloc = &vloc_in;
    sf = &sf_in;

    HR_TV->set_zero();
    this->cal_V_TV();
#ifdef __EXX
    if( GlobalC::exx_info.info_global.cal_exx )
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            Vxc_fromRI_d->cal_exx_ions(ucell_in);
        }
        else
        {
            Vxc_fromRI_c->cal_exx_ions(ucell_in);
        }
    }
#endif

    std::cout << "\n\n\n******\ndo rdmft_esolver.update_ion() successfully\n******\n\n\n" << std::endl;
}


template <typename TK, typename TR>
void RDMFT<TK, TR>::update_elec(UnitCell& ucell,
                                const ModuleBase::matrix& occ_number_in, 
                                const psi::Psi<TK>& wfc_in, const Charge* charge_in)
{
    // update occ_number, wg, wk_fun_occNum
    occ_number = (occ_number_in);
    wg = (occ_number);

    for(int ik=0; ik < wg.nr; ++ik)
    {
        for(int inb=0; inb < wg.nc; ++inb)
        {
            wg(ik, inb) *= kv->wk[ik];
            wk_fun_occNum(ik, inb) = kv->wk[ik] * occNum_func(occ_number(ik, inb), 2, XC_func_rdmft, alpha_power);
        }
    }

    // update wfc
    TK* pwfc_in = &wfc_in(0, 0, 0);
    TK* pwfc = &wfc(0, 0, 0);
    for(int i=0; i<wfc.size(); ++i) { pwfc[i] = pwfc_in[i];
}

    // update charge
    this->update_charge(ucell);

    // "default" = "pbe"
    // if(  !only_exx_type || this->cal_E_type != 1 )
    if( this->cal_E_type != 1 )
    {
        // the second cal_E_type need the complete pot to get effctive_V to calEband and so on.
        this->pelec->pot->update_from_charge(charge, &ucell);
    }

    this->cal_V_hartree();
    this->cal_V_XC(ucell);
    // this->cal_Hk_Hpsi();

    std::cout << "\n******\n" << "update elec in rdmft successfully" << "\n******\n" << std::endl;
}


// this code is copying from function ElecStateLCAO<TK>::psiToRho(), in elecstate_lcao.cpp
template <typename TK, typename TR>
void RDMFT<TK, TR>::update_charge(UnitCell& ucell)
{
    if( PARAM.inp.gamma_only )
    {
        // calculate DMK and DMR
        elecstate::DensityMatrix<TK, double> DM_gamma_only(ParaV, nspin);
        elecstate::cal_dm_psi(ParaV, wg, wfc, DM_gamma_only);
        DM_gamma_only.init_DMR(this->gd, &ucell);
        DM_gamma_only.cal_DMR();

        for (int is = 0; is < nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(charge->rho[is], charge->nrxx);
        }

        GG->transfer_DM2DtoGrid(DM_gamma_only.get_DMR_vector());
        Gint_inout inout(charge->rho, Gint_Tools::job_type::rho, nspin);
        GG->cal_gint(&inout);

        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            // for (int is = 0; is < nspin; is++)
            // {
            //     ModuleBase::GlobalFunc::ZEROS(charge->kin_r[is], charge->nrxx);
            // }
            // Gint_inout inout1(charge->kin_r, Gint_Tools::job_type::tau);
            // GG->cal_gint(&inout1);
            elecstate::lcao_cal_tau_gamma(GG, charge);
        }

        charge->renormalize_rho();
    }
    else
    {
        // calculate DMK and DMR
        elecstate::DensityMatrix<TK, double> DM(ParaV, nspin, kv->kvec_d, nk_total);
        elecstate::cal_dm_psi(ParaV, wg, wfc, DM);
        DM.init_DMR(this->gd, &ucell);
        DM.cal_DMR();

        for (int is = 0; is < nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(charge->rho[is], charge->nrxx);
        }

        GK->transfer_DM2DtoGrid(DM.get_DMR_vector());
        Gint_inout inout(charge->rho, Gint_Tools::job_type::rho, nspin);
        GK->cal_gint(&inout);

        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            // for (int is = 0; is < nspin; is++)
            // {
            //     ModuleBase::GlobalFunc::ZEROS(charge->kin_r[is], charge->nrxx);
            // }
            // Gint_inout inout1(charge->kin_r, Gint_Tools::job_type::tau);
            // GK->cal_gint(&inout1);
            elecstate::lcao_cal_tau_k(GK, charge);
        }

        charge->renormalize_rho();
    }

    // charge density symmetrization
    // this->pelec->calculate_weights();
    // this->pelec->calEBand();
    Symmetry_rho srho;
    for (int is = 0; is < nspin; is++)
    {
        srho.begin(is, *(this->charge), rho_basis, ucell.symm);
    }

}


template <typename TK, typename TR>
void RDMFT<TK, TR>::update_occNumber(const ModuleBase::matrix& occ_number_in)
{
    occ_number = (occ_number_in);
    wg = (occ_number);
    for(int ik=0; ik < wg.nr; ++ik)
    {
        for(int inb=0; inb < wg.nc; ++inb)
        {
            wg(ik, inb) *= kv->wk[ik];
            wk_fun_occNum(ik, inb) = kv->wk[ik] * occNum_func(occ_number(ik, inb), 2, XC_func_rdmft, alpha_power);
        }
    }
}


template <typename TK, typename TR>
void RDMFT<TK, TR>::update_wg(const ModuleBase::matrix& wg_in)
{
    wg = (wg_in);
    occ_number = (wg);
    for(int ik=0; ik < wg.nr; ++ik)
    {
        for(int inb=0; inb < wg.nc; ++inb)
        {
            occ_number(ik, inb) /= kv->wk[ik];
            wk_fun_occNum(ik, inb) = kv->wk[ik] * occNum_func(occ_number(ik, inb), 2, XC_func_rdmft, alpha_power);
        }
    }
}


template class RDMFT<double, double>;
template class RDMFT<std::complex<double>, double>;
template class RDMFT<std::complex<double>, std::complex<double>>;

}



