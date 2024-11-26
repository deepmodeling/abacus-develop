//==========================================================
// Author: Jingang Han
// DATE : 2024-03-11
//==========================================================

#include "module_rdmft/rdmft_tools.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
// used by class Veff_rdmft
#include "module_base/tool_title.h"
#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/pot_local.h"
#include "module_elecstate/potentials/pot_xc.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <cassert>


namespace rdmft
{


template <>
void conj_psi<double>(psi::Psi<double>& wfc) {}


template <>
void HkPsi<double>(const Parallel_Orbitals* ParaV, const double& HK, const double& wfc, double& H_wfc)
{
    const int one_int = 1;
    const double one_double = 1.0;
    const double zero_double = 0.0;
    const char N_char = 'N';
    const char C_char = 'C';

#ifdef __MPI
    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    //because wfc(bands, basis'), H(basis, basis'), we do wfc*H^T(in the perspective of cpp, not in fortran). And get H_wfc(bands, basis) is correct.
    pdgemm_( &C_char, &N_char, &nbasis, &nbands, &nbasis, &one_double, &HK, &one_int, &one_int, ParaV->desc,
        &wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_double, &H_wfc, &one_int, &one_int, ParaV->desc_wfc );
#endif

}


template <>
void cal_bra_op_ket<double>(const Parallel_Orbitals* ParaV, const Parallel_2D& para_Eij_in,
                                const double& wfc, const double& H_wfc, std::vector<double>& Dmn)
{
    const int one_int = 1;
    const double one_double = 1.0;
    const double zero_double = 0.0;
    const char N_char = 'N';
    const char T_char = 'T';

    const int nrow_bands = para_Eij_in.get_row_size();
    const int ncol_bands = para_Eij_in.get_col_size();

#ifdef __MPI
    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    pdgemm_( &T_char, &N_char, &nbands, &nbands, &nbasis, &one_double, &wfc, &one_int, &one_int, ParaV->desc_wfc,
            &H_wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_double, &Dmn[0], &one_int, &one_int, para_Eij_in.desc );
#endif
}


// occNum_wfcHwfc = occNum*wfcHwfc + occNum_wfcHwfc
// When symbol = 0, 1, 2, 3, 4, occNum = occNum, 0.5*occNum, g(occNum), 0.5*g(occNum), d_g(occNum)/d_occNum respectively. Default symbol=0.
void occNum_Mul_wfcHwfc(const ModuleBase::matrix& occ_number, 
                            const ModuleBase::matrix& wfcHwfc, 
                            ModuleBase::matrix& occNum_wfcHwfc,
                            int symbol, 
                            const std::string XC_func_rdmft, 
                            const double alpha)
{
    for(int ir=0; ir<occ_number.nr; ++ ir)
    {
        for(int ic=0; ic<occ_number.nc; ++ic) { occNum_wfcHwfc(ir, ic) += occNum_func(occ_number(ir, ic), symbol, XC_func_rdmft, alpha) * wfcHwfc(ir, ic);
}
    } 
}


// for the gradient of Etotal with respect to occupation numbers
void add_occNum(const K_Vectors& kv, 
                    const ModuleBase::matrix& occ_number, 
                    const ModuleBase::matrix& wfcHwfc_TV_in, 
                    const ModuleBase::matrix& wfcHwfc_hartree_in,
                    const ModuleBase::matrix& wfcHwfc_dft_XC_in, 
                    const ModuleBase::matrix& wfcHwfc_exx_XC_in, 
                    ModuleBase::matrix& occNum_wfcHwfc, 
                    const std::string XC_func_rdmft, 
                    const double alpha)
{ 
    occNum_wfcHwfc.zero_out();
    occNum_Mul_wfcHwfc(occ_number, wfcHwfc_exx_XC_in, occNum_wfcHwfc, 4, XC_func_rdmft, alpha);
    occNum_wfcHwfc+=(wfcHwfc_TV_in);
    occNum_wfcHwfc+=(wfcHwfc_hartree_in);
    occNum_wfcHwfc+=(wfcHwfc_dft_XC_in);

    // consider W_k for dE/d_occNum
    for(int ik=0; ik<occ_number.nr; ++ik)
    {
        for(int inb=0; inb<occ_number.nc; ++inb) { occNum_wfcHwfc(ik, inb) *= kv.wk[ik];
}
    } 
}


// do wk*g(occNum)*wfcHwfc and add for TV, hartree, XC. This function just use once, so it can be replace and delete
void add_wfcHwfc(const ModuleBase::matrix& wg, 
                    const ModuleBase::matrix& wk_fun_occNum, 
                    const ModuleBase::matrix& wfcHwfc_TV_in, 
                    const ModuleBase::matrix& wfcHwfc_hartree_in,
                    const ModuleBase::matrix& wfcHwfc_XC_in, 
                    ModuleBase::matrix& occNum_wfcHwfc, 
                    const std::string XC_func_rdmft, 
                    const double alpha)
{
    occNum_wfcHwfc.zero_out();
    occNum_Mul_wfcHwfc(wg, wfcHwfc_TV_in, occNum_wfcHwfc);
    occNum_Mul_wfcHwfc(wg, wfcHwfc_hartree_in, occNum_wfcHwfc, 1);
    occNum_Mul_wfcHwfc(wk_fun_occNum, wfcHwfc_XC_in, occNum_wfcHwfc, 1);
}


// give certain occNum_wfcHwfc, get the corresponding energy
double getEnergy(const ModuleBase::matrix& occNum_wfcHwfc)
{
    double energy = 0.0;
    for(int ir=0; ir<occNum_wfcHwfc.nr; ++ ir)
    {
        for(int ic=0; ic<occNum_wfcHwfc.nc; ++ic) { energy += occNum_wfcHwfc(ir, ic);
}
    }
    return energy;
}


// for HF, Muller and power functional, g(eta) = eta, eta^0.5, eta^alpha respectively.
// when symbol = 0, 1, 2, 3, 4, 5, return eta, 0.5*eta, g(eta), 0.5*g(eta), d_g(eta)/d_eta, 1.0 respectively.
// Default symbol=0, XC_func_rdmft="HF", alpha=0.656
double occNum_func(const double eta, const int symbol, const std::string XC_func_rdmft, double alpha)
{
    // if( XC_func_rdmft == "hf" || XC_func_rdmft == "default" || XC_func_rdmft == "pbe0" ) alpha = 1.0;
    // else if( XC_func_rdmft == "muller" ) alpha = 0.5;
    // else if( XC_func_rdmft == "power" || XC_func_rdmft == "wp22" || XC_func_rdmft == "cwp22" ) ;
    // else alpha = 1.0;
    if( XC_func_rdmft == "power" || XC_func_rdmft == "wp22" || XC_func_rdmft == "cwp22" ) { ; }
    else if( XC_func_rdmft == "muller" ) { alpha = 0.5; }
    else { alpha = 1.0; }

    assert(symbol <= 5);

    if( symbol==0 ) { return eta;
    } else if ( symbol==1 ) { return 0.5*eta;
    } else if ( symbol==2 ) { return std::pow(eta, alpha);
    } else if ( symbol==3 ) { return 0.5*std::pow(eta, alpha);
    } else if ( symbol==4 ) { return alpha*std::pow(eta, alpha-1.0);
    } else if ( symbol==5 ) { return 1.0;
    }
    
}




template class Veff_rdmft<double, double>;

template class Veff_rdmft<std::complex<double>, double>;

template class Veff_rdmft<std::complex<double>, std::complex<double>>;

// this part of the code is copying from class Veff
// initialize_HR()
template <typename TK, typename TR>
void Veff_rdmft<TK, TR>::initialize_HR(const UnitCell* ucell_in,
                                        Grid_Driver* GridD)
{
    ModuleBase::TITLE("Veff", "initialize_HR");
    ModuleBase::timer::tick("Veff", "initialize_HR");

    this->nspin = PARAM.inp.nspin;
    auto* paraV = this->hR->get_paraV();// get parallel orbitals from HR
    // TODO: if paraV is nullptr, AtomPair can not use paraV for constructor, I will repair it in the future.

    for (int iat1 = 0; iat1 < ucell_in->nat; iat1++)
    {
        auto tau1 = ucell_in->get_tau(iat1);
        int T1 = 0;
        int I1 = 0;
        ucell_in->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell_in, tau1, T1, I1, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T2 = adjs.ntype[ad1];
            const int I2 = adjs.natom[ad1];
            const int iat2 = ucell_in->itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad1];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (ucell_in->cal_dtau(iat1, iat2, R_index2).norm() * ucell_in->lat0
                < orb_cutoff_[T1] + orb_cutoff_[T2])
            {
                hamilt::AtomPair<TR> tmp(iat1, iat2, R_index2, paraV);
                this->hR->insert_pair(tmp);
            }
        }
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    this->hR->allocate(nullptr, true);

    ModuleBase::timer::tick("Veff", "initialize_HR");
}


// this part of the code is copying from class Veff and do some modifications.
template<typename TK, typename TR>
void Veff_rdmft<TK, TR>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::tick("Veff", "contributeHR");

    this->GK->reset_spin(this->current_spin);

    double* vr_eff_rdmft = nullptr;

    // calculate v_hartree(r) or v_local(r) or v_xc(r)
    if( potential_ == "hartree" )
    {   
        ModuleBase::matrix v_matrix_hartree(this->nspin, charge_->nrxx);
        elecstate::PotHartree potH(rho_basis_);
        potH.cal_v_eff(charge_, ucell, v_matrix_hartree);

        for(int is=0; is<this->nspin; ++is)
        {
            // use pointer to attach v(r) for current spin
            vr_eff_rdmft = &v_matrix_hartree(is, 0);

            // do grid integral calculation to get HR
            Gint_inout inout(vr_eff_rdmft, is, Gint_Tools::job_type::vlocal);
            this->GK->cal_gint(&inout);
        }
    }
    else if( potential_ == "local" )
    {   
        double vlocal_of_0 = 0.0;
        ModuleBase::matrix v_matrix_local(1, charge_->nrxx);
        elecstate::PotLocal potL(vloc_, sf_, rho_basis_, vlocal_of_0);
        potL.cal_fixed_v( &v_matrix_local(0, 0) );

        // use pointer to attach v(r)
        vr_eff_rdmft = &v_matrix_local(0, 0);

        // do grid integral calculation to get HR
        Gint_inout inout(vr_eff_rdmft, 0, Gint_Tools::job_type::vlocal);
        this->GK->cal_gint(&inout);
    }
    else if( potential_ == "xc" )
    {
        // meta-gga type has not been considered yet !!!

        ModuleBase::matrix vofk = *vloc_;
        vofk.zero_out();
        ModuleBase::matrix v_matrix_XC(this->nspin, charge_->nrxx);
        elecstate::PotXC potXC(rho_basis_, etxc, vtxc, &vofk);
        potXC.cal_v_eff(charge_, ucell, v_matrix_XC);

        // if need meta-GGA, go to study veff_lcao.cpp and modify the code
        for(int is=0; is<this->nspin; ++is)
        {
            // use pointer to attach v(r) for current spin
            vr_eff_rdmft = &v_matrix_XC(is, 0);

            // do grid integral calculation to get HR
            Gint_inout inout(vr_eff_rdmft, is, Gint_Tools::job_type::vlocal);
            this->GK->cal_gint(&inout);
        }
    }
    else
    {
        std::cout << "\n\n!!!!!!\n there may be something wrong when use class Veff_rdmft\n\n!!!!!!\n";
    }

    // get HR for 2D-block parallel format
    // this->GK->transfer_pvpR(this->hR);
    this->GK->transfer_pvpR(this->hR,this->ucell,this->gd);

    if(this->nspin == 2) { this->current_spin = 1 - this->current_spin; }

    ModuleBase::timer::tick("Veff", "contributeHR");
    return;
}

// this part of the code is copying from class Veff and do some modifications.
// special case of gamma-only
template<>
void Veff_rdmft<double, double>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::tick("Veff", "contributeHR");

    // this->GK->reset_spin(this->current_spin);

    double* vr_eff_rdmft = nullptr;

    // calculate v_hartree(r) or V_local(r) or v_xc(r)
    if( potential_ == "hartree" )
    {   
        ModuleBase::matrix v_matrix_hartree(this->nspin, charge_->nrxx);
        elecstate::PotHartree potH(rho_basis_);
        potH.cal_v_eff(charge_, ucell, v_matrix_hartree);

        for(int is=0; is<this->nspin; ++is)
        {
            // use pointer to attach v(r) for current spin
            vr_eff_rdmft = &v_matrix_hartree(is, 0);

            // do grid integral calculation to get HR
            Gint_inout inout(vr_eff_rdmft, is, Gint_Tools::job_type::vlocal);
            this->GG->cal_gint(&inout);
        }
    }
    else if( potential_ == "local" )
    {   
        double vlocal_of_0 = 0.0;
        ModuleBase::matrix v_matrix_local(1, charge_->nrxx);
        elecstate::PotLocal potL(vloc_, sf_, rho_basis_, vlocal_of_0);
        potL.cal_fixed_v( &v_matrix_local(0, 0) );

        // use pointer to attach v(r)
        vr_eff_rdmft = &v_matrix_local(0, 0);

        // do grid integral calculation to get HR
        Gint_inout inout(vr_eff_rdmft, 0, Gint_Tools::job_type::vlocal);

        // because in gamma_only, cal_gint would not set hRGint zero first
        // so must use cal_vlocal(), and in rdmft_test.h, calculate V_hartree->contributeHR() first

        this->GG->cal_vlocal(&inout, false);  // cal_gint ???
    }
    else if( potential_ == "xc" )
    {
        // meta-gga type has not been considered yet !!!

        ModuleBase::matrix vofk = *vloc_;
        vofk.zero_out();
        ModuleBase::matrix v_matrix_XC(this->nspin, charge_->nrxx);
        elecstate::PotXC potXC(rho_basis_, etxc, vtxc, &vofk);
        potXC.cal_v_eff(charge_, ucell, v_matrix_XC);
        
        for(int is=0; is<this->nspin; ++is)
        {
            // use pointer to attach v(r) for current spin
            vr_eff_rdmft = &v_matrix_XC(is, 0);

            // do grid integral calculation to get HR
            Gint_inout inout(vr_eff_rdmft, is, Gint_Tools::job_type::vlocal);
            this->GG->cal_gint(&inout);
        }
    }
    else
    {
        std::cout << "\n\n!!!!!!\n there may be something wrong when use class Veff_rdmft\n\n!!!!!!\n";
    }

    // get HR for 2D-block parallel format
    this->GG->transfer_pvpR(this->hR,this->ucell);

    this->new_e_iteration = false;

    if(this->nspin == 2) this->current_spin = 1 - this->current_spin;

    return;
}



}



