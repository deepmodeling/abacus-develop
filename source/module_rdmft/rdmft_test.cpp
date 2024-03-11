//==========================================================
// Author: Jingang Han
// DATE : 2024-03-07
//==========================================================
#include "rdmft_test.h"

#include "module_base/blas_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"
#include "module_psi/psi.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"

#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>

// used by class Veff_rdmft
//#include "veff_lcao.h"
//#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/pot_local.h"
#include "module_elecstate/potentials/pot_xc.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

// for test use dgemm_
#include "module_base/matrix.h"
#include "module_base/blas_connector.h"

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"  //test

namespace rdmft
{


// // for test use dgemm_
// void printResult_dgemm()
// {
//     std::vector<double> MA = {1.1, 1.2, 3.0, 2.0, 2.1, 4.3, 3.5, 3.6, 6.0};
//     std::vector<double> MB = {0.7, 0.6, 1.0, 1.1, 2.0, 3.1};
//     std::vector<double> MC(3*2);

//     const int one_int = 1;
//     const double one_double = 1.0;
//     const double zero_double = 0.0;
//     const char N_char = 'N';
//     const char T_char = 'T';

//     int M_three = 3;
//     int N_two = 2;
//     int K_three = 3;

//     MC[0] = 1.1;

//     dgemm_(&N_char, &N_char, &N_two, &M_three, &K_three, &one_double, &MB[0], &N_two, &MA[0], &M_three, &zero_double, &MC[0], &N_two);
//     std::cout << "\n\n******\n";
//     printMatrix_vector(M_three, M_three, MA, "matrixA");
//     printMatrix_vector(M_three, N_two, MB, "matrixB");
//     printMatrix_vector(M_three, N_two, MC, "matrixC");
//     std::cout << "\n******\n\n";
// }



template class Veff_rdmft<double, double>;

template class Veff_rdmft<std::complex<double>, double>;

template class Veff_rdmft<std::complex<double>, std::complex<double>>;

// this part of the code is copying from class Veff
// initialize_HR()
template <typename TK, typename TR>
void Veff_rdmft<TK, TR>::initialize_HR(const UnitCell* ucell_in,
                                        Grid_Driver* GridD,
                                        const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("Veff", "initialize_HR");
    ModuleBase::timer::tick("Veff", "initialize_HR");

    for (int iat1 = 0; iat1 < ucell_in->nat; iat1++)
    {
        auto tau1 = ucell_in->get_tau(iat1);
        int T1, I1;
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
            const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (ucell_in->cal_dtau(iat1, iat2, R_index2).norm() * ucell_in->lat0
                < orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut())
            {
                hamilt::AtomPair<TR> tmp(iat1, iat2, R_index2.x, R_index2.y, R_index2.z, paraV);
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

    this->GK->reset_spin(GlobalV::CURRENT_SPIN);

    double* vr_eff_rdmft = nullptr;

    // calculate v_hartree(r) or v_local(r) or v_xc(r)
    if( potential_ == "hartree" )
    {   
        ModuleBase::matrix v_matrix_hartree(GlobalV::NSPIN, charge_->nrxx);
        elecstate::PotHartree potH(&rho_basis_);
        potH.cal_v_eff(charge_, ucell_, v_matrix_hartree);

        for(int is=0; is<GlobalV::NSPIN; ++is)
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
        ModuleBase::matrix v_matrix_local(1, charge_->nrxx);
        elecstate::PotLocal potL(&vloc_, &sf_, &rho_basis_);
        potL.cal_fixed_v( &v_matrix_local(0, 0) );

        // use pointer to attach v(r)
        vr_eff_rdmft = &v_matrix_local(0, 0);

        // do grid integral calculation to get HR
        Gint_inout inout(vr_eff_rdmft, 0, Gint_Tools::job_type::vlocal);
        this->GK->cal_gint(&inout);
    }
    // else if( potential_ == "XC" )
    // {
    //     elecstate::PotXC potXC();
    //     potXC.cal_v_eff(charge_, ucell_, v_matrix);
    //     ...
    // }
    else
    {
        std::cout << "\n\n!!!!!!\n there may be something wrong when use class Veff_rdmft\n\n!!!!!!\n";
    }

    // get HR for 2D-block parallel format
    this->GK->transfer_pvpR(this->hR);

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

    // this->GK->reset_spin(GlobalV::CURRENT_SPIN);

    double* vr_eff_rdmft = nullptr;

    // calculate v_hartree(r) or V_local(r) or v_xc(r)
    if( potential_ == "hartree" )
    {   
        ModuleBase::matrix v_matrix_hartree(GlobalV::NSPIN, charge_->nrxx);
        elecstate::PotHartree potH(&rho_basis_);
        potH.cal_v_eff(charge_, ucell_, v_matrix_hartree);

        for(int is=0; is<GlobalV::NSPIN; ++is)
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
        ModuleBase::matrix v_matrix_local(1, charge_->nrxx);
        elecstate::PotLocal potL(&vloc_, &sf_, &rho_basis_);
        potL.cal_fixed_v( &v_matrix_local(0, 0) );

        // use pointer to attach v(r)
        vr_eff_rdmft = &v_matrix_local(0, 0);

        // do grid integral calculation to get HR
        Gint_inout inout(vr_eff_rdmft, 0, Gint_Tools::job_type::vlocal);

        // because in gamma_only, cal_gint would not set hRGint zero first
        // so must use cal_vlocal(), and in rdmft_test.h, calculate V_hartree->contributeHR() first
        this->GG->cal_vlocal(&inout, this->LM, false);
    }
    // else if( potential_ == "XC" )
    // {
    //     elecstate::PotXC potXC();
    //     potXC.cal_v_eff(charge_, ucell_, v_matrix);
    //     ...
    // }
    else
    {
        std::cout << "\n\n!!!!!!\n there may be something wrong when use class Veff_rdmft\n\n!!!!!!\n";
    }

    // get HR for 2D-block parallel format
    this->GG->transfer_pvpR(this->hR);

    this->new_e_iteration = false;

    return;
}




template <>
void conj_psi<double>(psi::Psi<double>& wfc) {}


template <>
void HkPsi<double>(const Parallel_Orbitals* ParaV, const double& HK, const double& wfc, double& H_wfc)
{
    const int one_int = 1;
    const double one_double = 1.0;
    const double zero_double = 0.0;
    const char N_char = 'N';
    const char T_char = 'T';

    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    //because wfc(bands, basis'), H(basis, basis'), we do wfc*H^T(in the perspective of cpp, not in fortran). And get H_wfc(bands, basis) is correct.
    pdgemm_( &T_char, &N_char, &nbasis, &nbands, &nbasis, &one_double, &HK, &one_int, &one_int, ParaV->desc,
        &wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_double, &H_wfc, &one_int, &one_int, ParaV->desc_wfc );

}


template <>
void psiDotPsi<double>(const Parallel_Orbitals* ParaV, const Parallel_2D& para_Eij_in,
                        const double& wfc, const double& H_wfc, std::vector<double>& Dmn, double* wfcHwfc)
{
    const int one_int = 1;
    const double one_double = 1.0;
    const double zero_double = 0.0;
    const char N_char = 'N';
    const char T_char = 'T';

    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    const int nrow_bands = para_Eij_in.get_row_size();
    const int ncol_bands = para_Eij_in.get_col_size();

    pdgemm_( &T_char, &N_char, &nbands, &nbands, &nbasis, &one_double, &wfc, &one_int, &one_int, ParaV->desc_wfc,
            &H_wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_double, &Dmn[0], &one_int, &one_int, para_Eij_in.desc );

    for(int i=0; i<nrow_bands; ++i)
    {
        int i_global = para_Eij_in.local2global_row(i);
        for(int j=0; j<ncol_bands; ++j)
        {
            int j_global = para_Eij_in.local2global_col(j);
            if(i_global==j_global)
            {
                wfcHwfc[j_global] = std::real( Dmn[i*ncol_bands+j] );
            }
        }
    }
}


// occNum_wfcHwfc = occNum*wfcHwfc + occNum_wfcHwfc
// When symbol = 0, 1, 2, 3, 4, occNum = occNum, 0.5*occNum, g(occNum), 0.5*g(occNum), d_g(occNum)/d_occNum respectively. Default symbol=0.
void occNum_Mul_wfcHwfc(const ModuleBase::matrix& occ_number, const ModuleBase::matrix& wfcHwfc, ModuleBase::matrix& occNum_wfcHwfc,
                        int symbol, const std::string XC_func_rdmft, const double alpha)
{
    for(int ir=0; ir<occ_number.nr; ++ ir)
    {
        for(int ic=0; ic<occ_number.nc; ++ic) occNum_wfcHwfc(ir, ic) += occNum_func(occ_number(ir, ic), symbol, XC_func_rdmft, alpha) * wfcHwfc(ir, ic);
    } 
}


// Default symbol = 0 for the gradient of Etotal with respect to occupation numbers
// symbol = 1 for the relevant calculation of Etotal
void add_occNum(const ModuleBase::matrix& occ_number, const ModuleBase::matrix& wfcHwfc_TV_in, const ModuleBase::matrix& wfcHwfc_hartree_in,
                const ModuleBase::matrix& wfcHwfc_XC_in, ModuleBase::matrix& occNum_wfcHwfc, const std::string XC_func_rdmft, const double alpha, int symbol)
{
    
    occNum_wfcHwfc.zero_out();
    
    if( symbol==0 )
    {
        occNum_Mul_wfcHwfc(occ_number, wfcHwfc_XC_in, occNum_wfcHwfc, 4, XC_func_rdmft, alpha);
        occNum_wfcHwfc+=(wfcHwfc_TV_in);
        occNum_wfcHwfc+=(wfcHwfc_hartree_in);
    }
    // else if( symbol==1 )
    // {
    //     occNum_Mul_wfcHwfc(wg, wfcHwfc_TV_in, occNum_wfcHwfc);
    //     occNum_Mul_wfcHwfc(wg, wfcHwfc_hartree_in, occNum_wfcHwfc, 1);
    //     occNum_Mul_wfcHwfc(wg, wfcHwfc_XC_in, occNum_wfcHwfc, 3, XC_func_rdmft, alpha);
    // }
    else std::cout << "\n\n\n!!!!!!\nthere are something wrong when calling rdmft_test() and calculation add_occNum()\n!!!!!!\n\n\n"; 
}


// do wk*g(occNum)*wfcHwfc and add for TV, hartree, XC. This function just use once, so it can be replace and delete
void add_wfcHwfc(const std::vector<double>& wk_in, const ModuleBase::matrix& occ_number, const ModuleBase::matrix& wfcHwfc_TV_in, const ModuleBase::matrix& wfcHwfc_hartree_in,
                const ModuleBase::matrix& wfcHwfc_XC_in, ModuleBase::matrix& occNum_wfcHwfc, const std::string XC_func_rdmft, const double alpha)
{
    /****** delete when function rdmft_cal() -> class RDMFT ******/
    ModuleBase::matrix wg(occ_number);
    for(int ik=0; ik < wg.nr; ++ik)
    {
        for(int inb=0; inb < wg.nc; ++inb) wg(ik, inb) *= wk_in[ik];
    }
    ModuleBase::matrix wk_fun_occNum(occ_number.nr, occ_number.nc, true);
    for(int ik=0; ik<wg.nr; ++ik)
    {
        for(int inb=0; inb<wg.nc; ++inb) wk_fun_occNum(ik, inb) = wk_in[ik] * occNum_func(occ_number(ik, inb), 2, XC_func_rdmft, alpha);
    }
    /****** delete when function rdmft_cal() -> class RDMFT ******/

    occNum_wfcHwfc.zero_out();
    occNum_Mul_wfcHwfc(wg, wfcHwfc_TV_in, occNum_wfcHwfc);
    occNum_Mul_wfcHwfc(wg, wfcHwfc_hartree_in, occNum_wfcHwfc, 1);
    occNum_Mul_wfcHwfc(wk_fun_occNum, wfcHwfc_XC_in, occNum_wfcHwfc, 1);
}


// give certain occNum_wfcHwfc, get the corresponding energy
double sum_getEnergy(const ModuleBase::matrix& occNum_wfcHwfc)
{
    double energy = 0.0;
    for(int ir=0; ir<occNum_wfcHwfc.nr; ++ ir)
    {
        for(int ic=0; ic<occNum_wfcHwfc.nc; ++ic) energy += occNum_wfcHwfc(ir, ic);
    }
    return energy;
}


// for HF, Muller and power functional, g(eta) = eta, eta^0.5, eta^alpha respectively.
// when symbol = 0, 1, 2, 3, 4, 5, return eta, 0.5*eta, g(eta), 0.5*g(eta), d_g(eta)/d_eta, 1.0 respectively.
// Default symbol=0, XC_func_rdmft="HF", alpha=0.656
double occNum_func(double eta, int symbol, const std::string XC_func_rdmft, double alpha)
{
    if( XC_func_rdmft == "HF" ) alpha = 1.0;
    else if( XC_func_rdmft == "Muller" ) alpha = 0.5;
    else if( XC_func_rdmft == "power" ) ;

    if( symbol==0 ) return eta;
    else if ( symbol==1 ) return 0.5*eta;
    else if ( symbol==2 ) return std::pow(eta, alpha);
    else if ( symbol==3 ) return 0.5*std::pow(eta, alpha);
    else if ( symbol==4 ) return alpha*std::pow(eta, alpha-1.0);
    else if ( symbol==5 ) return 1.0;
    else 
    {
        std::cout << "\n!!!!!!\nThere may be some errors when calling wg_fun()\n!!!!!!\n";
        return eta ;
    }
}



template <typename TK, typename TR, typename T_Gint>
RDMFT<TK, TR, T_Gint>::RDMFT()
{

}

template <typename TK, typename TR, typename T_Gint>
RDMFT<TK, TR, T_Gint>::~RDMFT()
{
    
}









}


