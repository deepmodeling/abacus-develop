//==========================================================
// Author: Jingang Han
// DATE : 2024-03-11
//==========================================================
#include "rdmft.h"
#include "module_rdmft/rdmft_tools.h"

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


template <typename TK, typename TR>
RDMFT<TK, TR>::RDMFT()
{

}

template <typename TK, typename TR>
RDMFT<TK, TR>::~RDMFT()
{
    delete HR_TV;
    delete HR_hartree;
    delete HR_XC;

    delete charge;
    delete Vxc_fromRI_d;
    delete Vxc_fromRI_c;

    delete V_ekinetic_potential;
    delete V_nonlocal;
    delete V_local;
    delete V_hartree;
    delete V_XC;

}

template <typename TK, typename TR>
void RDMFT<TK, TR>::init(Gint_Gamma* GG_in, Gint_k* GK_in, Parallel_Orbitals* ParaV_in, UnitCell* ucell_in, K_Vectors* kv_in, std::string XC_func_rdmft_in)
{
    GG = GG_in;
    GK = GK_in;
    ParaV = ParaV_in;
    ucell = ucell_in;
    kv = kv_in;
    nk_total = kv->nkstot_full;
    XC_func_rdmft = XC_func_rdmft_in;
    
    std::cout << "\n\n******\n" << "test class RDMFT and do rdmft_esolver.init()" << "\n******\n\n" << std::endl;

    // create desc[] and something about MPI to Eij(nbands*nbands)
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    para_Eij.set_block_size(GlobalV::NB2D);
    para_Eij.set_proc_dim(GlobalV::DSIZE);
    para_Eij.comm_2D = ParaV->comm_2D;
    para_Eij.blacs_ctxt = ParaV->blacs_ctxt;
    para_Eij.set_local2global( GlobalV::NBANDS, GlobalV::NBANDS, ofs_running, ofs_warning );
    para_Eij.set_desc( GlobalV::NBANDS, GlobalV::NBANDS, para_Eij.get_row_size(), false );

    // 
    occ_number.create(nk_total, GlobalV::NBANDS);
    wg.create(nk_total, GlobalV::NBANDS);
    wk_fun_occNum.create(nk_total, GlobalV::NBANDS);
    occNum_wfcHamiltWfc.create(nk_total, GlobalV::NBANDS);
    Etotal_n_k.create(nk_total, GlobalV::NBANDS);
    wfcHwfc_TV.create(nk_total, GlobalV::NBANDS);
    wfcHwfc_hartree.create(nk_total, GlobalV::NBANDS);
    wfcHwfc_XC.create(nk_total, GlobalV::NBANDS);

    // 
    wfc.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);   // test ParaV->nrow
    occNum_HamiltWfc.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);
    H_wfc_TV.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);
    H_wfc_hartree.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);
    H_wfc_XC.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);

    // 
    HK_TV.resize( ParaV->get_row_size()*ParaV->get_col_size() );
    HK_hartree.resize( ParaV->get_row_size()*ParaV->get_col_size() );
    HK_XC.resize( ParaV->get_row_size()*ParaV->get_col_size() );
    Eij_TV.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );
    Eij_hartree.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );
    Eij_XC.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );

    std::cout << "\n\n******\n" << "malloc for many xxx" << "\n******\n\n" << std::endl;

    // 
    HR_TV = new hamilt::HContainer<TR>(*ucell, ParaV);
    HR_hartree = new hamilt::HContainer<TR>(*ucell, ParaV);
    HR_XC = new hamilt::HContainer<TR>(*ucell, ParaV);

    // 
    // Vxc_fromRI_d = new Exx_LRI<double>(GlobalC::exx_info.info_ri);
    // Vxc_fromRI_c = new Exx_LRI<std::complex<double>>(GlobalC::exx_info.info_ri);

    std::cout << "\n\n******\n" << "malloc for HR" << "\n******\n\n" << std::endl;

}



template <typename TK, typename TR>
    template <typename T_Gint>
void RDMFT<TK, TR>::update_charge(const K_Vectors& kv, T_Gint& G_in, const Parallel_Orbitals* ParaV, const ModuleBase::matrix& wg,
                    const psi::Psi<TK>& wfc,Local_Orbital_Charge& loc, Charge& charge)  // loc can be deleted in the future
{
    if( GlobalV::GAMMA_ONLY_LOCAL )
    {
        // calculate DMK and DMR
        elecstate::DensityMatrix<TK, double> DM_gamma_only(ParaV, GlobalV::NSPIN);
        elecstate::cal_dm_psi(ParaV, wg, wfc, DM_gamma_only);
        DM_gamma_only.cal_DMR();
        DM_gamma_only.init_DMR(&GlobalC::GridD, &GlobalC::ucell);

        //this code is copying from function ElecStateLCAO<TK>::psiToRho(), in elecstate_lcao.cpp
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(charge.rho[is], charge.nrxx);
        }

        G_in.transfer_DM2DtoGrid(DM_gamma_only.get_DMR_vector());
        //double** invaild_ptr = nullptr;
        Gint_inout inout(loc.DM_R, charge.rho, Gint_Tools::job_type::rho);
        G_in.cal_gint(&inout);

        charge.renormalize_rho();
    }
    else
    {
        // calculate DMK and DMR
        elecstate::DensityMatrix<TK, double> DM(&kv, ParaV, GlobalV::NSPIN);
        elecstate::cal_dm_psi(ParaV, wg, wfc, DM);
        DM.init_DMR(&GlobalC::GridD, &GlobalC::ucell);
        DM.cal_DMR();

        // this code is copying from function ElecStateLCAO<TK>::psiToRho(), in elecstate_lcao.cpp
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(charge.rho[is], charge.nrxx);
        }

        G_in.transfer_DM2DtoGrid(DM.get_DMR_vector());
        //double** invaild_ptr = nullptr;   // use invaild_ptr replace loc.DM_R in the future
        Gint_inout inout(loc.DM_R, charge.rho, Gint_Tools::job_type::rho);  // what is Local_Orbital_Charge& loc_in? ///////////////
        G_in.cal_gint(&inout);

        charge.renormalize_rho();
    }
}





template class RDMFT<double, double>;
template class RDMFT<std::complex<double>, double>;
template class RDMFT<std::complex<double>, std::complex<double>>;



}


