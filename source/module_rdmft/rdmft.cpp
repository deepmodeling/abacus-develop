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
void RDMFT<TK, TR>::init(Gint_Gamma* GG_in, Gint_k* GK_in, Parallel_Orbitals* ParaV_in, UnitCell* ucell_in,
                                    K_Vectors* kv_in, std::string XC_func_rdmft_in, double alpha_power_in)
{
    GG = GG_in;
    GK = GK_in;
    ParaV = ParaV_in;
    ucell = ucell_in;
    kv = kv_in;
    nk_total = kv->nkstot_full;
    XC_func_rdmft = XC_func_rdmft_in;
    alpha_power = alpha_power_in;
    
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

    // set zero ( std::vector, ModuleBase::matrix will automatically be set to zero )
    wfc.zero_out();
    occNum_HamiltWfc.zero_out();
    H_wfc_TV.zero_out();
    H_wfc_hartree.zero_out();
    H_wfc_XC.zero_out();
    HR_TV->set_zero();
    HR_hartree->set_zero();
    HR_XC->set_zero();

    if( GlobalV::GAMMA_ONLY_LOCAL )
    {
        HR_TV->fix_gamma();
        HR_hartree->fix_gamma();
        HR_XC->fix_gamma();
    }

    // 
    // Vxc_fromRI_d = new Exx_LRI<double>(GlobalC::exx_info.info_ri);
    // Vxc_fromRI_c = new Exx_LRI<std::complex<double>>(GlobalC::exx_info.info_ri);

    std::cout << "\n\n******\n" << "malloc for HR" << "\n******\n\n" << std::endl;

}



template <typename TK, typename TR>
void RDMFT<TK, TR>::update_charge(ModuleBase::matrix& occ_number_in, const psi::Psi<TK>& wfc_in)
{
    // update occ_number, wg, wk_fun_occNum
    occ_number = (occ_number_in);
    // occ_number.zero_out();
    // occ_number+=(occ_number_in);
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
    for(int i=0; i<wfc.size(); ++i) pwfc[i] = pwfc_in[i];

    // update charge
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
            ModuleBase::GlobalFunc::ZEROS(charge->rho[is], charge->nrxx);
        }

        GG->transfer_DM2DtoGrid(DM_gamma_only.get_DMR_vector());
        //double** invaild_ptr = nullptr;   // use invaild_ptr replace loc.DM_R in the future
        Gint_inout inout(loc->DM_R, charge->rho, Gint_Tools::job_type::rho);
        GG->cal_gint(&inout);

        charge->renormalize_rho();
    }
    else
    {
        // calculate DMK and DMR
        elecstate::DensityMatrix<TK, double> DM(kv, ParaV, GlobalV::NSPIN);
        elecstate::cal_dm_psi(ParaV, wg, wfc, DM);
        DM.init_DMR(&GlobalC::GridD, &GlobalC::ucell);
        DM.cal_DMR();

        // this code is copying from function ElecStateLCAO<TK>::psiToRho(), in elecstate_lcao.cpp
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(charge->rho[is], charge->nrxx);
        }

        GK->transfer_DM2DtoGrid(DM.get_DMR_vector());
        //double** invaild_ptr = nullptr;   // use invaild_ptr replace loc.DM_R in the future
        Gint_inout inout(loc->DM_R, charge->rho, Gint_Tools::job_type::rho);  // what is Local_Orbital_Charge& loc_in? ///////////////
        GK->cal_gint(&inout);

        charge->renormalize_rho();
    }
}


template <typename TK, typename TR>
void RDMFT<TK, TR>::get_DM_XC(std::vector< std::vector<TK> >& DM_XC)
{
    // get wk_funEta_wfc = wk*g(eta)*conj(wfc)
    psi::Psi<TK> wk_funEta_wfc(wfc);
    conj_psi(wk_funEta_wfc);
    occNum_MulPsi(ParaV, wk_fun_occNum, wk_funEta_wfc, 0);

    // get the special DM_XC used in constructing V_XC
    for(int ik=0; ik<wfc.get_nk(); ++ik)
    {
        // after this, be careful with wfc.get_pointer(), we can use &wfc(ik,inbn,inbs) instead
        wfc.fix_k(ik);
        wk_funEta_wfc.fix_k(ik);
        TK* DM_Kpointer = DM_XC[ik].data();
#ifdef __MPI
        elecstate::psiMulPsiMpi(wk_funEta_wfc, wfc, DM_Kpointer, ParaV->desc_wfc, ParaV->desc);
#else
        elecstate::psiMulPsi(wk_funEta_wfc, wfc, DM_Kpointer);
#endif            
    }
}


template <typename TK, typename TR>
void RDMFT<TK, TR>::get_V_TV(LCAO_Matrix* LM_in)
{
    LM = LM_in;
    
    V_ekinetic_potential = new hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>(
        LM,
        kv->kvec_d,
        HR_TV,
        &HK_TV,
        &GlobalC::ucell,
        &GlobalC::GridD,
        ParaV
    );

    V_nonlocal = new hamilt::NonlocalNew<hamilt::OperatorLCAO<TK, TR>>(
        LM,
        kv->kvec_d,
        HR_TV,
        &HK_TV,
        &GlobalC::ucell,
        &GlobalC::GridD,
        ParaV
    );
}


template <typename TK, typename TR>
void RDMFT<TK, TR>::get_V_hartree_local(LCAO_Matrix* LM_in, const ModulePW::PW_Basis& rho_basis_in, const ModuleBase::matrix& vloc_in, const ModuleBase::ComplexMatrix& sf_in)
{
    LM = LM_in;
    if( GlobalV::GAMMA_ONLY_LOCAL )
    {
        V_local = new rdmft::Veff_rdmft<TK,TR>(
            GG,
            loc,
            LM,
            kv->kvec_d,
            charge,
            HR_TV,
            &HK_TV,
            &GlobalC::ucell,
            &GlobalC::GridD,
            ParaV,
            &rho_basis_in,
            &vloc_in,
            &sf_in,
            "local"
        );

        V_hartree = new rdmft::Veff_rdmft<TK,TR>(
            GG,
            loc,
            LM,
            kv->kvec_d,
            charge,
            HR_hartree,
            &HK_hartree,
            &GlobalC::ucell,
            &GlobalC::GridD,
            ParaV,
            &rho_basis_in,
            &vloc_in,
            &sf_in,
            "hartree"
        );
    }
    else
    {
        V_local = new rdmft::Veff_rdmft<TK,TR>(
            GK,
            loc,
            LM,
            kv->kvec_d,
            charge,
            HR_TV,
            &HK_TV,
            &GlobalC::ucell,
            &GlobalC::GridD,
            ParaV,
            &rho_basis_in,
            &vloc_in,
            &sf_in,
            "local"
        );

        V_hartree = new rdmft::Veff_rdmft<TK,TR>(
            GK,
            loc,
            LM,
            kv->kvec_d,
            charge,
            HR_hartree,
            &HK_hartree,
            &GlobalC::ucell,
            &GlobalC::GridD,
            ParaV,
            &rho_basis_in,
            &vloc_in,
            &sf_in,
            "hartree"
        );
    }
}


template <typename TK, typename TR>
void RDMFT<TK, TR>::get_V_XC()
{
    std::vector< std::vector<TK> > DM_XC(nk_total, std::vector<TK>(ParaV->nloc));
    std::vector< const std::vector<TK>* > DM_XC_pointer(nk_total);
    for(int ik=0; ik<nk_total; ++ik) DM_XC_pointer[ik] = &DM_XC[ik];
    
    get_DM_XC(DM_XC);

    if (GlobalC::exx_info.info_ri.real_number)
    {
        // transfer the DM_XC to appropriate format
        std::vector<std::map<int,std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>>> Ds_XC_d = 
            RI_2D_Comm::split_m2D_ktoR<double>(*kv, DM_XC_pointer, *ParaV);

        // provide the Ds_XC to Vxc_fromRI(V_XC)
        Vxc_fromRI_d = new Exx_LRI<double>(GlobalC::exx_info.info_ri);
        Vxc_fromRI_d->init(MPI_COMM_WORLD, *kv);
        Vxc_fromRI_d->cal_exx_ions();
        Vxc_fromRI_d->cal_exx_elec(Ds_XC_d, *ParaV);

        // when we doing V_XC.contributeHk(ik), we get HK_XC constructed by the special DM_XC
        V_XC = new hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>>(
            LM,
            HR_XC,
            &HK_XC,
            *kv,
            &Vxc_fromRI_d->Hexxs
        );
    }
    else
    {
        // transfer the DM_XC to appropriate format
        std::vector<std::map<int,std::map<std::pair<int,std::array<int,3>>,RI::Tensor<std::complex<double>>>>> Ds_XC_c = 
            RI_2D_Comm::split_m2D_ktoR<std::complex<double>>(*kv, DM_XC_pointer, *ParaV);

        // provide the Ds_XC to Vxc_fromRI(V_XC)
        Vxc_fromRI_c = new Exx_LRI<std::complex<double>>(GlobalC::exx_info.info_ri);
        Vxc_fromRI_c->init(MPI_COMM_WORLD, *kv);
        Vxc_fromRI_c->cal_exx_ions();
        Vxc_fromRI_c->cal_exx_elec(Ds_XC_c, *ParaV);

        // when we doing V_XC.contributeHk(ik), we get HK_XC constructed by the special DM_XC
        V_XC = new hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>>(
            LM,
            HR_XC,
            &HK_XC,
            *kv,
            nullptr,
            &Vxc_fromRI_c->Hexxs
        );
    }
}









template class RDMFT<double, double>;
template class RDMFT<std::complex<double>, double>;
template class RDMFT<std::complex<double>, std::complex<double>>;



}


