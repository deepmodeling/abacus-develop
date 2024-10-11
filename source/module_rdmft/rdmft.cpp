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

#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_cell/module_symmetry/symmetry.h"
// #include "module_hamilt_general/module_xc/xc_functional.h"

#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>

// used by class Veff_rdmft
//#include "veff_lcao.h"
//#include "module_base/timer.h"
// #include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/veff_lcao.h"
// #include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_base/tool_title.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/pot_local.h"
#include "module_elecstate/potentials/pot_xc.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_elecstate/elecstate_lcao.h"

// for test use dgemm_
#include "module_base/matrix.h"
#include "module_base/blas_connector.h"

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"  //test

#include "module_elecstate/module_charge/symmetry_rho.h"




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
    // delete HR_XC;
    delete HR_exx_XC;
    // delete HR_local;
    delete HR_dft_XC;

    delete Vxc_fromRI_d;
    delete Vxc_fromRI_c;

    delete V_ekinetic_potential;
    delete V_nonlocal;
    delete V_local;
    delete V_hartree;
    // delete V_XC;
    delete V_exx_XC;
    delete V_dft_XC;
    delete V_hartree_XC;
}

template <typename TK, typename TR>
void RDMFT<TK, TR>::init(Gint_Gamma& GG_in, Gint_k& GK_in, Parallel_Orbitals& ParaV_in, UnitCell& ucell_in,
                                    K_Vectors& kv_in, elecstate::ElecState& pelec_in, LCAO_Orbitals& orb_in, TwoCenterBundle& two_center_bundle_in, std::string XC_func_rdmft_in, double alpha_power_in)
{
    GG = &GG_in;
    GK = &GK_in;
    ParaV = &ParaV_in;
    ucell = &ucell_in;
    kv = &kv_in;
    charge = pelec_in.charge;
    pelec = &pelec_in;
    orb = &orb_in;
    two_center_bundle = &two_center_bundle_in;
    XC_func_rdmft = XC_func_rdmft_in;
    alpha_power = alpha_power_in;

    // if (ModuleSymmetry::Symmetry::symm_flag == -1) nk_total = kv->nkstot_full;
    // else nk_total = kv->nks;

    nk_total = ModuleSymmetry::Symmetry::symm_flag == -1 ? kv->nkstot_full: kv->nks;
    nspin = PARAM.inp.nspin;

    // XC_func_rdmft = "power"; // just for test 
    // alpha_power = 0.525;
    std::cout << "\n\n\n******\nXC-functional in rdmft: " << XC_func_rdmft << "\n******\n\n\n" << std::endl;
    // XC_func_rdmft = "hf";
    // std::cout << "\n\n\n******\nXC-functional in rdmft: " << XC_func_rdmft << "\n******\n\n\n" << std::endl;
    // std::cout << "\n\n\n******\nXC-functional in GlobalC::atom: " << GlobalC::ucell.atoms[0].ncpp.xc_func << "\n******\n\n\n" << std::endl;
    // if( XC_func_rdmft == "default" ) XC_func_rdmft = "default";

    // create desc[] and something about MPI to Eij(nbands*nbands)
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    // para_Eij.set_block_size(GlobalV::NB2D);
    // para_Eij.set_proc_dim(GlobalV::DSIZE);
    // para_Eij.comm_2D = ParaV->comm_2D;
    // para_Eij.blacs_ctxt = ParaV->blacs_ctxt;
    // para_Eij.set_local2global( GlobalV::NBANDS, GlobalV::NBANDS, ofs_running, ofs_warning );
    // para_Eij.set_desc( GlobalV::NBANDS, GlobalV::NBANDS, para_Eij.get_row_size(), false );
    para_Eij.set(nbands_total, nbands_total, PARAM.inp.nb2d, ParaV->blacs_ctxt);

    // // learn from "module_hamilt_lcao/hamilt_lcaodft/LCAO_init_basis.cpp"
    // int try_nb = para_Eij.init(GlobalV::NBANDS, GlobalV::NBANDS, PARAM.inp.nb2d, DIAG_WORLD); // DIAG_WORLD is wrong
    // try_nb += para_Eij.set_nloc_wfc_Eij(GlobalV::NBANDS, ofs_running, ofs_warning);
    // if (try_nb != 0)
    // {
    //     para_Eij.set(GlobalV::NBANDS, GlobalV::NBANDS, 1, para_Eij.blacs_ctxt);
    //     try_nb = para_Eij.set_nloc_wfc_Eij(GlobalV::NBANDS, GlobalV::ofs_running, GlobalV::ofs_warning);
    // }
    // para_Eij.set_desc_wfc_Eij(GlobalV::NBANDS, GlobalV::NBANDS, para_Eij.nrow);


    // 
    occ_number.create(nk_total, nbands_total);
    wg.create(nk_total, nbands_total);
    wk_fun_occNum.create(nk_total, nbands_total);
    occNum_wfcHamiltWfc.create(nk_total, nbands_total);
    Etotal_n_k.create(nk_total, nbands_total);
    wfcHwfc_TV.create(nk_total, nbands_total);
    wfcHwfc_hartree.create(nk_total, nbands_total);
    wfcHwfc_XC.create(nk_total, nbands_total);
    wfcHwfc_exx_XC.create(nk_total, nbands_total);
    wfcHwfc_dft_XC.create(nk_total, nbands_total);

    // 
    wfc.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);   // test ParaV->nrow
    occNum_HamiltWfc.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);
    H_wfc_TV.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);
    H_wfc_hartree.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);
    H_wfc_XC.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);
    H_wfc_exx_XC.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);
    H_wfc_dft_XC.resize(nk_total, ParaV->ncol_bands, ParaV->nrow);

    // 
    // HK_TV.resize( ParaV->get_row_size()*ParaV->get_col_size() );
    // HK_hartree.resize( ParaV->get_row_size()*ParaV->get_col_size() );
    // HK_exx_XC.resize( ParaV->get_row_size()*ParaV->get_col_size() );
    // HK_dft_XC.resize( ParaV->get_row_size()*ParaV->get_col_size() );
    //
    hsk_TV = new hamilt::HS_Matrix_K<TK>(ParaV, true);
    hsk_hartree = new hamilt::HS_Matrix_K<TK>(ParaV, true);
    hsk_dft_XC = new hamilt::HS_Matrix_K<TK>(ParaV, true);
    hsk_exx_XC = new hamilt::HS_Matrix_K<TK>(ParaV, true);

    HK_XC.resize( ParaV->get_row_size()*ParaV->get_col_size() );
    HK_local.resize( ParaV->get_row_size()*ParaV->get_col_size() );
    // HK_RDMFT_pass.resize(nk_total, ParaV->get_row_size(), ParaV->get_col_size());
    // HK_XC_pass.resize(nk_total, ParaV->get_row_size(), ParaV->get_col_size());


    Eij_TV.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );
    Eij_hartree.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );
    Eij_XC.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );
    Eij_exx_XC.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );

    // 
    HR_TV = new hamilt::HContainer<TR>(*ucell, ParaV);
    HR_hartree = new hamilt::HContainer<TR>(*ucell, ParaV);
    // HR_XC = new hamilt::HContainer<TR>(*ucell, ParaV);
    HR_exx_XC = new hamilt::HContainer<TR>(*ucell, ParaV);
    HR_dft_XC = new hamilt::HContainer<TR>(*ucell, ParaV);
    // HR_local = new hamilt::HContainer<TR>(*ucell, ParaV);

    // set zero ( std::vector, ModuleBase::matrix will automatically be set to zero )
    wfc.zero_out();
    occNum_HamiltWfc.zero_out();
    H_wfc_TV.zero_out();
    H_wfc_hartree.zero_out();
    H_wfc_XC.zero_out();
    H_wfc_exx_XC.zero_out();
    H_wfc_dft_XC.zero_out();
    HR_TV->set_zero();         // HR->set_zero() might be delete here, test on Gamma_only in the furure 
    HR_hartree->set_zero();
    // HR_XC->set_zero();
    HR_exx_XC->set_zero();
    HR_dft_XC->set_zero();
    // HR_local->set_zero();

    if( GlobalC::exx_info.info_global.cal_exx )
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            Vxc_fromRI_d = new Exx_LRI<double>(GlobalC::exx_info.info_ri);
            Vxc_fromRI_d->init(MPI_COMM_WORLD, *kv, *orb);
        }
        else
        {
            Vxc_fromRI_c = new Exx_LRI<std::complex<double>>(GlobalC::exx_info.info_ri);
            Vxc_fromRI_c->init(MPI_COMM_WORLD, *kv, *orb);
        }
    }

    if( PARAM.inp.gamma_only )
    {
        HR_TV->fix_gamma();
        HR_hartree->fix_gamma();
        // HR_XC->fix_gamma();
        HR_exx_XC->fix_gamma();
        HR_dft_XC->fix_gamma();
    }

}


template <typename TK, typename TR>
void RDMFT<TK, TR>::update_ion(UnitCell& ucell_in, ModulePW::PW_Basis& rho_basis_in,
                                ModuleBase::matrix& vloc_in, ModuleBase::ComplexMatrix& sf_in)
{
    ucell = &ucell_in;
    // LM = &LM_in;
    rho_basis = &rho_basis_in;
    vloc = &vloc_in;
    sf = &sf_in;
    // loc = &loc_in;


    HR_TV->set_zero();
    this->cal_V_TV();

    if( GlobalC::exx_info.info_global.cal_exx )
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            Vxc_fromRI_d->cal_exx_ions();
        }
        else
        {
            Vxc_fromRI_c->cal_exx_ions();
        }
    }

    std::cout << "\n\n\n******\ndo rdmft_esolver.update_ion() successfully\n******\n\n\n" << std::endl;
}


template <typename TK, typename TR>
void RDMFT<TK, TR>::update_elec(const ModuleBase::matrix& occ_number_in, const psi::Psi<TK>& wfc_in, const Charge* charge_in)
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
    for(int i=0; i<wfc.size(); ++i) pwfc[i] = pwfc_in[i];

    // update charge
    this->update_charge();

    // "default" = "pbe"
    // if(  (XC_func_rdmft != "hf" && XC_func_rdmft != "muller" && XC_func_rdmft != "power") || this->cal_E_type != 1 )
    if( this->cal_E_type != 1 )
    {
        // the second cal_E_type need the complete pot to get effctive_V to calEband and so on.
        this->pelec->pot->update_from_charge(charge, ucell);
    }

    this->cal_V_hartree();
    this->cal_V_XC();
    this->cal_Hk_Hpsi();
}


// this code is copying from function ElecStateLCAO<TK>::psiToRho(), in elecstate_lcao.cpp
template <typename TK, typename TR>
void RDMFT<TK, TR>::update_charge()
{
    if( PARAM.inp.gamma_only )
    {
        // calculate DMK and DMR
        elecstate::DensityMatrix<TK, double> DM_gamma_only(ParaV, nspin);
        elecstate::cal_dm_psi(ParaV, wg, wfc, DM_gamma_only);
        DM_gamma_only.init_DMR(&GlobalC::GridD, &GlobalC::ucell);
        DM_gamma_only.cal_DMR();

        for (int is = 0; is < nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(charge->rho[is], charge->nrxx);
        }

        GG->transfer_DM2DtoGrid(DM_gamma_only.get_DMR_vector());
        Gint_inout inout(charge->rho, Gint_Tools::job_type::rho);
        GG->cal_gint(&inout);

        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            // for (int is = 0; is < nspin; is++)
            // {
            //     ModuleBase::GlobalFunc::ZEROS(charge->kin_r[is], charge->nrxx);
            // }
            // Gint_inout inout1(charge->kin_r, Gint_Tools::job_type::tau);
            // GG->cal_gint(&inout1);
            this->pelec->cal_tau(wfc);
        }

        charge->renormalize_rho();
    }
    else
    {
        // calculate DMK and DMR
        elecstate::DensityMatrix<TK, double> DM(kv, ParaV, nspin);
        elecstate::cal_dm_psi(ParaV, wg, wfc, DM);
        DM.init_DMR(&GlobalC::GridD, &GlobalC::ucell);
        DM.cal_DMR();

        for (int is = 0; is < nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(charge->rho[is], charge->nrxx);
        }

        GK->transfer_DM2DtoGrid(DM.get_DMR_vector());
        Gint_inout inout(charge->rho, Gint_Tools::job_type::rho);
        GK->cal_gint(&inout);

        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            // for (int is = 0; is < nspin; is++)
            // {
            //     ModuleBase::GlobalFunc::ZEROS(charge->kin_r[is], charge->nrxx);
            // }
            // Gint_inout inout1(charge->kin_r, Gint_Tools::job_type::tau);
            // GK->cal_gint(&inout1);
            this->pelec->cal_tau(wfc);
        }

        charge->renormalize_rho();
    }

    /*********** what's this? When we use PBE-functional, this can't be deleted *************/
    // this->pelec->calculate_weights();
    // this->pelec->calEBand();
    Symmetry_rho srho;
    for (int is = 0; is < nspin; is++)
    {
        srho.begin(is, *(this->charge), rho_basis, GlobalC::ucell.symm);
    }
    /*********** what's this? When we use PBE-functional, this can't be deleted *************/

    // // what this?
    // if(GlobalV::VL_IN_H)
    // {
    //     // update Gint_K
    //     if (!PARAM.inp.gamma_only)
    //     {
    //         this->UHM.GK.renew();
    //     }
    //     // update real space Hamiltonian
    //     // this->p_hamilt->refresh();
    // }

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


// template <typename TK, typename TR>
// void RDMFT<TK, TR>::update_wg(const ModuleBase::matrix& wg_in)
// {
//     wg = (wg_in);
//     occ_number = (wg);
//     for(int ik=0; ik < wg.nr; ++ik)
//     {
//         for(int inb=0; inb < wg.nc; ++inb)
//         {
//             occ_number(ik, inb) /= kv->wk[ik];
//             wk_fun_occNum(ik, inb) = kv->wk[ik] * occNum_func(occ_number(ik, inb), 2, XC_func_rdmft, alpha_power);
//         }
//     }
// }


template <typename TK, typename TR>
void RDMFT<TK, TR>::get_DM_XC(std::vector< std::vector<TK> >& DM_XC)
{
    // get wk_funEta_wfc = wk*g(eta)*conj(wfc)
    psi::Psi<TK> wk_funEta_wfc(wfc);
    conj_psi(wk_funEta_wfc);
    occNum_MulPsi(ParaV, wk_fun_occNum, wk_funEta_wfc, 0);

    // get the special DM_XC used in constructing V_exx_XC
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
void RDMFT<TK, TR>::cal_V_TV()
{
    HR_TV->set_zero();
    
    V_ekinetic_potential = new hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>(
        hsk_TV,
        kv->kvec_d,
        HR_TV,
        &GlobalC::ucell,
        orb->cutoffs(),
        &GlobalC::GridD,
        two_center_bundle->kinetic_orb.get()
    );

    V_nonlocal = new hamilt::NonlocalNew<hamilt::OperatorLCAO<TK, TR>>(
        hsk_TV,
        kv->kvec_d,
        HR_TV,
        &GlobalC::ucell,
        orb->cutoffs(),
        &GlobalC::GridD,
        two_center_bundle->overlap_orb_beta.get()
    );

    if( PARAM.inp.gamma_only )
    {
        V_local = new rdmft::Veff_rdmft<TK,TR>(
            GG,
            hsk_TV,
            kv->kvec_d,
            this->pelec->pot,
            HR_TV,
            &GlobalC::ucell,
            orb->cutoffs(),
            &GlobalC::GridD,

            charge,
            rho_basis,
            vloc,
            sf,
            "local"
        );
    }
    else
    {
        V_local = new rdmft::Veff_rdmft<TK,TR>(
            GK,
            hsk_TV,
            kv->kvec_d,
            this->pelec->pot,
            HR_TV,
            &GlobalC::ucell,
            orb->cutoffs(),
            &GlobalC::GridD,

            charge,
            rho_basis,
            vloc,
            sf,
            "local"
        );
    }

    // update HR_TV in ion-step, now HR_TV has the HR of V_ekinetic + V_nonlcao + V_local
    V_ekinetic_potential->contributeHR();
    V_nonlocal->contributeHR();
    V_local->contributeHR();

}


template <typename TK, typename TR>
void RDMFT<TK, TR>::cal_V_hartree()
{
    HR_hartree->set_zero();

    if( PARAM.inp.gamma_only )
    {
        V_hartree = new rdmft::Veff_rdmft<TK,TR>(
            GG,
            hsk_hartree,
            kv->kvec_d,
            this->pelec->pot,
            HR_hartree,
            &GlobalC::ucell,
            orb->cutoffs(),
            &GlobalC::GridD,

            charge,
            rho_basis,
            vloc,
            sf,
            "hartree"
        );
    }
    else
    {
        // this can be optimized, use potHartree.update_from_charge()
        V_hartree = new rdmft::Veff_rdmft<TK,TR>(
            GK,
            hsk_hartree,
            kv->kvec_d,
            this->pelec->pot,
            HR_hartree,
            &GlobalC::ucell,
            orb->cutoffs(),
            &GlobalC::GridD,

            charge,
            rho_basis,
            vloc,
            sf,
            "hartree"
        );
    }

    // in gamma only, must calculate HR_hartree before HR_local
    // HR_exx_XC get from another way, so don't need to do this 
    V_hartree->contributeHR();

    // // update HR_local in e-step, now HR_TV has the HR of V_ekinetic + V_nonlcao + V_local, 
    // V_local->contributeHR();
    // HR_local->add(*HR_TV);  // now HR_local has the HR of V_ekinetic + V_nonlcao + V_local

}


template <typename TK, typename TR>
void RDMFT<TK, TR>::cal_V_XC()
{
    HR_dft_XC->set_zero();
    HR_exx_XC->set_zero();

    std::vector< std::vector<TK> > DM_XC(nk_total, std::vector<TK>(ParaV->nloc));
    std::vector< const std::vector<TK>* > DM_XC_pointer(nk_total);
    for(int ik=0; ik<nk_total; ++ik) DM_XC_pointer[ik] = &DM_XC[ik];
    
    get_DM_XC(DM_XC);

    // //test

    DM_XC_pass = DM_XC;

    // elecstate::DensityMatrix<TK, double> DM_test(kv, ParaV, nspin);
    // elecstate::cal_dm_psi(ParaV, wg, wfc, DM_test);
    // DM_test.init_DMR(&GlobalC::GridD, &GlobalC::ucell);
    // DM_test.cal_DMR();

    // // compare DM_XC and DM get in update_charge(or ABACUS)
    // std::cout << "\n\ntest DM_XC - DM in ABACUS: \n" << std::endl;
    // double DM_XC_minus_DMtest = 0.0;
    // for(int ik=0; ik<nk_total; ++ik)
    // {
    //     TK* dmk_pointer = DM_test.get_DMK_pointer(ik);
    //     for(int iloc=0; iloc<ParaV->nloc; ++iloc)
    //     {
    //         double test = std::abs(DM_XC[ik][iloc] - dmk_pointer[iloc]);
    //         DM_XC_minus_DMtest += test;
    //         if( test > 1e-16 )
    //         {
    //             std::cout << "\nik, iloc, minus[ik][iloc]: " << ik << " " << iloc << " " << test << std::endl; 
    //         }
    //     }
    // }
    // std::cout << "\nsum of DM_XC - DM in ABACUS: " << DM_XC_minus_DMtest << std::endl;

    if( XC_func_rdmft != "hf" && XC_func_rdmft != "muller" && XC_func_rdmft != "power" )
    {
        if( PARAM.inp.gamma_only )
        {
            // this can be optimized, use potXC.update_from_charge()
            V_dft_XC = new rdmft::Veff_rdmft<TK,TR>(
                GG,
                hsk_dft_XC,
                kv->kvec_d,
                this->pelec->pot,
                HR_dft_XC,
                &GlobalC::ucell,
                orb->cutoffs(),
                &GlobalC::GridD,

                charge,
                rho_basis,
                vloc,
                sf,
                "xc",
                &etxc,
                &vtxc
            );
            // V_dft_XC->contributeHR();
        }
        else
        {   
            // this can be optimized, use potXC.update_from_charge()
            V_dft_XC = new rdmft::Veff_rdmft<TK,TR>(
                GK,
                hsk_dft_XC,
                kv->kvec_d,
                this->pelec->pot,
                HR_dft_XC,
                &GlobalC::ucell,
                orb->cutoffs(),
                &GlobalC::GridD,

                charge,
                rho_basis,
                vloc,
                sf,
                "xc",
                &etxc,
                &vtxc
            );
            // V_dft_XC->contributeHR();
        }
        V_dft_XC->contributeHR();
    }
    
    if(GlobalC::exx_info.info_global.cal_exx)
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            // transfer the DM_XC to appropriate format
            std::vector<std::map<int,std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>>> Ds_XC_d = 
                RI_2D_Comm::split_m2D_ktoR<double>(*kv, DM_XC_pointer, *ParaV, nspin);

            // provide the Ds_XC to Vxc_fromRI(V_exx_XC)
            // Vxc_fromRI_d = new Exx_LRI<double>(GlobalC::exx_info.info_ri);
            // Vxc_fromRI_d->init(MPI_COMM_WORLD, *kv);
            // Vxc_fromRI_d->cal_exx_ions();
            Vxc_fromRI_d->cal_exx_elec(Ds_XC_d, *ParaV);

            // when we doing V_exx_XC.contributeHk(ik), we get HK_XC constructed by the special DM_XC
            V_exx_XC = new hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>>(
                hsk_exx_XC,
                HR_exx_XC,
                *kv,
                &Vxc_fromRI_d->Hexxs
            );
        }
        else
        {
            // transfer the DM_XC to appropriate format
            std::vector<std::map<int,std::map<std::pair<int,std::array<int,3>>,RI::Tensor<std::complex<double>>>>> Ds_XC_c = 
                RI_2D_Comm::split_m2D_ktoR<std::complex<double>>(*kv, DM_XC_pointer, *ParaV, nspin);

            // provide the Ds_XC to Vxc_fromRI(V_exx_XC)
            // Vxc_fromRI_c = new Exx_LRI<std::complex<double>>(GlobalC::exx_info.info_ri);
            // Vxc_fromRI_c->init(MPI_COMM_WORLD, *kv);
            // Vxc_fromRI_c->cal_exx_ions();
            Vxc_fromRI_c->cal_exx_elec(Ds_XC_c, *ParaV);

            // when we doing V_exx_XC.contributeHk(ik), we get HK_XC constructed by the special DM_XC
            V_exx_XC = new hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>>(
                hsk_exx_XC,
                HR_exx_XC,
                *kv,
                nullptr,
                &Vxc_fromRI_c->Hexxs
            );
        }
        V_exx_XC->contributeHR();
    }
}


template <typename TK, typename TR>
void RDMFT<TK, TR>::cal_Hk_Hpsi()
{
    /****** get occNum_wfcHamiltWfc, occNum_HamiltWfc ******/
    // HK_RDMFT_pass.reset();

    // double XC_minus_XC = 0.0;
    // std::cout << "\n\ntest V_exx_XC in rdmft.cpp: " << std::endl;
    // HK_XC_pass.reset();

    //calculate Hwfc, wfcHwfc for each potential
    for(int ik=0; ik<nk_total; ++ik)
    {
        hsk_TV->set_zero_hk();
        hsk_hartree->set_zero_hk();
        set_zero_vector(HK_XC);

        // get the HK with ik-th k vector, the result is stored in HK_TV, HK_hartree and HK_XC respectively
        V_local->contributeHk(ik);
        V_hartree->contributeHk(ik);

        // get H(k) * wfc
        HkPsi( ParaV, hsk_TV->get_hk()[0], wfc(ik, 0, 0), H_wfc_TV(ik, 0, 0));
        HkPsi( ParaV, hsk_hartree->get_hk()[0], wfc(ik, 0, 0), H_wfc_hartree(ik, 0, 0));

        // get wfc * H(k)_wfc
        psiDotPsi( ParaV, para_Eij, wfc(ik, 0, 0), H_wfc_TV(ik, 0, 0), Eij_TV, &(wfcHwfc_TV(ik, 0)) );
        psiDotPsi( ParaV, para_Eij, wfc(ik, 0, 0), H_wfc_hartree(ik, 0, 0), Eij_hartree, &(wfcHwfc_hartree(ik, 0)) );

        if(GlobalC::exx_info.info_global.cal_exx)
        {
            // set_zero_vector(HK_exx_XC);
            hsk_exx_XC->set_zero_hk();

            V_exx_XC->contributeHk(ik);
            HkPsi( ParaV, hsk_exx_XC->get_hk()[0], wfc(ik, 0, 0), H_wfc_exx_XC(ik, 0, 0));
            psiDotPsi( ParaV, para_Eij, wfc(ik, 0, 0), H_wfc_exx_XC(ik, 0, 0), Eij_exx_XC, &(wfcHwfc_exx_XC(ik, 0)) );
            
            for(int iloc=0; iloc<HK_XC.size(); ++iloc) HK_XC[iloc] += hsk_exx_XC->get_hk()[iloc];
        }

        if( XC_func_rdmft != "hf" && XC_func_rdmft != "muller" && XC_func_rdmft != "power" )
        {
            // set_zero_vector(HK_dft_XC);
            hsk_dft_XC->set_zero_hk();

            V_dft_XC->contributeHk(ik);
            HkPsi( ParaV, hsk_dft_XC->get_hk()[0], wfc(ik, 0, 0), H_wfc_dft_XC(ik, 0, 0));
            psiDotPsi( ParaV, para_Eij, wfc(ik, 0, 0), H_wfc_dft_XC(ik, 0, 0), Eij_exx_XC, &(wfcHwfc_dft_XC(ik, 0)) );
            
            for(int iloc=0; iloc<HK_XC.size(); ++iloc) HK_XC[iloc] += hsk_dft_XC->get_hk()[iloc];
        }
        // elseif()

        // // store HK_RDMFT
        // for(int ir=0; ir<HK_RDMFT_pass.nr; ++ir)
        // {
        //     for(int ic=0; ic<HK_RDMFT_pass.nc; ++ic)
        //     {
        //         HK_RDMFT_pass[ik](ir, ic) = HK_TV[ic * ParaV->get_col_size() + ir]
        //                                 + HK_hartree[ic * ParaV->get_col_size() + ir]
        //                                 + HK_XC[ic * ParaV->get_col_size() + ir];
        //         // HK_XC_pass[ik](ir, ic) = HK_XC[ic * ParaV->get_col_size() + ir];
        //     }
        // }


        // using them to the gradient of Etotal is not correct when do hybrid calculation, it's correct just for exx-type functional
        // HkPsi( ParaV, HK_XC[0], wfc(ik, 0, 0), H_wfc_XC(ik, 0, 0));
        // psiDotPsi( ParaV, para_Eij, wfc(ik, 0, 0), H_wfc_XC(ik, 0, 0), Eij_XC, &(wfcHwfc_XC(ik, 0)) );
        
        // let H(k)=0 to storing next one, H(k+1)
        // set_zero_vector(HK_TV);
        // set_zero_vector(HK_hartree);

    }

    // std::cout << "\n\nsum of XC_minus_XC: " << XC_minus_XC << "\n\n" << std::endl;

}


template <typename TK, typename TR>
double RDMFT<TK, TR>::cal_E_gradient()
{
    /****** get occNum_wfcHamiltWfc, occNum_HamiltWfc and Etotal ******/

    // !this would transfer the value of H_wfc_TV, H_wfc_hartree, H_wfc_XC --> occNum_H_wfc
    // get the gradient of energy with respect to the wfc, i.e., Wk_occNum_HamiltWfc
    add_psi(ParaV, kv, occ_number, H_wfc_TV, H_wfc_hartree, H_wfc_dft_XC, H_wfc_exx_XC, occNum_HamiltWfc, XC_func_rdmft, alpha_power);

    // get the gradient of energy with respect to the natural occupation numbers, i.e., Wk_occNum_wfcHamiltWfc
    add_occNum(*kv, occ_number, wfcHwfc_TV, wfcHwfc_hartree, wfcHwfc_dft_XC, wfcHwfc_exx_XC, occNum_wfcHamiltWfc, XC_func_rdmft, alpha_power);

    // get the total energy
    // add_wfcHwfc(kv->wk, occ_number, wfcHwfc_TV, wfcHwfc_hartree, wfcHwfc_XC, Etotal_n_k, XC_func_rdmft, alpha_power);
    // add_wfcHwfc(wg, wk_fun_occNum, wfcHwfc_TV, wfcHwfc_hartree, wfcHwfc_XC, Etotal_n_k, XC_func_rdmft, alpha_power);
    // E_RDMFT[3] = getEnergy(Etotal_n_k);
    // Parallel_Reduce::reduce_all(E_RDMFT[3]);

    return E_RDMFT[3];

    /****** get occNum_wfcHamiltWfc, occNum_HamiltWfc and Etotal ******/

}


// cal_type = 2 just support XC-functional without exx
template <typename TK, typename TR>
void RDMFT<TK, TR>::cal_Energy(const int cal_type)
{
    double E_Ewald = pelec->f_en.ewald_energy;
    double E_entropy = pelec->f_en.demet;
    double E_descf = pelec->f_en.descf = 0.0;
    // double E_descf = 0.0;
    double E_xc_KS = pelec->f_en.etxc - pelec->f_en.etxcc;
    double E_exx_KS = pelec->f_en.exx;
    double E_deband_KS = pelec->f_en.deband;
    double E_deband_harris_KS = pelec->f_en.deband_harris;

    if( cal_type == 1 )
    {
        // for E_TV
        ModuleBase::matrix ETV_n_k(wg.nr, wg.nc, true);
        occNum_Mul_wfcHwfc(wg, wfcHwfc_TV, ETV_n_k, 0);
        E_RDMFT[0] = getEnergy(ETV_n_k);

        // for Ehartree
        ModuleBase::matrix Ehartree_n_k(wg.nr, wg.nc, true);
        occNum_Mul_wfcHwfc(wg, wfcHwfc_hartree, Ehartree_n_k, 1);
        E_RDMFT[1] = getEnergy(Ehartree_n_k);

        // for Exc
        if( GlobalC::exx_info.info_global.cal_exx )
        {
            // E_RDMFT[2] = 0.0;
            ModuleBase::matrix Exc_n_k(wg.nr, wg.nc, true);
            // because we have got wk_fun_occNum, we can use symbol=1 realize it
            occNum_Mul_wfcHwfc(wk_fun_occNum, wfcHwfc_exx_XC, Exc_n_k, 1);
            E_RDMFT[2] = getEnergy(Exc_n_k);
            Parallel_Reduce::reduce_all(E_RDMFT[2]);

            // // test
            std::cout << "\n\n\n******\nE_exx-type in rdmft: " << E_RDMFT[2] << "\n******\n\n" << std::endl;
            std::cout << "\n\n\n******\nE_dft-xc in rdmft: " << etxc << "\n******\n\n" << std::endl;

            // if E_XC is hybrid functional
            E_RDMFT[2] += etxc;
        }
        else
        {
            E_RDMFT[2] = etxc;
        }

        // add up the results obtained by all processors, or we can do reduce_all(wfcHwfc_) before add_wg() used for Etotal to replace it
        Parallel_Reduce::reduce_all(E_RDMFT[0]);
        Parallel_Reduce::reduce_all(E_RDMFT[1]);

        Etotal = E_RDMFT[0] + E_RDMFT[1] + E_RDMFT[2] + E_Ewald + E_entropy + E_descf;

        // temp
        E_RDMFT[3] = E_RDMFT[0] + E_RDMFT[1] + E_RDMFT[2];
    }
    else
    {
        this->pelec->f_en.deband  = this->pelec->cal_delta_eband();
        E_descf = pelec->f_en.descf = 0.0;
        this->pelec->cal_energies(2);
        Etotal = this->pelec->f_en.etot;

        // if( GlobalC::exx_info.info_global.cal_exx )
        // {
        //     ModuleBase::matrix Exc_n_k(wg.nr, wg.nc, true);
        //     // because we have got wk_fun_occNum, we can use symbol=1 realize it
        //     occNum_Mul_wfcHwfc(wk_fun_occNum, wfcHwfc_XC, Exc_n_k, 1);
        //     E_RDMFT[2] = getEnergy(Exc_n_k);
        //     Parallel_Reduce::reduce_all(E_RDMFT[2]);

        //     // test
        //     Etotal -= E_RDMFT[2];
        // }
    }

    // print results
    if( GlobalC::exx_info.info_global.cal_exx )
    {
    std::cout << "\n\nfrom class RDMFT: \nXC_fun: " << XC_func_rdmft << ",   alpha_power:" << alpha_power << std::endl;

    std::cout << std::fixed << std::setprecision(10) << "******\nE(TV + Hartree + XC) by RDMFT:   " << E_RDMFT[3] << "\n\nETV_RDMFT:      " 
                << E_RDMFT[0] << "\nEhartree_RDMFT: " << E_RDMFT[1] << "\nExc_RDMFT:      " << E_RDMFT[2] << "\nE_Ewald:        " << E_Ewald
                << "\nE_entropy(-TS): " << E_entropy << "\nE_descf:        " << E_descf << "\n\nEtotal_RDMFT:   " << Etotal << "\nExc_dft:         " << E_xc_KS 
                << "\nE_exx_KS:         " << E_exx_KS <<"\n******\n\n" << std::endl;

    std::cout << "\nE_deband_KS:  " << E_deband_KS << "\nE_deband_harris_KS:  " << E_deband_harris_KS << "\n\n" << std::endl;
    }
    else
    {
        std::cout << "\n\nfrom class RDMFT: \nXC_fun: " << XC_func_rdmft << std::endl;

        std::cout << std::fixed << std::setprecision(10) << "******\nE(TV + Hartree + XC) by RDMFT:   " << E_RDMFT[3] << "\n\nETV_RDMFT:      " 
                    << E_RDMFT[0] << "\nE_hartree_RDMFT: " << E_RDMFT[1] << "\nEex_PBE_RDMFT:   " << E_RDMFT[2] << "\nE_Ewald:        " << E_Ewald
                    << "\nE_entropy(-TS): " << E_entropy << "\nE_descf:        " << E_descf << "\n\nEtotal_RDMFT:   " << Etotal << "\nExc_dft:         " << E_xc_KS 
                    << "\nE_exx_KS:         " << E_exx_KS <<"\n******\n\n" << std::endl;

        std::cout << "\netxc:  " << etxc << "\nvtxc:  " << vtxc << "\n";

        std::cout << "\nE_deband_KS:  " << E_deband_KS << "\nE_deband_harris_KS:  " << E_deband_harris_KS << "\n\n" << std::endl;
    }

    ModuleBase::timer::tick("rdmftTest", "RDMFT_E&Egradient");

}


template <typename TK, typename TR>
double RDMFT<TK, TR>::Run(ModuleBase::matrix& E_gradient_occNum, psi::Psi<TK>& E_gradient_wfc)
{
    // this->cal_V_hartree();
    // this->cal_V_XC();
    // this->cal_Hk_Hpsi();
    this->cal_E_gradient();
    this->cal_Energy(this->cal_E_type);

    E_gradient_occNum = (occNum_wfcHamiltWfc);
    
    TK* pwfc = &occNum_HamiltWfc(0, 0, 0);
    TK* pwfc_out = &E_gradient_wfc(0, 0, 0);
    for(int i=0; i<wfc.size(); ++i) pwfc_out[i] = pwfc[i];

    // test
    // rdmft::printMatrix_pointer(E_gradient_occNum.nr, E_gradient_occNum.nc, &E_gradient_occNum(0, 0), "E_gradient_occNum");
    // rdmft::printMatrix_pointer(occ_number.nr, occ_number.nc, &occ_number(0, 0), "occ_number");
    // rdmft::printMatrix_pointer(wfcHwfc_TV.nr, wfcHwfc_TV.nc, &wfcHwfc_TV(0, 0), "wfcHwfc_TV");
    // rdmft::printMatrix_pointer(wfcHwfc_hartree.nr, wfcHwfc_hartree.nc, &wfcHwfc_hartree(0, 0), "wfcHwfc_hartree");
    // rdmft::printMatrix_pointer(wfcHwfc_XC.nr, wfcHwfc_XC.nc, &wfcHwfc_XC(0, 0), "wfcHwfc_XC");
    // rdmft::printMatrix_pointer(E_gradient_wfc.get_nbands(), E_gradient_wfc.get_nbasis(), &E_gradient_wfc(0, 0, 0), "E_gradient_wfc(ik=0)");
    // rdmft::printMatrix_pointer(E_gradient_wfc.get_nbands(), E_gradient_wfc.get_nbasis(), &E_gradient_wfc(2, 0, 0), "E_gradient_wfc(ik=2)");
    // test

    // // test
    // std::cout << "\n\n rdmft_solver->wg(ik, inband): " << std::endl;
    // for(int ik=0; ik < this->wg.nr; ++ik)
    // {
    //     for(int inb=0; inb < this->wg.nc; ++inb) std::cout << this->wg(ik, inb) << " ";
    //     std::cout << "\n" << std::endl;
    // }

    // return E_RDMFT[3];
    return Etotal;
}



template class RDMFT<double, double>;
template class RDMFT<std::complex<double>, double>;
template class RDMFT<std::complex<double>, std::complex<double>>;



}


