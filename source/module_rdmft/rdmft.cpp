//==========================================================
// Author: Jingang Han
// DATE : 2024-03-11
//==========================================================

#include "rdmft.h"
#include "module_rdmft/rdmft_tools.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/parallel_reduce.h"
#include "module_cell/module_symmetry/symmetry.h"


#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <algorithm>



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
    delete HR_exx_XC;
    delete HR_dft_XC;
    // delete HR_local;
    delete hsk_TV;
    delete hsk_hartree;
    delete hsk_dft_XC;
    delete hsk_exx_XC;
#ifdef __EXX
    delete Vxc_fromRI_d;
    delete Vxc_fromRI_c;
#endif
    delete V_ekinetic_potential;
    delete V_nonlocal;
    delete V_local;
    delete V_hartree;
    delete V_exx_XC;
    delete V_dft_XC;
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

    nspin = PARAM.inp.nspin;
    nbands_total = PARAM.inp.nbands;
    nk_total = ModuleSymmetry::Symmetry::symm_flag == -1 ? kv->get_nkstot_full(): kv->get_nks();
    nk_total *= nspin;
    only_exx_type = ( XC_func_rdmft == "hf" || XC_func_rdmft == "muller" || XC_func_rdmft == "power" );

    // create desc[] and something about MPI to Eij(nbands*nbands)
#ifdef __MPI
    para_Eij.set(nbands_total, nbands_total, ParaV->nb, ParaV->blacs_ctxt); // maybe in default, PARAM.inp.nb2d = 0, can't be used
#endif

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
    hsk_TV = new hamilt::HS_Matrix_K<TK>(ParaV, true);
    hsk_hartree = new hamilt::HS_Matrix_K<TK>(ParaV, true);
    hsk_dft_XC = new hamilt::HS_Matrix_K<TK>(ParaV, true);
    hsk_exx_XC = new hamilt::HS_Matrix_K<TK>(ParaV, true);

    HK_XC.resize( ParaV->get_row_size()*ParaV->get_col_size() );
    // HK_RDMFT_pass.resize(nk_total, ParaV->get_row_size(), ParaV->get_col_size());
    // HK_XC_pass.resize(nk_total, ParaV->get_row_size(), ParaV->get_col_size());


    Eij_TV.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );
    Eij_hartree.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );
    Eij_XC.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );
    Eij_exx_XC.resize( para_Eij.get_row_size()*para_Eij.get_col_size() );

    // 
    HR_TV = new hamilt::HContainer<TR>(*ucell, ParaV);
    HR_hartree = new hamilt::HContainer<TR>(*ucell, ParaV);
    HR_exx_XC = new hamilt::HContainer<TR>(*ucell, ParaV);
    HR_dft_XC = new hamilt::HContainer<TR>(*ucell, ParaV);
    // HR_local = new hamilt::HContainer<TR>(*ucell, ParaV);

    wfc.zero_out();
    occNum_HamiltWfc.zero_out();
    H_wfc_TV.zero_out();
    H_wfc_hartree.zero_out();
    H_wfc_XC.zero_out();
    H_wfc_exx_XC.zero_out();
    H_wfc_dft_XC.zero_out();
    
    HR_TV->set_zero();         // HR->set_zero() might be delete here, test on Gamma_only in the furure 
    HR_hartree->set_zero();
    HR_exx_XC->set_zero();
    HR_dft_XC->set_zero();
    // HR_local->set_zero();

#ifdef __EXX
    if( GlobalC::exx_info.info_global.cal_exx )
    {
        // if the irreducible k-points can change with symmetry during cell-relax, it should be moved back to update_ion()
        exx_spacegroup_symmetry = (PARAM.inp.nspin < 4 && ModuleSymmetry::Symmetry::symm_flag == 1);
        if (exx_spacegroup_symmetry)
        {
            const std::array<int, 3>& period = RI_Util::get_Born_vonKarmen_period(*kv);
            this->symrot_exx.find_irreducible_sector(ucell->symm, ucell->atoms, ucell->st,
                    RI_Util::get_Born_von_Karmen_cells(period), period, ucell->lat);
            this->symrot_exx.cal_Ms(*kv, *ucell, *ParaV);
        }

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
#endif

    if( PARAM.inp.gamma_only )
    {
        HR_TV->fix_gamma();
        HR_hartree->fix_gamma();
        HR_exx_XC->fix_gamma();
        HR_dft_XC->fix_gamma();
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
        std::fill(HK_XC.begin(), HK_XC.end(), 0.0);

        // get the HK with ik-th k vector, the result is stored in HK_TV, HK_hartree and HK_XC respectively
        V_local->contributeHk(ik);
        V_hartree->contributeHk(ik);

        // get H(k) * wfc
        HkPsi( ParaV, hsk_TV->get_hk()[0], wfc(ik, 0, 0), H_wfc_TV(ik, 0, 0));
        HkPsi( ParaV, hsk_hartree->get_hk()[0], wfc(ik, 0, 0), H_wfc_hartree(ik, 0, 0));

        // get wfc * H(k)_wfc
        cal_bra_op_ket( ParaV, para_Eij, wfc(ik, 0, 0), H_wfc_TV(ik, 0, 0), Eij_TV );
        cal_bra_op_ket( ParaV, para_Eij, wfc(ik, 0, 0), H_wfc_hartree(ik, 0, 0), Eij_hartree );
        _diagonal_in_serial( para_Eij, Eij_TV, &(wfcHwfc_TV(ik, 0)) );
        _diagonal_in_serial( para_Eij, Eij_hartree, &(wfcHwfc_hartree(ik, 0)) );

#ifdef __EXX
        if(GlobalC::exx_info.info_global.cal_exx)
        {
            hsk_exx_XC->set_zero_hk();

            V_exx_XC->contributeHk(ik);
            HkPsi( ParaV, hsk_exx_XC->get_hk()[0], wfc(ik, 0, 0), H_wfc_exx_XC(ik, 0, 0));
            cal_bra_op_ket( ParaV, para_Eij, wfc(ik, 0, 0), H_wfc_exx_XC(ik, 0, 0), Eij_exx_XC );
            _diagonal_in_serial( para_Eij, Eij_exx_XC, &(wfcHwfc_exx_XC(ik, 0)) );
            
            for(int iloc=0; iloc<HK_XC.size(); ++iloc) HK_XC[iloc] += hsk_exx_XC->get_hk()[iloc];
        }
#endif
        if( !only_exx_type )
        {
            hsk_dft_XC->set_zero_hk();

            V_dft_XC->contributeHk(ik);
            HkPsi( ParaV, hsk_dft_XC->get_hk()[0], wfc(ik, 0, 0), H_wfc_dft_XC(ik, 0, 0));
            cal_bra_op_ket( ParaV, para_Eij, wfc(ik, 0, 0), H_wfc_dft_XC(ik, 0, 0), Eij_XC );
            _diagonal_in_serial( para_Eij, Eij_XC, &(wfcHwfc_dft_XC(ik, 0)) );
            
            for(int iloc=0; iloc<HK_XC.size(); ++iloc) HK_XC[iloc] += hsk_dft_XC->get_hk()[iloc];
        }

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

    }

    // std::cout << "\n\nsum of XC_minus_XC: " << XC_minus_XC << "\n\n" << std::endl;

}


template <typename TK, typename TR>
double RDMFT<TK, TR>::cal_E_grad_wfc_occ_num()
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

    double E_exxType_rdmft = 0.0; // delete in the future

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
        E_RDMFT[2] = 0.0;
#ifdef __EXX
        if( GlobalC::exx_info.info_global.cal_exx )
        {
            ModuleBase::matrix Exc_n_k(wg.nr, wg.nc, true);
            // because we have got wk_fun_occNum, we can use symbol=1 realize it
            occNum_Mul_wfcHwfc(wk_fun_occNum, wfcHwfc_exx_XC, Exc_n_k, 1);
            E_RDMFT[2] = getEnergy(Exc_n_k);
            Parallel_Reduce::reduce_all(E_RDMFT[2]);
            E_exxType_rdmft = E_RDMFT[2];
        }
#endif
        E_RDMFT[2] += etxc;

        // add up the results obtained by all processors, or we can do reduce_all(wfcHwfc_) before add_wg() used for Etotal to replace it
        Parallel_Reduce::reduce_all(E_RDMFT[0]);
        Parallel_Reduce::reduce_all(E_RDMFT[1]);

        this->Etotal = E_RDMFT[0] + E_RDMFT[1] + E_RDMFT[2] + E_Ewald + E_entropy + E_descf;

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

//     // print results
//     std::cout << "\n\nfrom class RDMFT: \nXC_fun: " << XC_func_rdmft << std::endl;
// #ifdef __EXX
//     if( GlobalC::exx_info.info_global.cal_exx ) std::cout << "alpha_power: " << alpha_power << std::endl;
// #endif
//     std::cout << std::fixed << std::setprecision(10) 
//                 << "******\nE(TV + Hartree + XC) by RDMFT:   " << E_RDMFT[3] 
//                 << "\n\nE_TV_RDMFT:      " << E_RDMFT[0] 
//                 << "\nE_hartree_RDMFT: " << E_RDMFT[1] 
//                 << "\nExc_" << XC_func_rdmft << "_RDMFT:    " << E_RDMFT[2] 
//                 << "\nE_Ewald:         " << E_Ewald
//                 << "\nE_entropy(-TS):  " << E_entropy 
//                 << "\nE_descf:         " << E_descf
//                 << "\n\nEtotal_RDMFT:    " << Etotal 
//                 << "\n\nExc_ksdft:       " << E_xc_KS 
//                 << "\nE_exx_ksdft:     " << E_exx_KS 
//                 <<"\n******\n\n" << std::endl;

//     std::cout << "\netxc:  " << etxc << "\nvtxc:  " << vtxc << "\n";
//     std::cout << "\nE_deband_KS:  " << E_deband_KS << "\nE_deband_harris_KS:  " << E_deband_harris_KS << "\n\n" << std::endl;

    if( PARAM.inp.rdmft == true )
    {
        GlobalV::ofs_running << "\n\nfrom class RDMFT: \nXC_fun: " << XC_func_rdmft << std::endl;
#ifdef __EXX
        if( GlobalC::exx_info.info_global.cal_exx ) { GlobalV::ofs_running << "alpha_power: " << alpha_power << std::endl;
}
#endif
        // GlobalV::ofs_running << std::setprecision(12);
        // GlobalV::ofs_running << std::setiosflags(std::ios::right);
        GlobalV::ofs_running << std::fixed << std::setprecision(10)
                << "\n******\nE(TV + Hartree + XC) by RDMFT:   " << E_RDMFT[3] 
                << "\n\nE_TV_RDMFT:      " << E_RDMFT[0] 
                << "\nE_hartree_RDMFT: " << E_RDMFT[1] 
                << "\nExc_" << XC_func_rdmft << "_RDMFT:    " << E_RDMFT[2] 
                << "\nE_Ewald:         " << E_Ewald
                << "\nE_entropy(-TS):  " << E_entropy 
                << "\nE_descf:         " << E_descf 
                << "\n\nEtotal_RDMFT:    " << Etotal 
                << "\n\nExc_ksdft:       " << E_xc_KS 
                << "\nE_exx_ksdft:     " << E_exx_KS
                << "\nE_exxType_rdmft: " << E_exxType_rdmft
                <<"\n******\n" << std::endl;
    }
    std::cout << std::defaultfloat;

}


template <typename TK, typename TR>
double RDMFT<TK, TR>::run(ModuleBase::matrix& E_gradient_occNum, psi::Psi<TK>& E_gradient_wfc)
{
    ModuleBase::TITLE("RDMFT", "E_Egradient");
    ModuleBase::timer::tick("RDMFT", "E_Egradient");

    // this->cal_V_hartree();
    // this->cal_V_XC();
    this->cal_Hk_Hpsi();
    this->cal_E_grad_wfc_occ_num();
    this->cal_Energy(this->cal_E_type);
    // this->cal_Energy(2);

    E_gradient_occNum = (occNum_wfcHamiltWfc);
    
    TK* pwfc = &occNum_HamiltWfc(0, 0, 0);
    TK* pwfc_out = &E_gradient_wfc(0, 0, 0);
    for(int i=0; i<wfc.size(); ++i) { pwfc_out[i] = pwfc[i]; }

    ModuleBase::timer::tick("RDMFT", "E_Egradient");
    // return E_RDMFT[3];
    return Etotal;
}

template class RDMFT<double, double>;
template class RDMFT<std::complex<double>, double>;
template class RDMFT<std::complex<double>, std::complex<double>>;

}


