//==========================================================
// Author: Jingang Han
// DATE : 2024-10-30
//==========================================================

#include "rdmft.h"
#include "module_rdmft/rdmft_tools.h"
#include "module_psi/psi.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"

#ifdef __EXX
#include "module_ri/RI_2D_Comm.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
#endif
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/ekinetic_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/nonlocal_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/veff_lcao.h"

namespace rdmft
{

//! @brief This file is to get each potential matrix in NAOs, 'V' represents the potential matrix


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
            nspin,
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
            nspin,
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
            nspin,
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
            nspin,
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
void RDMFT<TK, TR>::cal_V_XC(const UnitCell& ucell)
{
    // // //test
    // DM_XC_pass = DM_XC;

    // elecstate::DensityMatrix<TK, double> DM_test(ParaV, nspin, kv->kvec_d, nk_total);
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

    if( !only_exx_type )
    {
        HR_dft_XC->set_zero();
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
                nspin,
                charge,
                rho_basis,
                vloc,
                sf,
                "xc",
                &etxc,
                &vtxc
            );
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
                nspin,
                charge,
                rho_basis,
                vloc,
                sf,
                "xc",
                &etxc,
                &vtxc
            );
        }
        V_dft_XC->contributeHR();
    }

#ifdef __EXX
    if(GlobalC::exx_info.info_global.cal_exx)
    {
        HR_exx_XC->set_zero();

        std::vector< std::vector<TK> > DM_XC(nk_total, std::vector<TK>(ParaV->nloc));
        get_DM_XC(DM_XC);
        // get DM_XC of all k points
        if( exx_spacegroup_symmetry )
        {
            DM_XC = symrot_exx.restore_dm(*this->kv, DM_XC, *ParaV); // class vector could be auto resize()
        }
        std::vector< const std::vector<TK>* > DM_XC_pointer(DM_XC.size());
        for(int ik=0; ik<DM_XC.size(); ++ik) { DM_XC_pointer[ik] = &DM_XC[ik];
}

        if (GlobalC::exx_info.info_ri.real_number)
        {
            // transfer the DM_XC to appropriate format
            std::vector<std::map<int,std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>>> 
                Ds_XC_d = std::is_same<TK, double>::value //gamma_only_local
                ? RI_2D_Comm::split_m2D_ktoR<double>(ucell,*kv, DM_XC_pointer, *ParaV, nspin)
                : RI_2D_Comm::split_m2D_ktoR<double>(ucell,*kv, DM_XC_pointer, *ParaV, nspin, this->exx_spacegroup_symmetry);

            // provide the Ds_XC to Vxc_fromRI(V_exx_XC)
            if (this->exx_spacegroup_symmetry && GlobalC::exx_info.info_global.exx_symmetry_realspace)
            {
                Vxc_fromRI_d->cal_exx_elec(Ds_XC_d, ucell,*ParaV, &this->symrot_exx);
            }
            else
            {
                Vxc_fromRI_d->cal_exx_elec(Ds_XC_d, ucell,*ParaV);
            }

            // when we doing V_exx_XC.contributeHk(ik), we get HK_XC constructed by the special DM_XC
            V_exx_XC = new hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>>(
                hsk_exx_XC,
                HR_exx_XC,
                ucell,
                *kv,
                &Vxc_fromRI_d->Hexxs,
                nullptr,
                hamilt::Add_Hexx_Type::k
            );
        }
        else
        {
            // transfer the DM_XC to appropriate format
            std::vector<std::map<int,std::map<std::pair<int,std::array<int,3>>,RI::Tensor<std::complex<double>>>>> 
                Ds_XC_c = std::is_same<TK, double>::value //gamma_only_local
                ? RI_2D_Comm::split_m2D_ktoR<std::complex<double>>(ucell,*kv, DM_XC_pointer, *ParaV, nspin)
                : RI_2D_Comm::split_m2D_ktoR<std::complex<double>>(ucell,*kv, DM_XC_pointer, *ParaV, nspin, this->exx_spacegroup_symmetry);

            // // provide the Ds_XC to Vxc_fromRI(V_exx_XC)
            if (this->exx_spacegroup_symmetry && GlobalC::exx_info.info_global.exx_symmetry_realspace)
            {
                Vxc_fromRI_c->cal_exx_elec(Ds_XC_c, ucell,*ParaV, &this->symrot_exx);
            }
            else
            {
                Vxc_fromRI_c->cal_exx_elec(Ds_XC_c, ucell,*ParaV);
            }

            // when we doing V_exx_XC.contributeHk(ik), we get HK_XC constructed by the special DM_XC
            V_exx_XC = new hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>>(
                hsk_exx_XC,
                HR_exx_XC,
                ucell,
                *kv,
                nullptr,
                &Vxc_fromRI_c->Hexxs,
                hamilt::Add_Hexx_Type::k
            );
        }
        // use hamilt::Add_Hexx_Type::k, not R, contributeHR() should be skipped
        // V_exx_XC->contributeHR();
    }
#endif

}




template class RDMFT<double, double>;
template class RDMFT<std::complex<double>, double>;
template class RDMFT<std::complex<double>, std::complex<double>>;

}
