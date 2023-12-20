#ifndef RDMFT_TEST_H
#define RDMFT_TEST_H

#include "module_base/matrix.h"
//#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_base/blas_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_basis/module_ao/parallel_2d.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#include "module_hamilt_general/operator.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/ekinetic_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/nonlocal_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/veff_lcao.h"

// used by class Veff_rdmft
#include "module_base/timer.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
//#include "module_hamilt_lcao/module_gint/gint_gamma.h"
//#include "module_hamilt_lcao/module_gint/gint_k.h"
//#include "operator_lcao.h"
//#include "module_cell/module_neighbor/sltk_grid_driver.h"
//#include "module_cell/unitcell.h"


#include <iostream>
#include <type_traits>
#include <complex>



namespace hamilt
{


// for test use dgemm_
void printResult_dgemm();

//for print matrix
template <typename TK>
void printMatrix_pointer(int M, int N, TK& matrixA, std::string nameA)
{
    std::cout << "\n" << nameA << ": \n";
    for(int i=0; i<M; ++i)
    {
        for(int j=0; j<N; ++j)
        {
            std::cout << *(&matrixA+i*N+j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


template <typename TK>
void printMatrix_vector(int M, int N, std::vector<TK>& matrixA, std::string nameA)
{
    std::cout << "\n" << nameA << ": \n";
    for(int i=0; i<M; ++i)
    {
        for(int j=0; j<N; ++j)
        {
            std::cout << matrixA[i*N+j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}




template <typename TK, typename TR>
class Veff_rdmft : public OperatorLCAO<TK, TR>
{
  public:
    Veff_rdmft(Gint_k* GK_in,
                      Local_Orbital_Charge* loc_in,
                      LCAO_Matrix* LM_in,
                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                      elecstate::Potential* pot_in,
                      hamilt::HContainer<TR>* hR_in,
                      std::vector<TK>* hK_in,
                      const UnitCell* ucell_in,
                      Grid_Driver* GridD_in,
                      const Parallel_Orbitals* paraV)
        : GK(GK_in),
          loc(loc_in),
          pot(pot_in),
          OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
    {
        this->cal_type = lcao_gint;

        this->initialize_HR(ucell_in, GridD_in, paraV);
        GK_in->initialize_pvpR(*ucell_in, GridD_in);
    }
    Veff_rdmft(Gint_Gamma* GG_in,
                          Local_Orbital_Charge* loc_in,
                          LCAO_Matrix* LM_in,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                          elecstate::Potential* pot_in,
                          hamilt::HContainer<TR>* hR_in,
                          std::vector<TK>* hK_in,
                          const UnitCell* ucell_in,
                          Grid_Driver* GridD_in,
                          const Parallel_Orbitals* paraV
                          )
        : GG(GG_in), loc(loc_in), pot(pot_in),
        OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
    {
        this->cal_type = lcao_gint;

        this->initialize_HR(ucell_in, GridD_in, paraV);

        GG_in->initialize_pvpR(*ucell_in, GridD_in);
    }

    ~Veff_rdmft(){};

    virtual void contributeHR() override;


  private:
    // used for k-dependent grid integration.
    Gint_k* GK = nullptr;

    // used for gamma only algorithms.
    Gint_Gamma* GG = nullptr;

    // Charge calculating method in LCAO base and contained grid base calculation: DM_R, DM, pvpR_reduced
    Local_Orbital_Charge* loc = nullptr;

    elecstate::Potential* pot = nullptr;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the electronic kinetic matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(const UnitCell* ucell_in, Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);

};









double wg_func(double wg, int symbol = 0);


template <typename TK>
void set_zero_HK(std::vector<TK>& HK)
{
    for(int i=0; i<HK.size(); ++i) HK[i] = 0.0;
}


// wfc and H_wfc need to be k_firest and provide wfc(ik, 0, 0) and H_wfc(ik, 0, 0)
// psi::Psi<TK> psi's (), std::vector<TK> HK's [] operator overloading return TK
template <typename TK>
void HkPsi(const Parallel_Orbitals* ParaV, const Parallel_2D& para_wfc_in, const TK& HK, const TK& wfc, TK& H_wfc)
{

    const int one_int = 1;
    //const double one_double = 1.0, zero_double = 0.0;
    const std::complex<double> one_complex = {1.0, 0.0};
    const std::complex<double> zero_complex = {0.0, 0.0};
    const char N_char = 'N';
    const char T_char = 'T';

    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    //test
    std::cout << "\n\n\n******\nH_wfc[0], Hwfc[1]: " << H_wfc << " " << *(&H_wfc+1) <<"\n";

    //because wfc(bands, basis'), H(basis, basis'), we do wfc*H^T(in the perspective of cpp, not in fortran). And get H_wfc(bands, basis) is correct.
    // pzgemm_( &N_char, &T_char, &nbands, &nbasis, &nbasis, &one_complex, wfc, &one_int, &one_int, ParaV->desc_wfc,
    //     HK, &one_int, &one_int, ParaV->desc, &zero_complex, H_wfc, &one_int, &one_int, ParaV->desc_wfc );
    // pzgemm_( &T_char, &N_char, &nbasis, &nbands, &nbasis, &one_complex, &HK, &one_int, &one_int, ParaV->desc,
    //     &wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_complex, &H_wfc, &one_int, &one_int, ParaV->desc_wfc );

    // zgemm_( &T_char, &N_char, &nbasis, &nbands, &nbasis, &one_complex, &HK, &nbasis, &wfc, &nbasis, &zero_complex, &H_wfc, &nbasis );

    for(int ib=0; ib<nbands; ++ib)
    {
        for(int ibs2=0; ibs2<nbasis; ++ibs2)
        {
            for(int ibs1=0; ibs1<nbasis; ++ibs1)
            {
                *(&H_wfc+ib*nbasis+ibs2) += *(&HK+ibs2*nbasis+ibs1) * *(&wfc+ib*nbasis+ibs1);
            }
        }
    }




    //test
    std::cout << "after HkPsi(),\nH_wfc[0], Hwfc[1]: " << H_wfc << " " << *(&H_wfc+1) <<"\n******\n\n\n";


    /*
    if ( std::is_same<TK, double>::value )
    {
        pdgemm_( &N_char, &T_char, &nbands, &nbasis, &nbasis, &one_double, &wfc, &one_int, &one_int, ParaV->desc_wfc,
                &HK, &one_int, &one_int, ParaV->desc, &zero_double, &H_wfc, &one_int, &one_int, ParaV->desc_wfc );
    }
    else if ( std::is_same<TK, std::complex<double>>::value )
    {
        pzgemm_( &N_char, &T_char, &nbands, &nbasis, &nbasis, &one_complex, &wfc, &one_int, &one_int, ParaV->desc_wfc,
                &HK, &one_int, &one_int, ParaV->desc, &zero_complex, &H_wfc, &one_int, &one_int, ParaV->desc_wfc );
    }
    else std::cout << "\n\n******\nthere maybe something wrong when calling rdmft_cal() and use pdemm_/pzgemm_\n******\n\n";
    */
}


template <>
void HkPsi<double>(const Parallel_Orbitals* ParaV, const Parallel_2D& para_wfc_in, const double& HK, const double& wfc, double& H_wfc);


// ModuleBase::matrix wfcHwfc(ik, 0)
template <typename TK>
void psiDotPsi(const Parallel_Orbitals* ParaV, const Parallel_2D& para_wfc_in, const Parallel_2D& para_Eij_in,
                const TK& wfc, const TK& H_wfc, std::vector<TK>& Dmn, double* wfcHwfc)
{
    const int one_int = 1;
    //const double one_double = 1.0, zero_double = 0.0;
    const std::complex<double> one_complex = {1.0, 0.0};
    const std::complex<double> zero_complex = {0.0, 0.0};
    const char N_char = 'N';
    const char C_char = 'C'; 

    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    // const int nrow_bands = para_Eij_in.get_row_size();
    // const int ncol_bands = para_Eij_in.get_col_size();
    const int nrow_bands = para_Eij_in.get_row_size();
    const int ncol_bands = para_Eij_in.get_col_size();

    //test
    std::cout << "\n\n\n******\nwfcHwfc[0], wfcHwfc[1]: " << wfcHwfc[0] << " " << wfcHwfc[1] <<"\n";
    
    // in parallel_orbitals.h, there has int desc_Eij[9] which used for Eij in TDDFT, nbands*nbands. Just proper here.
    // std::vector<TK> Dmn(nrow_bands*ncol_bands);
    // pzgemm_( &C_char, &N_char, &nbands, &nbands, &nbasis, &one_complex, &wfc, &one_int, &one_int, ParaV->desc_wfc,
    //         &H_wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_complex, &Dmn[0], &one_int, &one_int, para_Eij_in.desc );
    // zgemm_( &C_char, &N_char, &nbands, &nbands, &nbasis, &one_complex,  &wfc, &nbasis, &H_wfc, &nbasis, &zero_complex, &Dmn[0], &nbands );

    set_zero_HK(Dmn);
    for(int ib1=0; ib1<nbands; ++ib1)
    {
        for(int ib2=0; ib2<nbands; ++ib2)
        {
            for(int ibas=0; ibas<nbasis; ++ibas)
            {
                Dmn[ib1*nbands+ib2] += std::conj( *(&wfc+ib1*nbasis+ibas) ) * (*(&H_wfc+ib2*nbasis+ibas));
            }
        }
    }


    
    for(int i=0; i<nrow_bands; ++i)
    {
        int i_global = para_Eij_in.local2global_row(i);
        for(int j=0; j<ncol_bands; ++j)
        {
            int j_global = para_Eij_in.local2global_col(j);
            if(i_global==j_global)
            {
                wfcHwfc[j_global] = std::real( Dmn[i*ncol_bands+j] ); // need to be sure imag(Dmn[i*nrow_bands+j]) == 0
                // double testEnn = std::abs( std::imag( Dmn[i*ncol_bands+j] ) );
                // if( testEnn>1e-16 )std::cout << "\n\nimag(Enn)!=0? imag(Enn)= " << testEnn << "\n\n";
            }
        }
    }

}

template <>
void psiDotPsi<double>(const Parallel_Orbitals* ParaV, const Parallel_2D& para_wfc_in, const Parallel_2D& para_Eij_in,
                        const double& wfc, const double& H_wfc, std::vector<double>& Dmn, double* wfcHwfc);


// realize wg_wfc = wg * wfc. Calling this function and we can get wfc = wg*wfc.
template <typename TK>
void wgMulPsi(const Parallel_Orbitals* ParaV, const Parallel_2D& para_wfc_in, const ModuleBase::matrix& wg, psi::Psi<TK>& wfc, int symbol = 0)
{
    const int nk_local = wfc.get_nk();
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    for (int ik = 0; ik < nk_local; ++ik)
    {
        for (int ib_local = 0; ib_local < nbands_local; ++ib_local)
        {
            const double wg_local = wg_func( wg(ik, ParaV->local2global_col(ib_local)), symbol);  // JG：得到确定n，k下的占据数eta_nk(的函数)
            TK* wfc_pointer = &(wfc(ik, ib_local, 0));        // JG：这里k_first=true,将wfc确定k，n，再找到基组指标miu=0处的指针，对指针做加减法就可以得到不同原子基miu处的wfc元素
            BlasConnector::scal(nbasis_local, wg_local, wfc_pointer, 1);        // JG：确定k，n的wg_wfc有nbasis_local个不同的原子基（miu），将它们都与确定n，k下的占据数wg_local相乘
        }
    }
}


// add psi with eta and g(eta)
template <typename TK>
void add_psi(const Parallel_Orbitals* ParaV, const Parallel_2D& para_wfc_in, const ModuleBase::matrix& wg, 
                psi::Psi<TK>& psi_TV, psi::Psi<TK>& psi_hartree, psi::Psi<TK>& psi_XC, psi::Psi<TK>& wg_Hpsi)
{
    const int nk = psi_TV.get_nk();
    const int nbn_local = psi_TV.get_nbands();
    const int nbs_local = psi_TV.get_nbasis();
    wgMulPsi(ParaV, para_wfc_in, wg, psi_TV);
    wgMulPsi(ParaV, para_wfc_in, wg, psi_hartree);
    wgMulPsi(ParaV, para_wfc_in, wg, psi_XC, 2);

    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    for(int ik=0; ik<nk; ++ik)
    {
        for(int inbn=0; inbn<nbn_local; ++inbn)
        {
            TK* pwg_Hpsi = &( wg_Hpsi(ik, inbn, 0) );
            for(int inbs=0; inbs<nbs_local; ++inbs)
            {
                pwg_Hpsi[inbs] = psi_TV(ik, inbn, inbs) + psi_hartree(ik, inbn, inbs) + psi_XC(ik, inbn, inbs);
            }
        }
    }

}


// wg_wfcHwfc = wg*wfcHwfc + wg_wfcHwfc
// When symbol = 0, 1, 2, 3, 4, wg = wg, 0.5*wg, g(wg), 0.5*g(wg), d_g(wg)/d_ewg respectively. Default symbol=0.
void wgMul_wfcHwfc(const ModuleBase::matrix& wg, const ModuleBase::matrix& wfcHwfc, ModuleBase::matrix& wg_wfcHwfc, int symbol = 0);


// Default symbol = 0 for the gradient of Etotal with respect to occupancy
// symbol = 1 for the relevant calculation of Etotal
void add_wg(const ModuleBase::matrix& wg, const ModuleBase::matrix& wfcHwfc_TV_in, const ModuleBase::matrix& wfcHwfc_hartree_in,
                const ModuleBase::matrix& wfcHwfc_XC_in, ModuleBase::matrix& wg_wfcHwfc, int symbol = 0);


//give certain wg_wfcHwfc, get the corresponding energy
double sumWg_getEnergy(const ModuleBase::matrix& wg_wfcHwfc);


// for test add a function and call it in module_elecstate/elecstate_lcao.cpp
// !!!just used for k-dependent grid integration. For gamma only algorithms, transfer Gint_k& GK_in to Gint_Gamma GG_in and use it in Veff<OperatorLCAO<TK, TR>>
template <typename TK, typename TR>
double rdmft_cal(LCAO_Matrix* LM_in,
                        Parallel_Orbitals* ParaV,
                        const ModuleBase::matrix& wg,
                        const psi::Psi<TK>& wfc,
                        ModuleBase::matrix& wg_wfcHamiltWfc,
                        psi::Psi<TK>& wg_HamiltWfc,
                        const K_Vectors& kv_in,
                        Gint_k& GK_in,
                        Local_Orbital_Charge& loc_in,
                        const Charge& chg,
                        elecstate::Potential& pot_in)  // delete pot_in parameter later
{
    ModuleBase::TITLE("hamilt_lcao", "RDMFT_E&Egradient");
    ModuleBase::timer::tick("hamilt_lcao", "RDMFT_E&Egradient");

    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    
    // create desc[] and something about MPI to wfc(nbands*nbasis)
    // para_wfc.desc[2] describe row of global matrix(NBANDS here),para_wfc.desc[3] describe col of global matrix (NLOCAL here)
    Parallel_2D para_wfc;
    para_wfc.set_block_size(GlobalV::NB2D);
    para_wfc.set_proc_dim(GlobalV::DSIZE);
    para_wfc.comm_2D = ParaV->comm_2D;
    para_wfc.blacs_ctxt = ParaV->blacs_ctxt;
    para_wfc.set_local2global( GlobalV::NLOCAL, GlobalV::NBANDS, ofs_running, ofs_warning );
    para_wfc.set_desc( GlobalV::NLOCAL, GlobalV::NBANDS, para_wfc.get_row_size(), false );

    // create desc[] and something about MPI to Eij(nbands*nbands)
    Parallel_2D para_Eij;
    para_Eij.set_block_size(GlobalV::NB2D);
    para_Eij.set_proc_dim(GlobalV::DSIZE);
    para_Eij.comm_2D = ParaV->comm_2D;
    para_Eij.blacs_ctxt = ParaV->blacs_ctxt;
    para_Eij.set_local2global( GlobalV::NBANDS, GlobalV::NBANDS, ofs_running, ofs_warning );
    para_Eij.set_desc( GlobalV::NBANDS, GlobalV::NBANDS, para_Eij.get_row_size(), false );

    // initialization
    const int nk_total = wfc.get_nk();
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();
    const std::vector<ModuleBase::Vector3<double>> kvec_d_in = kv_in.kvec_d;

    //hK_in nk*nbasis*nbasis
    hamilt::HContainer<TR> HR_TV(GlobalC::ucell, ParaV);
    hamilt::HContainer<TR> HR_hartree(GlobalC::ucell, ParaV);
    hamilt::HContainer<TR> HR_XC(GlobalC::ucell, ParaV);
    std::vector<TK> HK_TV(ParaV->get_row_size()*ParaV->get_col_size());
    std::vector<TK> HK_hartree(ParaV->get_row_size()*ParaV->get_col_size());
    std::vector<TK> HK_XC(ParaV->get_row_size()*ParaV->get_col_size());
    
    //set zero ( std::vector will automatically be set to zero )
    HR_TV.set_zero();
    HR_hartree.set_zero();
    HR_XC.set_zero();

    // get every Hamiltion matrix

    OperatorLCAO<TK, TR>* V_ekinetic_potential = new EkineticNew<OperatorLCAO<TK, TR>>(
        LM_in,
        kvec_d_in,
        &HR_TV,
        &HK_TV,
        &GlobalC::ucell,
        &GlobalC::GridD,
        ParaV
    );
    
    OperatorLCAO<TK, TR>* V_nonlocal = new NonlocalNew<OperatorLCAO<TK, TR>>(
        LM_in,
        kvec_d_in,
        &HR_TV,
        &HK_TV,
        &GlobalC::ucell,
        &GlobalC::GridD,
        ParaV
    );

    // use class Veff<OperatorLCAO<TK, TR>> get the local potential
    std::vector<std::string> pot_register_localV;
    pot_register_localV.push_back("local");
    pot_in.pot_register(pot_register_localV);
    OperatorLCAO<TK, TR>* V_local = new Veff<OperatorLCAO<TK, TR>>(
        &GK_in,
        &loc_in,
        LM_in,
        kvec_d_in,
        &pot_in,
        &HR_TV,
        &HK_TV,
        &GlobalC::ucell,
        &GlobalC::GridD,
        ParaV
    );
    // OperatorLCAO<TK, TR>* V_local_new = new OperatorLCAO<TK, TR>(
    //     LM_in,
    //     kvec_d_in,
    //     &HR_TV,
    //     &HK_TV
    // );

    // use class Veff<OperatorLCAO<TK, TR>> get the hartree potential
    std::vector<std::string> pot_register_hartree;
    pot_register_hartree.push_back("hartree");
    pot_in.pot_register(pot_register_hartree);
    OperatorLCAO<TK, TR>* V_hartree = new Veff<OperatorLCAO<TK, TR>>(
        &GK_in,
        &loc_in,
        LM_in,
        kvec_d_in,
        &pot_in,
        &HR_hartree,
        &HK_hartree,
        &GlobalC::ucell,
        &GlobalC::GridD,
        ParaV
    );
    OperatorLCAO<TK, TR>* V_hartree_new = new OperatorLCAO<TK, TR>(
        LM_in,
        kvec_d_in,
        &HR_hartree,
        &HK_hartree
    );


    OperatorLCAO<TK, TR>* V_XC = new OperatorEXX<OperatorLCAO<TK, TR>>(
        LM_in,
        &HR_XC,
        &HK_XC,
        kv_in
    );

    // // // add V_ekinetic, V_nonlocal and V_local, so we get V_ekinetic_potential = T + V_nonlocal + V_local
    // V_ekinetic_potential->add(V_nonlocal);
    // V_ekinetic_potential->add(V_local);

    // now HR_TV has the HR of V_ekinetic + V_nonlcao + V_local, 
    V_ekinetic_potential->contributeHR();
    V_nonlocal->contributeHR();
    V_local->contributeHR();

    // HR_hartree has the HR of V_hartree. HR_XC get from another way, so don't need to do this 
    V_hartree->contributeHR();

    //************test**************//
    TR* pHR_TV = HR_TV.data(0, 0);
    TR* pHR_hartree = HR_hartree.data(1, 0);
    TR* pHR_XC = HR_XC.data(1, 0);

    std::cout << "\n\n\n******\n";
    for(int i=0; i<10; ++i)
    {
        std::cout << "HR_TV, HR_hartree, HR_XC: " << pHR_TV[i] << " " << pHR_hartree[i] << " " << pHR_XC[i] << "\n";
    }
    std::cout << "******\n\n\n";
    //************test**************//

    //prepare for actual calculation
    //wg is global matrix, wg.nr*wg.nc = nk_total*nbands_global
    // ModuleBase::matrix wg_forEtotal(wg.nr, wg.nc, true);
    // ModuleBase::matrix wfcHwfc_TV(wg.nr, wg.nc, true);
    // ModuleBase::matrix wfcHwfc_hartree(wg.nr, wg.nc, true);
    // ModuleBase::matrix wfcHwfc_XC(wg.nr, wg.nc, true);

    ModuleBase::matrix wg_forEtotal(nk_total, GlobalV::NBANDS, true);
    ModuleBase::matrix wfcHwfc_TV(nk_total, GlobalV::NBANDS, true);
    ModuleBase::matrix wfcHwfc_hartree(nk_total, GlobalV::NBANDS, true);
    ModuleBase::matrix wfcHwfc_XC(nk_total, GlobalV::NBANDS, true);

    // let the 2d-block of H_wfc is same to wfc, so we can use desc_wfc and 2d-block messages of wfc to describe H_wfc
    psi::Psi<TK> H_wfc_TV(nk_total, nbands_local, nbasis_local);
    psi::Psi<TK> H_wfc_hartree(nk_total, nbands_local, nbasis_local);
    psi::Psi<TK> H_wfc_XC(nk_total, nbands_local, nbasis_local);

    // set zero
    TK* pH_wfc_TV = H_wfc_TV.get_pointer();
    TK* pH_wfc_hartree = H_wfc_hartree.get_pointer();
    TK* pH_wfc_XC = H_wfc_XC.get_pointer();
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024)
    #endif
    for(int i=0; i<H_wfc_TV.size(); ++i)
    {
        pH_wfc_TV[i] = 0.0;
        pH_wfc_hartree[i] = 0.0;
        pH_wfc_XC[i] = 0.0;
    }

    //************test**************//
    int wg_half_number = 0;
    int wg_zeroEight_num = 0;
    int wg_zeroThree_num = 0;
    for(int ik=0; ik<wg.nr; ++ik)
    {
        for(int ib=0; ib<wg.nc; ++ib)
        {
            if(ib%4 == 0) std::cout << "\n\n\n******\nwg: " << wg(ik, ib) << "\n******\n\n\n";
            if( wg(ik, ib) < 0 ) std::cout << "\n\n\n******\nsomething wrong in wg!!!\n******\n\n\n";
            if( wg(ik, ib) < 0.8 ) ++wg_zeroEight_num;
            if( wg(ik, ib) < 0.5 ) ++wg_half_number;
            if( wg(ik, ib) < 0.3 ) ++wg_zeroThree_num;
        }
    }
    std::cout << "\n\n\n******\nwg_size, wg<0.8, wg<0.5, wg<0.3: " << wg.nr*wg.nc << " " << wg_zeroEight_num << " " 
                << wg_half_number << " " << wg_zeroThree_num << "\n******\n\n\n";

    // int H_col_size = ParaV->get_col_size();
    // int wfc_col_size = para_wfc.get_col_size();

    // std::cout << 


    //************test**************//

    /*
    //test desc_wfc and desc_Eij
    std::cout << "\n\n\n******\n";
    for(int i=0 ; i<9; ++i)
    {
        std::cout << "Parav->desc_wfc[" << i << "], " << "para_wfc.desc[" << i << "]: " << ParaV->desc_wfc[i] << " " << para_wfc.desc[i] << "\n";
    }
    std::cout << "\n******\n\n\n";
    std::cout << "\n\n\n******\n";
    for(int i=0 ; i<9; ++i)
    {
        std::cout << "Parav->desc_Eij[" << i << "], " << "para_Eij.desc[" << i << "]: " << ParaV->desc_Eij[i] << " " << para_Eij.desc[i] << "\n";
    }
    std::cout << "\n******\n\n\n";
    */

    // just for temperate. in the future when realize psiDotPsi() without pzgemm_/pdgemm_,we don't need it
    // const int nrow_bands = ParaV->nrow_bands;
    // const int ncol_bands = ParaV->ncol_bands;
    // std::vector<TK> Eij(ncol_bands*ncol_bands);   //////////ncol*ncol? why
    const int nrow_bands = para_Eij.get_row_size();
    const int ncol_bands = para_Eij.get_col_size();
    std::vector<TK> Eij(nrow_bands*ncol_bands);

    //calculate wg_wfcHamiltWfc, wg_HamiltWfc and Etotal
    for(int ik=0; ik<nk_total; ++ik)
    {
        // get the HK with ik-th k vector
        V_ekinetic_potential->contributeHk(ik);
        // V_nonlocal->contributeHk(ik);
        // V_local_new->contributeHk(ik);
        V_hartree_new->contributeHk(ik);    // because contributeHk() in class Veff is {}, so we get a OperatorLCAO* class object V_hartree_new to do contributeHk()
        V_XC->contributeHk(ik);

        std::cout << "\n\nik= " << ik << "\n\n";

        if(ik==3)
        {
            std::cout << "\n\n\n******\n";
            for(int i=0; i<nbasis_local; ++i)
            {
                std::cout << "HK_TV, HK_hartree, HK_XC: " << HK_TV[i] << " " << HK_hartree[i] << " " << HK_XC[i] << "\n";
            }
            std::cout << "******\n\n\n";
        }

        //
        HkPsi( ParaV, para_wfc, HK_TV[0], wfc(ik, 0, 0), H_wfc_TV(ik, 0, 0));
        HkPsi( ParaV, para_wfc, HK_hartree[0], wfc(ik, 0, 0), H_wfc_hartree(ik, 0, 0));
        HkPsi( ParaV, para_wfc, HK_XC[0], wfc(ik, 0, 0), H_wfc_XC(ik, 0, 0));

        std::cout << "\n\n\nHkPsi pass!\n\n\n";
        
        // something wrong
        psiDotPsi( ParaV, para_wfc, para_Eij, wfc(ik, 0, 0), H_wfc_TV(ik, 0, 0), Eij, &(wfcHwfc_TV(ik, 0)));
        psiDotPsi( ParaV, para_wfc, para_Eij, wfc(ik, 0, 0), H_wfc_hartree(ik, 0, 0), Eij, &(wfcHwfc_hartree(ik, 0)));
        psiDotPsi( ParaV, para_wfc, para_Eij, wfc(ik, 0, 0), H_wfc_XC(ik, 0, 0), Eij, &(wfcHwfc_XC(ik, 0)));

        std::cout << "\n\n\npsiDotPsi pass!\n\n\n";
        
        set_zero_HK(HK_TV);
        set_zero_HK(HK_hartree);
        set_zero_HK(HK_XC);
    }

    //test
    std::cout << "\n\n\n******\n";
    printMatrix_pointer(nk_total, GlobalV::NBANDS, wfcHwfc_TV(0, 0), "wfcHwfc_TV");
    printMatrix_pointer(nk_total, GlobalV::NBANDS, wfcHwfc_hartree(0, 0), "wfcHwfc_hartree");
    printMatrix_pointer(nk_total, GlobalV::NBANDS, wfcHwfc_XC(0, 0), "wfcHwfc_XC");
    std::cout << "\n******\n\n\n";

    //this would transfer the value of H_wfc_TV, H_wfc_hartree, H_wfc_XC
    add_psi(ParaV, para_wfc, wg, H_wfc_TV, H_wfc_hartree, H_wfc_XC, wg_HamiltWfc);
    add_wg(wg, wfcHwfc_TV, wfcHwfc_hartree, wfcHwfc_XC, wg_wfcHamiltWfc);
    add_wg(wg, wfcHwfc_TV, wfcHwfc_hartree, wfcHwfc_XC, wg_forEtotal, 1);

    double Etotal_RDMFT = sumWg_getEnergy(wg_forEtotal);

    std::cout << "\n\n\nEtotal_RDMFT pass!\n\n\n"; /////////////////////////

    //for E_TV
    // ModuleBase::matrix wg_forETV(wg.nr, wg.nc, true);
    ModuleBase::matrix wg_forETV(nk_total, GlobalV::NBANDS, true);
    wgMul_wfcHwfc(wg, wfcHwfc_TV, wg_forETV, 0);
    double ETV_RDMFT = sumWg_getEnergy(wg_forETV);

    //for Ehartree
    // ModuleBase::matrix wg_forEhartree(wg.nr, wg.nc, true);
    ModuleBase::matrix wg_forEhartree(nk_total, GlobalV::NBANDS, true);
    wgMul_wfcHwfc(wg, wfcHwfc_hartree, wg_forEhartree, 1);
    double Ehartree_RDMFT = sumWg_getEnergy(wg_forEhartree);

    //for Exc
    // ModuleBase::matrix wg_forExc(wg.nr, wg.nc, true);
    ModuleBase::matrix wg_forExc(nk_total, GlobalV::NBANDS, true);
    wgMul_wfcHwfc(wg, wfcHwfc_XC, wg_forExc, 3);
    double Exc_RDMFT = sumWg_getEnergy(wg_forExc);

    std::cout << "\n\n\n******\nEtotal_RDMFT:   " << Etotal_RDMFT << "\nETV_RDMFT: " << ETV_RDMFT << "\nEhartree_RDMFT: " 
                << Ehartree_RDMFT << "\nExc_RDMFT:      " << Exc_RDMFT << "\n******\n\n\n";

    ModuleBase::timer::tick("hamilt_lcao", "RDMFT_E&Egradient");
    
    return Etotal_RDMFT;
    // return 1.0;

}



}

#endif
