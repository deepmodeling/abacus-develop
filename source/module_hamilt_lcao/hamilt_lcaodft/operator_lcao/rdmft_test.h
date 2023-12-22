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
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/pot_local.h"
#include "module_elecstate/potentials/pot_xc.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"


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



// this part of the code is copying from class Veff and do some modifications.
template <typename TK, typename TR>
class Veff_rdmft : public OperatorLCAO<TK, TR>
{
  public:
    Veff_rdmft(Gint_k* GK_in,
                      Local_Orbital_Charge* loc_in,
                      LCAO_Matrix* LM_in,
                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                      const Charge* charge_in,
                      hamilt::HContainer<TR>* hR_in,
                      std::vector<TK>* hK_in,
                      const UnitCell* ucell_in,
                      Grid_Driver* GridD_in,
                      const Parallel_Orbitals* paraV,
                      const ModulePW::PW_Basis& rho_basis_in,
                      const ModuleBase::matrix& vloc_in,
                      const ModuleBase::ComplexMatrix& sf_in,
                      std::string& potential_in)
        : GK(GK_in),
          loc(loc_in),
          charge_(charge_in),
          ucell_(ucell_in),
          rho_basis_(rho_basis_in),
          vloc_(vloc_in),
          sf_(sf_in),
          potential_(potential_in),
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
                          const Charge* charge_in,
                          hamilt::HContainer<TR>* hR_in,
                          std::vector<TK>* hK_in,
                          const UnitCell* ucell_in,
                          Grid_Driver* GridD_in,
                          const Parallel_Orbitals* paraV,
                          const ModulePW::PW_Basis& rho_basis_in,
                          const ModuleBase::matrix& vloc_in,
                          const ModuleBase::ComplexMatrix& sf_in,  
                          std::string& potential_in
                          )
        : GG(GG_in), 
          loc(loc_in), 
          charge_(charge_in),
          ucell_(ucell_in),
          rho_basis_(rho_basis_in),
          vloc_(vloc_in),
          sf_(sf_in),
          potential_(potential_in),
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

    // add by jghan
    const UnitCell* ucell_;

    const Charge* charge_;

    std::string potential_;

    const ModulePW::PW_Basis rho_basis_;

    const ModuleBase::matrix vloc_;

    const ModuleBase::ComplexMatrix sf_;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the electronic kinetic matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(const UnitCell* ucell_in, Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);

};


// template class Veff_rdmft<double, double>;

// template class Veff_rdmft<std::complex<double>, double>;

// template class Veff_rdmft<std::complex<double>, std::complex<double>>;






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
    pzgemm_( &T_char, &N_char, &nbasis, &nbands, &nbasis, &one_complex, &HK, &one_int, &one_int, ParaV->desc,
        &wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_complex, &H_wfc, &one_int, &one_int, ParaV->desc_wfc );

    // zgemm_( &T_char, &N_char, &nbasis, &nbands, &nbasis, &one_complex, &HK, &nbasis, &wfc, &nbasis, &zero_complex, &H_wfc, &nbasis );

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
                const TK& wfc, const TK& H_wfc, std::vector<TK>& Dmn, double* wfcHwfc, std::string& test_rank_file, std::ofstream& outFile)
{
    const int one_int = 1;
    //const double one_double = 1.0, zero_double = 0.0;
    const std::complex<double> one_complex = {1.0, 0.0};
    const std::complex<double> zero_complex = {0.0, 0.0};
    const char N_char = 'N';
    const char C_char = 'C'; 

    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    const int nrow_bands = para_Eij_in.get_row_size();
    const int ncol_bands = para_Eij_in.get_col_size();
    // const int nrow_bands = ParaV->nrow_bands;
    // const int ncol_bands = ParaV->ncol_bands;

    //test
    // std::cout << "\n\n\n******\nwfcHwfc[0], wfcHwfc[1]: " << wfcHwfc[0] << " " << wfcHwfc[1] <<"\n";
    // set_zero_HK(Dmn);
    
    // in parallel_orbitals.h, there has int desc_Eij[9] which used for Eij in TDDFT, nbands*nbands. Just proper here.
    // std::vector<TK> Dmn(nrow_bands*ncol_bands);
    // zgemm_( &C_char, &N_char, &nbands, &nbands, &nbasis, &one_complex,  &wfc, &nbasis, &H_wfc, &nbasis, &zero_complex, &Dmn[0], &nbands );
    pzgemm_( &C_char, &N_char, &nbands, &nbands, &nbasis, &one_complex, &wfc, &one_int, &one_int, ParaV->desc_wfc,
            &H_wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_complex, &Dmn[0], &one_int, &one_int, para_Eij_in.desc );

    // int dsizeNow;
    // int rankNow;
    // MPI_Comm_size(MPI_COMM_WORLD, &dsizeNow);
    // MPI_Comm_rank(MPI_COMM_WORLD, &rankNow);
    // if( dsizeNow == 2 && rankNow == 1 )   // test for 2 processors
    // {
    //     for(int i=0; i<nrow_bands; ++i)
    //     {
    //         int i_global = para_Eij_in.local2global_row(i);
    //         for(int j=0; j<ncol_bands; ++j)
    //         {
    //             int j_global = para_Eij_in.local2global_col(j);
    //             if(i_global==j_global && i_global==1)
    //             {
    //                 // wfcHwfc[j_global] = std::real( Dmn[i*ncol_bands+j] ); // need to be sure imag(Dmn[i*nrow_bands+j]) == 0
    //                 wfcHwfc[j_global] = std::real( Dmn[0+1] );
    //                 // double testEnn = std::abs( std::imag( Dmn[i*ncol_bands+j] ) );
    //                 // if( testEnn>1e-16 )std::cout << "\n\nimag(Enn)!=0? imag(Enn)= " << testEnn << "\n\n";
    //                 outFile << "\n(i, j) of Dij: " << i << " " << j << "\n(i, i) of global Eij, i: " << i_global << "\n";
    //             }
    //             else if(i_global==j_global && i_global==3)
    //             {
    //                 // wfcHwfc[j_global] = std::real( Dmn[i*ncol_bands+j] ); // need to be sure imag(Dmn[i*nrow_bands+j]) == 0
    //                 wfcHwfc[j_global] = std::real( Dmn[4*ncol_bands+0] );
    //                 // double testEnn = std::abs( std::imag( Dmn[i*ncol_bands+j] ) );
    //                 // if( testEnn>1e-16 )std::cout << "\n\nimag(Enn)!=0? imag(Enn)= " << testEnn << "\n\n";
    //                 outFile << "\n(i, j) of Dij: " << i << " " << j << "\n(i, i) of global Eij, i: " << i_global << "\n";
    //             }
    //         }
    //     }
    // }
    // else  // generally
    // {
    for(int i=0; i<nrow_bands; ++i)
    {
        int i_global = para_Eij_in.local2global_row(i);
        for(int j=0; j<ncol_bands; ++j)
        {
            int j_global = para_Eij_in.local2global_col(j);
            if(i_global==j_global)
            {   
                // wfcHwfc = &wg(ik, 0)
                wfcHwfc[j_global] = std::real( Dmn[j*nrow_bands+i] ); // need to be sure imag(Dmn[i*nrow_bands+j]) == 0
                // double testEnn = std::abs( std::imag( Dmn[i*ncol_bands+j] ) );
                // if( testEnn>1e-16 )std::cout << "\n\nimag(Enn)!=0? imag(Enn)= " << testEnn << "\n\n";
                outFile << "\n(i, j) of Dij: " << i << " " << j << "\n(i, i) of global Eij, i: " << i_global << "\n";
            }
        }
    }
    // }

    // test
    //std::cout << "after psiDotPsi()\nwfcHwfc[0], wfcHwfc[1]: " << wfcHwfc[0] << " " << wfcHwfc[1] <<"\n******\n\n\n";

    // outFile <<

    outFile << "\n" << "Dmn" << ": \n";
    for(int i=0; i<nrow_bands; ++i)
    {
        for(int j=0; j<ncol_bands; ++j)
        {
            outFile << Dmn[i*ncol_bands+j] << " ";
        }
        outFile << "\n";
    }
    outFile << "\n";
    // outFile.close();
}

template <>
void psiDotPsi<double>(const Parallel_Orbitals* ParaV, const Parallel_2D& para_wfc_in, const Parallel_2D& para_Eij_in,
                        const double& wfc, const double& H_wfc, std::vector<double>& Dmn, double* wfcHwfc, std::string& test_rank_file, std::ofstream& outFile);


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
// !!!just used for k-dependent grid integration. For gamma only algorithms, transfer Gint_k& GK_in to Gint_Gamma& GG_in and use it in Veff<OperatorLCAO<TK, TR>>
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
                        const Charge& charge_in,
                        const ModulePW::PW_Basis& rho_basis_in,
                        const ModuleBase::matrix& vloc_in,
                        const ModuleBase::ComplexMatrix& sf_in)  // delete pot_in parameter later
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

    std::string local_pot = "local";
    OperatorLCAO<TK, TR>* V_local = new Veff_rdmft<TK,TR>(
        &GK_in,
        &loc_in,
        LM_in,
        kvec_d_in,
        &charge_in,
        &HR_TV,
        &HK_TV,
        &GlobalC::ucell,
        &GlobalC::GridD,
        ParaV,
        rho_basis_in,
        vloc_in,
        sf_in,
        local_pot
    );

    std::string hartree_pot = "hartree";
    OperatorLCAO<TK, TR>* V_hartree = new Veff_rdmft<TK,TR>(
        &GK_in,
        &loc_in,
        LM_in,
        kvec_d_in,
        &charge_in,
        &HR_hartree,
        &HK_hartree,
        &GlobalC::ucell,
        &GlobalC::GridD,
        ParaV,
        rho_basis_in,
        vloc_in,
        sf_in,
        hartree_pot
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
    //wg is global matrix, wg.nr = nk_total, wg.nc = GlobalV::NBANDS
    ModuleBase::matrix wg_forEtotal(wg.nr, wg.nc, true);
    ModuleBase::matrix wfcHwfc_TV(wg.nr, wg.nc, true);
    ModuleBase::matrix wfcHwfc_hartree(wg.nr, wg.nc, true);
    ModuleBase::matrix wfcHwfc_XC(wg.nr, wg.nc, true);

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

    /********* sum energy of different mpi rank ***********/
    int mydsize;
    int my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mydsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    std::stringstream ss;
    ss << "processor_" << my_rank << ".txt";
    std::string test_rank_file = ss.str();

    std::ofstream outFile(test_rank_file);
    if (!outFile.is_open())
    {
        std::cerr << "Error opening file: " << test_rank_file << std::endl;
    }

    outFile << "\n******\nnumber of processors: " << mydsize << "\n******\n" ;
    outFile << "rank: " << my_rank << std::endl;
    /********* sum energy of different mpi rank ***********/

    // just for temperate. in the future when realize psiDotPsi() without pzgemm_/pdgemm_,we don't need it

    // const int ncol_bands = ParaV->ncol_bands;
    // std::vector<TK> Eij(ParaV->nloc_Eij);   //////////ncol*ncol? why
    const int nrow_bands = para_Eij.get_row_size();
    const int ncol_bands = para_Eij.get_col_size();
    std::vector<TK> Eij_TV(nrow_bands*ncol_bands);
    std::vector<TK> Eij_hartree(nrow_bands*ncol_bands);
    std::vector<TK> Eij_XC(nrow_bands*ncol_bands);

    //calculate wg_wfcHamiltWfc, wg_HamiltWfc and Etotal
    for(int ik=0; ik<nk_total; ++ik)
    {
        // get the HK with ik-th k vector
        V_ekinetic_potential->contributeHk(ik);
        V_hartree->contributeHk(ik);
        V_XC->contributeHk(ik);

        // if(ik==3)
        // {
        //     std::cout << "\n\n\n******\n";
        //     for(int i=0; i<nbasis_local; ++i)
        //     {
        //         std::cout << "HK_TV, HK_hartree, HK_XC: " << HK_TV[i] << " " << HK_hartree[i] << " " << HK_XC[i] << "\n";
        //     }
        //     std::cout << "******\n\n\n";
        // }

        //
        HkPsi( ParaV, para_wfc, HK_TV[0], wfc(ik, 0, 0), H_wfc_TV(ik, 0, 0));
        HkPsi( ParaV, para_wfc, HK_hartree[0], wfc(ik, 0, 0), H_wfc_hartree(ik, 0, 0));
        HkPsi( ParaV, para_wfc, HK_XC[0], wfc(ik, 0, 0), H_wfc_XC(ik, 0, 0));

        std::cout << "\n\n\nHkPsi pass!\n\n\n";
        
        // something wrong
        psiDotPsi( ParaV, para_wfc, para_Eij, wfc(ik, 0, 0), H_wfc_TV(ik, 0, 0), Eij_TV, &(wfcHwfc_TV(ik, 0)), test_rank_file, outFile);
        psiDotPsi( ParaV, para_wfc, para_Eij, wfc(ik, 0, 0), H_wfc_hartree(ik, 0, 0), Eij_hartree, &(wfcHwfc_hartree(ik, 0)), test_rank_file, outFile);
        psiDotPsi( ParaV, para_wfc, para_Eij, wfc(ik, 0, 0), H_wfc_XC(ik, 0, 0), Eij_XC, &(wfcHwfc_XC(ik, 0)), test_rank_file, outFile);

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
    //Parallel_Reduce::reduce_all(Etotal_RDMFT);


    std::cout << "\n\n\nEtotal_RDMFT pass!\n\n\n";

    //for E_TV
    ModuleBase::matrix wg_forETV(wg.nr, wg.nc, true);
    wgMul_wfcHwfc(wg, wfcHwfc_TV, wg_forETV, 0);
    double ETV_RDMFT = sumWg_getEnergy(wg_forETV);

    //for Ehartree
    ModuleBase::matrix wg_forEhartree(wg.nr, wg.nc, true);
    wgMul_wfcHwfc(wg, wfcHwfc_hartree, wg_forEhartree, 1);
    double Ehartree_RDMFT = sumWg_getEnergy(wg_forEhartree);

    //for Exc
    ModuleBase::matrix wg_forExc(wg.nr, wg.nc, true);
    wgMul_wfcHwfc(wg, wfcHwfc_XC, wg_forExc, 3);
    double Exc_RDMFT = sumWg_getEnergy(wg_forExc);

    Parallel_Reduce::reduce_all(Etotal_RDMFT);
    Parallel_Reduce::reduce_all(ETV_RDMFT);
    Parallel_Reduce::reduce_all(Ehartree_RDMFT);
    Parallel_Reduce::reduce_all(Exc_RDMFT);

    std::cout << "\n\n\n******\nin 0 processor\nEtotal_RDMFT:   " << Etotal_RDMFT << "\nETV_RDMFT: " << ETV_RDMFT << "\nEhartree_RDMFT: " 
                << Ehartree_RDMFT << "\nExc_RDMFT:      " << Exc_RDMFT << "\n******\n\n\n";

    ModuleBase::timer::tick("hamilt_lcao", "RDMFT_E&Egradient");


    /********* sum energy of different mpi rank ***********/
    // int mydsize;
    // int my_rank;
    // MPI_Comm_size(MPI_COMM_WORLD, &mydsize);
    // MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // std::stringstream ss;
    // ss << "processor_" << my_rank << ".txt";
    // std::string test_rank_file = ss.str();

    // std::ofstream outFile(test_rank_file);
    // if (!outFile.is_open())
    // {
    //     std::cerr << "Error opening file: " << test_rank_file << std::endl;
    // }

    outFile << "\n\n\n******\nEtotal_RDMFT:   " << Etotal_RDMFT << "\nETV_RDMFT: " << ETV_RDMFT << "\nEhartree_RDMFT: " 
                << Ehartree_RDMFT << "\nExc_RDMFT:      " << Exc_RDMFT << "\n******\n\n\n"; 

    outFile << "\n";

    outFile.close();
    /********* sum energy of different mpi rank ***********/


    
    return Etotal_RDMFT;
    // return 1.0;

}


}

#endif
