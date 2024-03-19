//==========================================================
// Author: Jingang Han
// DATE : 2024-03-11
//==========================================================
#ifndef RDMFT_TOOLS_H
#define RDMFT_TOOLS_H

#include "module_base/matrix.h"
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
#include "module_base/parallel_reduce.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_elecstate/module_dm/density_matrix.h"

#include "module_hamilt_general/operator.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/ekinetic_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/nonlocal_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/veff_lcao.h"

// used by Exx&LRI
#include "module_ri/RI_2D_Comm.h"
#include "module_ri/Exx_LRI.h"

// there are some operator reload to print data in different formats
#include "module_ri/test_code/test_function.h"

// used by class Veff_rdmft
// #include "module_base/timer.h"
// #include "module_elecstate/potentials/potential_new.h"
// #include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
// //#include "module_hamilt_lcao/module_gint/gint_gamma.h"
// //#include "module_hamilt_lcao/module_gint/gint_k.h"
// //#include "operator_lcao.h"
// //#include "module_cell/module_neighbor/sltk_grid_driver.h"
// //#include "module_cell/unitcell.h"
// #include "module_elecstate/potentials/H_Hartree_pw.h"
// #include "module_elecstate/potentials/pot_local.h"
// #include "module_elecstate/potentials/pot_xc.h"
// #include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"


#include <iostream>
#include <type_traits>
#include <complex>
#include <vector>
#include <iomanip>



namespace rdmft
{


//for print matrix
template <typename TK>
void printMatrix_pointer(int M, int N, const TK* matrixA, std::string nameA)
{
    std::cout << "\n" << nameA << ": \n";
    for(int i=0; i<M; ++i)
    {
        for(int j=0; j<N; ++j)
        {
            if( j%5 == 0 ) std::cout << "\n";
            std::cout << *(matrixA+i*N+j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}


template <typename TK>
void printMatrix_vector(int M, int N, const std::vector<TK>& matrixA, std::string nameA)
{
    std::cout << "\n" << nameA << ": \n";
    for(int i=0; i<M; ++i)
    {
        for(int j=0; j<N; ++j)
        {
            if( j%5 == 0 ) std::cout << "\n";
            std::cout << matrixA[i*N+j] << " ";
        }
        std::cout << "\n\n";
    }
    std::cout << std::endl;
}


// now support XC_func_rdmft = "HF", "Muller", "power" 
double occNum_func(double eta, int symbol = 0, const std::string XC_func_rdmft = "HF", const double alpha_power = 0.656);


template <typename TK>
void set_zero_vector(std::vector<TK>& HK)
{
    for(int i=0; i<HK.size(); ++i) HK[i] = 0.0;
}


template <typename TK>
void set_zero_psi(psi::Psi<TK>& wfc)
{
    TK* pwfc_in = &wfc(0, 0, 0);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024)
    #endif
    
    for(int i=0; i<wfc.size(); ++i) pwfc_in[i] = 0.0;
}


template <typename TK>
void conj_psi(psi::Psi<TK>& wfc)
{
    TK* pwfc = &wfc(0, 0, 0);
    for(int i=0; i<wfc.size(); ++i) pwfc[i] = std::conj( pwfc[i] );
}


template <>
void conj_psi<double>(psi::Psi<double>& wfc);


// wfc and H_wfc need to be k_firest and provide wfc(ik, 0, 0) and H_wfc(ik, 0, 0)
// psi::Psi<TK> psi's (), std::vector<TK> HK's [] operator overloading return TK
template <typename TK>
void HkPsi(const Parallel_Orbitals* ParaV, const TK& HK, const TK& wfc, TK& H_wfc)
{

    const int one_int = 1;
    //const double one_double = 1.0, zero_double = 0.0;
    const std::complex<double> one_complex = {1.0, 0.0};
    const std::complex<double> zero_complex = {0.0, 0.0};
    const char N_char = 'N';
    const char C_char = 'C';

    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    //because wfc(bands, basis'), H(basis, basis'), we do wfc*H^T(in the perspective of cpp, not in fortran). And get H_wfc(bands, basis) is correct.
    pzgemm_( &C_char, &N_char, &nbasis, &nbands, &nbasis, &one_complex, &HK, &one_int, &one_int, ParaV->desc,
        &wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_complex, &H_wfc, &one_int, &one_int, ParaV->desc_wfc );
}


template <>
void HkPsi<double>(const Parallel_Orbitals* ParaV, const double& HK, const double& wfc, double& H_wfc);


template <typename TK>
void psiDotPsi(const Parallel_Orbitals* ParaV, const Parallel_2D& para_Eij_in, const TK& wfc, const TK& H_wfc, std::vector<TK>& Dmn, double* wfcHwfc)
{
    const int one_int = 1;
    const std::complex<double> one_complex = {1.0, 0.0};
    const std::complex<double> zero_complex = {0.0, 0.0};
    const char N_char = 'N';
    const char C_char = 'C'; 

    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    const int nrow_bands = para_Eij_in.get_row_size();
    const int ncol_bands = para_Eij_in.get_col_size();

    pzgemm_( &C_char, &N_char, &nbands, &nbands, &nbasis, &one_complex, &wfc, &one_int, &one_int, ParaV->desc_wfc,
            &H_wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_complex, &Dmn[0], &one_int, &one_int, para_Eij_in.desc );

    for(int i=0; i<nrow_bands; ++i)
    {
        int i_global = para_Eij_in.local2global_row(i);
        for(int j=0; j<ncol_bands; ++j)
        {
            int j_global = para_Eij_in.local2global_col(j);
            if(i_global==j_global)
            {   
                // because the Dmn obtained from pzgemm_() is stored column-major
                wfcHwfc[j_global] = std::real( Dmn[i+j*nrow_bands] );
            }
        }
    }
}

template <>
void psiDotPsi<double>(const Parallel_Orbitals* ParaV, const Parallel_2D& para_wfc_in,
                        const double& wfc, const double& H_wfc, std::vector<double>& Dmn, double* wfcHwfc);


// realize occNum_wfc = occNum * wfc. Calling this function and we can get wfc = occNum*wfc.
template <typename TK>
void occNum_MulPsi(const Parallel_Orbitals* ParaV, const ModuleBase::matrix& occ_number, psi::Psi<TK>& wfc, int symbol = 0,
                const std::string XC_func_rdmft = "HF", const double alpha = 0.656)
{
    const int nk_local = wfc.get_nk();
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    const int nbasis = ParaV->desc[2];      // need to be deleted
    const int nbands = ParaV->desc_wfc[3];

    for (int ik = 0; ik < nk_local; ++ik)
    {
        for (int ib_local = 0; ib_local < nbands_local; ++ib_local)  // ib_local < nbands_local , some problem, ParaV->ncol_bands
        {
            const double occNum_local = occNum_func( occ_number(ik, ParaV->local2global_col(ib_local)), symbol, XC_func_rdmft, alpha);
            TK* wfc_pointer = &(wfc(ik, ib_local, 0));
            BlasConnector::scal(nbasis_local, occNum_local, wfc_pointer, 1);
        }
    }
}


// add psi with eta and g(eta)
template <typename TK>
void add_psi(const Parallel_Orbitals* ParaV, const ModuleBase::matrix& occ_number, psi::Psi<TK>& psi_TV, psi::Psi<TK>& psi_hartree,
                psi::Psi<TK>& psi_XC, psi::Psi<TK>& occNum_Hpsi, const std::string XC_func_rdmft = "HF", const double alpha = 0.656)
{
    const int nk = psi_TV.get_nk();
    const int nbn_local = psi_TV.get_nbands();
    const int nbs_local = psi_TV.get_nbasis();
    occNum_MulPsi(ParaV, occ_number, psi_TV);
    occNum_MulPsi(ParaV, occ_number, psi_hartree);
    occNum_MulPsi(ParaV, occ_number, psi_XC, 2, XC_func_rdmft, alpha);

    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    for(int ik=0; ik<nk; ++ik)
    {
        for(int inbn=0; inbn<nbn_local; ++inbn)
        {
            TK* p_occNum_Hpsi = &( occNum_Hpsi(ik, inbn, 0) );
            for(int inbs=0; inbs<nbs_local; ++inbs)
            {
                p_occNum_Hpsi[inbs] = psi_TV(ik, inbn, inbs) + psi_hartree(ik, inbn, inbs) + psi_XC(ik, inbn, inbs);
            }
        }
    }

}


// occNum_wfcHwfc = occNum*wfcHwfc + occNum_wfcHwfc
// When symbol = 0, 1, 2, 3, 4, occNum = occNum, 0.5*occNum, g(occNum), 0.5*g(occNum), d_g(occNum)/d_occNum respectively. Default symbol=0.
void occNum_Mul_wfcHwfc(const ModuleBase::matrix& occ_number, const ModuleBase::matrix& wfcHwfc, ModuleBase::matrix& occNum_wfcHwfc,
                        int symbol = 0, const std::string XC_func_rdmft = "HF", const double alpha = 0.656);


// Default symbol = 0 for the gradient of Etotal with respect to occupancy
// symbol = 1 for the relevant calculation of Etotal
void add_occNum(const ModuleBase::matrix& occ_number, const ModuleBase::matrix& wfcHwfc_TV_in, const ModuleBase::matrix& wfcHwfc_hartree_in,
            const ModuleBase::matrix& wfcHwfc_XC_in, ModuleBase::matrix& occNum_wfcHwfc, const std::string XC_func_rdmft = "HF", const double alpha = 0.656, int symbol = 0);


// do wk*g(occNum)*wfcHwfc and add for TV, hartree, XC. This function just use once, so it can be replace and delete
void add_wfcHwfc(const std::vector<double>& wk_in, const ModuleBase::matrix& occ_number, const ModuleBase::matrix& wfcHwfc_TV_in, const ModuleBase::matrix& wfcHwfc_hartree_in,
                const ModuleBase::matrix& wfcHwfc_XC_in, ModuleBase::matrix& occNum_wfcHwfc, const std::string XC_func_rdmft, const double alpha);


//give certain occNum_wfcHwfc, get the corresponding energy
double getEnergy(const ModuleBase::matrix& occNum_wfcHwfc);




















// this part of the code is copying from class Veff and do some modifications.
template <typename TK, typename TR>
class Veff_rdmft : public hamilt::OperatorLCAO<TK, TR>
{
  public:
    Veff_rdmft(Gint_k* GK_in,
                      LCAO_Matrix* LM_in,
                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                      const Charge* charge_in,
                      hamilt::HContainer<TR>* hR_in,
                      std::vector<TK>* hK_in,
                      const UnitCell* ucell_in,
                      Grid_Driver* GridD_in,
                      const Parallel_Orbitals* paraV,
                      const ModulePW::PW_Basis* rho_basis_in,
                      const ModuleBase::matrix* vloc_in,
                      const ModuleBase::ComplexMatrix* sf_in,
                      const std::string potential_in)
        : GK(GK_in),
          charge_(charge_in),
          ucell_(ucell_in),
          rho_basis_(rho_basis_in),
          vloc_(vloc_in),
          sf_(sf_in),
          potential_(potential_in),
          hamilt::OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
    {
        this->cal_type = hamilt::lcao_gint;

        this->initialize_HR(ucell_in, GridD_in, paraV);

        GK_in->initialize_pvpR(*ucell_in, GridD_in);
    }
    Veff_rdmft(Gint_Gamma* GG_in,
                          LCAO_Matrix* LM_in,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                          const Charge* charge_in,
                          hamilt::HContainer<TR>* hR_in,
                          std::vector<TK>* hK_in,
                          const UnitCell* ucell_in,
                          Grid_Driver* GridD_in,
                          const Parallel_Orbitals* paraV,
                          const ModulePW::PW_Basis* rho_basis_in,
                          const ModuleBase::matrix* vloc_in,
                          const ModuleBase::ComplexMatrix* sf_in,  
                          const std::string potential_in
                          )
        : GG(GG_in), 
          charge_(charge_in),
          ucell_(ucell_in),
          rho_basis_(rho_basis_in),
          vloc_(vloc_in),
          sf_(sf_in),
          potential_(potential_in),
          hamilt::OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
    {
        this->cal_type = hamilt::lcao_gint;

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

    elecstate::Potential* pot = nullptr;

    void initialize_HR(const UnitCell* ucell_in, Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);

    // add by jghan
    const UnitCell* ucell_;

    const Charge* charge_;

    std::string potential_;

    const ModulePW::PW_Basis* rho_basis_;

    const ModuleBase::matrix* vloc_;

    const ModuleBase::ComplexMatrix* sf_;

};

































}

#endif
