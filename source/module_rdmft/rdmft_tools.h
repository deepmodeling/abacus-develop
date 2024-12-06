//==========================================================
// Author: Jingang Han
// DATE : 2024-03-11
//==========================================================
#ifndef RDMFT_TOOLS_H
#define RDMFT_TOOLS_H

#include "module_psi/psi.h"
#include "module_base/matrix.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_base/blas_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_base/parallel_2d.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_elecstate/module_dm/density_matrix.h"

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hs_matrix_k.hpp"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"


#ifdef __EXX
#include "module_ri/RI_2D_Comm.h"
#include "module_ri/Exx_LRI.h"
// there are some operator reload to print data in different formats
#include "module_ri/test_code/test_function.h"
#endif

#include <iostream>
#include <type_traits>
#include <complex>
#include <vector>
#include <iomanip>



namespace rdmft
{


// //for print matrix
// template <typename TK>
// void printMatrix_pointer(int M, int N, const TK* matrixA, std::string nameA)
// {
//     std::cout << "\n" << nameA << ": \n";
//     for(int i=0; i<M; ++i)
//     {
//         for(int j=0; j<N; ++j)
//         {
//             if( j%5 == 0 ) { std::cout << "\n";
// }
//             std::cout << *(matrixA+i*N+j) << " ";
//         }
//         std::cout << "\n";
//     }
//     std::cout << std::endl;
// }


// template <typename TK>
// void printMatrix_vector(int M, int N, const std::vector<TK>& matrixA, std::string nameA)
// {
//     std::cout << "\n" << nameA << ": \n";
//     for(int i=0; i<M; ++i)
//     {
//         for(int j=0; j<N; ++j)
//         {
//             if( j%5 == 0 ) { std::cout << "\n";
// }
//             std::cout << matrixA[i*N+j] << " ";
//         }
//         std::cout << "\n\n";
//     }
//     std::cout << std::endl;
// }


//! now support XC_func_rdmft = "hf", "muller", "power", "pbe", "pbe0". "wp22" and "cwp22" is realizing.
// for the dft-xc-functional part of xc-functional, just use the default is right! Or don't use the function
double occNum_func(double eta, int symbol = 0, const std::string XC_func_rdmft = "hf", const double alpha_power = 1.0);


template <typename TK>
void conj_psi(psi::Psi<TK>& wfc)
{
    TK* pwfc = &wfc(0, 0, 0);
    for(int i=0; i<wfc.size(); ++i) { pwfc[i] = std::conj( pwfc[i] );
}
}


template <>
void conj_psi<double>(psi::Psi<double>& wfc);


// wfc and H_wfc need to be k_firest and provide wfc(ik, 0, 0) and H_wfc(ik, 0, 0)
//! implement matrix multiplication of Hk^dagger and psi
template <typename TK>
void HkPsi(const Parallel_Orbitals* ParaV, const TK& HK, const TK& wfc, TK& H_wfc)
{

    const int one_int = 1;
    //const double one_double = 1.0, zero_double = 0.0;
    const std::complex<double> one_complex = {1.0, 0.0};
    const std::complex<double> zero_complex = {0.0, 0.0};
    const char N_char = 'N';
    const char C_char = 'C';    // Using 'C' is consistent with the formula

#ifdef __MPI
    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    //because wfc(bands, basis'), H(basis, basis'), we do wfc*H^T(in the perspective of cpp, not in fortran). And get H_wfc(bands, basis) is correct.
    pzgemm_( &C_char, &N_char, &nbasis, &nbands, &nbasis, &one_complex, &HK, &one_int, &one_int, ParaV->desc,
        &wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_complex, &H_wfc, &one_int, &one_int, ParaV->desc_wfc );
#endif
}


template <>
void HkPsi<double>(const Parallel_Orbitals* ParaV, const double& HK, const double& wfc, double& H_wfc);


//! implement matrix multiplication of sum_mu conj(wfc(ik, m ,mu)) * op_wfc(ik, n, mu)
template <typename TK>
void cal_bra_op_ket(const Parallel_Orbitals* ParaV, const Parallel_2D& para_Eij_in, const TK& wfc, const TK& H_wfc, std::vector<TK>& Dmn)
{
    const int one_int = 1;
    const std::complex<double> one_complex = {1.0, 0.0};
    const std::complex<double> zero_complex = {0.0, 0.0};
    const char N_char = 'N';
    const char C_char = 'C';

    const int nrow_bands = para_Eij_in.get_row_size();
    const int ncol_bands = para_Eij_in.get_col_size();

#ifdef __MPI
    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    pzgemm_( &C_char, &N_char, &nbands, &nbands, &nbasis, &one_complex, &wfc, &one_int, &one_int, ParaV->desc_wfc,
            &H_wfc, &one_int, &one_int, ParaV->desc_wfc, &zero_complex, &Dmn[0], &one_int, &one_int, para_Eij_in.desc );
#endif
}


template <>
void cal_bra_op_ket<double>(const Parallel_Orbitals* ParaV, const Parallel_2D& para_Eij_in,
                                const double& wfc, const double& H_wfc, std::vector<double>& Dmn);


//! for Dmn that conforms to the 2d-block rule, get its diagonal elements
template <typename TK>
void _diagonal_in_serial(const Parallel_2D& para_Eij_in, const std::vector<TK>& Dmn, double* wfcHwfc)
{
    const int nrow_bands = para_Eij_in.get_row_size();
    const int ncol_bands = para_Eij_in.get_col_size();

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


//! realize occNum_wfc = occNum * wfc. Calling this function and we can get wfc = occNum*wfc.
template <typename TK>
void occNum_MulPsi(const Parallel_Orbitals* ParaV, const ModuleBase::matrix& occ_number, psi::Psi<TK>& wfc, int symbol = 0,
                const std::string XC_func_rdmft = "hf", const double alpha = 1.0)
{
    const int nk_local = wfc.get_nk();
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    // const int nbasis = ParaV->desc[2];      // need to be deleted
    // const int nbands = ParaV->desc_wfc[3];

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


//! add psi with eta and g(eta)
template <typename TK>
void add_psi(const Parallel_Orbitals* ParaV, 
                const K_Vectors* kv, 
                const ModuleBase::matrix& occ_number, 
                psi::Psi<TK>& psi_TV, 
                psi::Psi<TK>& psi_hartree,
                psi::Psi<TK>& psi_dft_XC, 
                psi::Psi<TK>& psi_exx_XC, 
                psi::Psi<TK>& occNum_Hpsi, 
                const std::string XC_func_rdmft = "hf", 
                const double alpha = 1.0)
{
    const int nk = psi_TV.get_nk();
    const int nbn_local = psi_TV.get_nbands();
    const int nbs_local = psi_TV.get_nbasis();
    occNum_MulPsi(ParaV, occ_number, psi_TV);
    occNum_MulPsi(ParaV, occ_number, psi_hartree);
    occNum_MulPsi(ParaV, occ_number, psi_dft_XC);
    occNum_MulPsi(ParaV, occ_number, psi_exx_XC, 2, XC_func_rdmft, alpha);

    // const int nbasis = ParaV->desc[2];
    // const int nbands = ParaV->desc_wfc[3];

    for(int ik=0; ik<nk; ++ik)
    {
        for(int inbn=0; inbn<nbn_local; ++inbn)
        {
            TK* p_occNum_Hpsi = &( occNum_Hpsi(ik, inbn, 0) );
            for(int inbs=0; inbs<nbs_local; ++inbs)
            {
                p_occNum_Hpsi[inbs] = psi_TV(ik, inbn, inbs) + psi_hartree(ik, inbn, inbs) + psi_dft_XC(ik, inbn, inbs) + psi_exx_XC(ik, inbn, inbs);
            }

            // test, consider the wk into psi or dE/d(wfc)
            BlasConnector::scal(nbs_local, kv->wk[ik], p_occNum_Hpsi, 1);
        }
    }

}

/**
 * @brief occNum_wfcHwfc = occNum*wfcHwfc + occNum_wfcHwfc
 * @param symbol: When symbol = 0, 1, 2, 3, 4, occNum = occNum, 0.5*occNum, g(occNum), 0.5*g(occNum), d_g(occNum)/d_occNum respectively.
 *                Default symbol=0.
*/
void occNum_Mul_wfcHwfc(const ModuleBase::matrix& occ_number, 
                            const ModuleBase::matrix& wfcHwfc, 
                            ModuleBase::matrix& occNum_wfcHwfc,
                            int symbol = 0, 
                            const std::string XC_func_rdmft = "hf", 
                            const double alpha = 1.0);


/**
 * @brief Default symbol = 0 for the gradient of Etotal with respect to occupancy,
 *        symbol = 1 for the relevant calculation of Etotal
*/
void add_occNum(const K_Vectors& kv, 
                    const ModuleBase::matrix& occ_number, 
                    const ModuleBase::matrix& wfcHwfc_TV_in, 
                    const ModuleBase::matrix& wfcHwfc_hartree_in,
                    const ModuleBase::matrix& wfcHwfc_dft_XC_in, 
                    const ModuleBase::matrix& wfcHwfc_exx_XC_in, 
                    ModuleBase::matrix& occNum_wfcHwfc, 
                    const std::string XC_func_rdmft = "hf", 
                    const double alpha = 1.0);


//! do wk*g(occNum)*wfcHwfc and add for TV, hartree, XC. This function just use once, so it can be replace and delete
void add_wfcHwfc(const ModuleBase::matrix& wg, 
                    const ModuleBase::matrix& wk_fun_occNum, 
                    const ModuleBase::matrix& wfcHwfc_TV_in, 
                    const ModuleBase::matrix& wfcHwfc_hartree_in,
                    const ModuleBase::matrix& wfcHwfc_XC_in, 
                    ModuleBase::matrix& occNum_wfcHwfc, 
                    const std::string XC_func_rdmft, 
                    const double alpha);


//! give certain occNum_wfcHwfc, get the corresponding energy
double getEnergy(const ModuleBase::matrix& occNum_wfcHwfc);






//! this part of the code is copying from class Veff and do some modifications.
template <typename TK, typename TR>
class Veff_rdmft : public hamilt::OperatorLCAO<TK, TR>
{
  public:
    /**
     * @brief Construct a new Veff object for multi-kpoint calculation
     * @param GK_in: the pointer of Gint_k object, used for grid integration
    */
    Veff_rdmft(Gint_k* GK_in,
                      hamilt::HS_Matrix_K<TK>* hsk_in,
                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                      elecstate::Potential* pot_in,
                      hamilt::HContainer<TR>* hR_in,
                      const UnitCell* ucell_in,
                      const std::vector<double>& orb_cutoff,
                      Grid_Driver* GridD_in,
                      const int& nspin,
                      const Charge* charge_in,
                      const ModulePW::PW_Basis* rho_basis_in,
                      const ModuleBase::matrix* vloc_in,
                      const ModuleBase::ComplexMatrix* sf_in,
                      const std::string potential_in,
                      double* etxc_in = nullptr,
                      double* vtxc_in = nullptr
                    )
        : GK(GK_in),
          orb_cutoff_(orb_cutoff),
          pot(pot_in),
          ucell(ucell_in),
          gd(GridD_in),
          hamilt::OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in),
          charge_(charge_in),
          rho_basis_(rho_basis_in),
          vloc_(vloc_in),
          sf_(sf_in),
          potential_(potential_in),
          etxc(etxc_in),
          vtxc(vtxc_in)
    {
        this->cal_type = hamilt::calculation_type::lcao_gint;

        this->initialize_HR(ucell_in, GridD_in);

        GK_in->initialize_pvpR(*ucell_in, GridD_in, nspin);
    }
    Veff_rdmft(Gint_Gamma* GG_in,
                      hamilt::HS_Matrix_K<TK>* hsk_in,
                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                      elecstate::Potential* pot_in,
                      hamilt::HContainer<TR>* hR_in,
                      const UnitCell* ucell_in,
                      const std::vector<double>& orb_cutoff,
                      Grid_Driver* GridD_in,
                      const int& nspin,
                      const Charge* charge_in,
                      const ModulePW::PW_Basis* rho_basis_in,
                      const ModuleBase::matrix* vloc_in,
                      const ModuleBase::ComplexMatrix* sf_in,  
                      const std::string potential_in,
                      double* etxc_in = nullptr,
                      double* vtxc_in = nullptr
                    )
        : GG(GG_in), 
          orb_cutoff_(orb_cutoff),
          pot(pot_in),
          hamilt::OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in),
          ucell(ucell_in),
          gd(GridD_in),
          charge_(charge_in),
          rho_basis_(rho_basis_in),
          vloc_(vloc_in),
          sf_(sf_in),
          potential_(potential_in),
          etxc(etxc_in),
          vtxc(vtxc_in)
    {
        this->cal_type = hamilt::calculation_type::lcao_gint;

        this->initialize_HR(ucell_in, GridD_in);

        GG_in->initialize_pvpR(*ucell_in, GridD_in, nspin);
    }

    ~Veff_rdmft<TK, TR>(){};

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|V_{eff}|phi_{\nu, R}>
     * the contribution of V_{eff} is calculated by the contribution of V_{H} and V_{XC} and V_{local pseudopotential} and so on.
     * grid integration is used to calculate the contribution Hamiltonian of effective potential
     */
    virtual void contributeHR() override;

    const UnitCell* ucell;

    Grid_Driver* gd;

  private:
    // used for k-dependent grid integration.
    Gint_k* GK = nullptr;

    // used for gamma only algorithms.
    Gint_Gamma* GG = nullptr;

    std::vector<double> orb_cutoff_;

    // Charge calculating method in LCAO base and contained grid base calculation: DM_R, DM, pvpR_reduced

    elecstate::Potential* pot = nullptr;

    int nspin = 1;
    int current_spin = 0;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the electronic kinetic matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(const UnitCell* ucell_in, Grid_Driver* GridD_in);


    // added by jghan

    const Charge* charge_;

    std::string potential_;

    const ModulePW::PW_Basis* rho_basis_;

    const ModuleBase::matrix* vloc_;

    const ModuleBase::ComplexMatrix* sf_;

    double* etxc;

    double* vtxc;

};




}

#endif
