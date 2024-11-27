//==========================================================
// Author: Jingang Han
// DATE : 2024-03-11
//==========================================================
#ifndef RDMFT_H
#define RDMFT_H

#include "module_parameter/parameter.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_psi/psi.h"
#include "module_base/matrix.h"

#include "module_basis/module_ao/parallel_2d.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_basis/module_nao/two_center_bundle.h"

#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hs_matrix_k.hpp"

#ifdef __EXX
#include "module_ri/Exx_LRI.h"
#include "module_ri/module_exx_symmetry/symmetry_rotation.h"
// there are some operator reload to print data in different formats
#include "module_ri/test_code/test_function.h"
#endif

#include <iostream>
#include <type_traits>
#include <complex>
#include <vector>
#include <iomanip>

//! Reduced Density Matrix Functional Theory (RDMFT)
namespace rdmft
{

template <typename TK, typename TR>
class RDMFT
{

  public:
    RDMFT();
    ~RDMFT();
    
    const Parallel_Orbitals* ParaV = nullptr;

    Parallel_2D para_Eij;

    //! obain Ewald and this->pelec->pot
    elecstate::ElecState* pelec = nullptr;

    //! update after ion step
    const K_Vectors* kv = nullptr; 

    int nk_total = 0;
    int nbands_total;
    int nspin = 1;
    std::string XC_func_rdmft;

    //! 0.656 for soilds, 0.525 for dissociation of H2, 0.55~0.58 for HEG
    double alpha_power = 0.656; 

    //! natrual occupation numbers and wavefunction
    ModuleBase::matrix occ_number;
    psi::Psi<TK> wfc;
    ModuleBase::matrix wg;
    ModuleBase::matrix wk_fun_occNum;

    //! gradients of total energy with respect to the natural occupation numbers and wfc
    ModuleBase::matrix occNum_wfcHamiltWfc;
    psi::Psi<TK> occNum_HamiltWfc;

    //! E_RDMFT[4] stores ETV, Ehartree, Exc, Etotal
    double E_RDMFT[4] = {0.0};
    double Etotal = 0.0;
    // std::vector<double> E_RDMFT(4);

    //! initialization of rdmft calculation
    void init(Gint_Gamma& GG_in, Gint_k& GK_in, Parallel_Orbitals& ParaV_in, UnitCell& ucell_in,
                        K_Vectors& kv_in, elecstate::ElecState& pelec_in, LCAO_Orbitals& orb_in, TwoCenterBundle& two_center_bundle_in, std::string XC_func_rdmft_in, double alpha_power_in);

    //! update in ion-step and get V_TV
    void update_ion(UnitCell& ucell_in, ModulePW::PW_Basis& rho_basis_in,
                        ModuleBase::matrix& vloc_in, ModuleBase::ComplexMatrix& sf_in);

    //! update in elec-step
    // Or we can use rdmft_solver.wfc/occ_number directly when optimizing, so that the update_elec() function does not require parameters.
    void update_elec(const ModuleBase::matrix& occ_number_in, const psi::Psi<TK>& wfc_in, const Charge* charge_in = nullptr);

    //! obtain the gradient of total energy with respect to occupation number and wfc
    double cal_E_grad_wfc_occ_num();

    void cal_Energy(const int cal_type = 1);

    //! update occ_number for optimization algorithms that depend on Hamilton
    void update_occNumber(const ModuleBase::matrix& occ_number_in);

    //! update occ_number for optimization algorithms that depend on Hamilton
    void update_wg(const ModuleBase::matrix& wg_in);

    //! do all calculation after update occNum&wfc, get Etotal and the gradient of energy with respect to the occNum&wfc
    double run(ModuleBase::matrix& E_gradient_occNum, psi::Psi<TK>& E_gradient_wfc);


  protected:

    //! get the special density matrix DM_XC(nk*nbasis_local*nbasis_local)
    void get_DM_XC(std::vector< std::vector<TK> >& DM_XC);

    void cal_V_TV();

    void cal_V_hartree();

    //! construct V_XC based on different XC_functional( i.e. RDMFT class member XC_func_rdmft)
    void cal_V_XC();

    //! get the total Hamilton in k-space
    void cal_Hk_Hpsi();
    
    void update_charge();

  private:

    //! Hamiltonian matrices in real space
    hamilt::HContainer<TR>* HR_TV = nullptr;
    hamilt::HContainer<TR>* HR_hartree = nullptr;
    hamilt::HContainer<TR>* HR_dft_XC = nullptr;
    hamilt::HContainer<TR>* HR_exx_XC = nullptr;
    // hamilt::HContainer<TR>* HR_local = nullptr;

    //! Hamiltonian matrices in reciprocal space
    hamilt::HS_Matrix_K<TK>* hsk_TV = nullptr;
    hamilt::HS_Matrix_K<TK>* hsk_hartree = nullptr;
    hamilt::HS_Matrix_K<TK>* hsk_dft_XC = nullptr;
    hamilt::HS_Matrix_K<TK>* hsk_exx_XC = nullptr;

    std::vector<TK> HK_XC;
    std::vector< std::vector<TK> > DM_XC_pass;
    // ModuleDirectMin::ProdStiefelVariable HK_RDMFT_pass;
    // ModuleDirectMin::ProdStiefelVariable HK_XC_pass;

    ModuleBase::matrix Etotal_n_k;
    ModuleBase::matrix wfcHwfc_TV;
    ModuleBase::matrix wfcHwfc_hartree;
    ModuleBase::matrix wfcHwfc_XC;
    ModuleBase::matrix wfcHwfc_exx_XC;
    ModuleBase::matrix wfcHwfc_dft_XC;

    psi::Psi<TK> H_wfc_TV;
    psi::Psi<TK> H_wfc_hartree;
    psi::Psi<TK> H_wfc_XC;
    psi::Psi<TK> H_wfc_exx_XC;
    psi::Psi<TK> H_wfc_dft_XC;

    std::vector<TK> Eij_TV;
    std::vector<TK> Eij_hartree;
    std::vector<TK> Eij_XC;
    std::vector<TK> Eij_exx_XC;

    hamilt::OperatorLCAO<TK, TR>* V_ekinetic_potential = nullptr;
    hamilt::OperatorLCAO<TK, TR>* V_nonlocal = nullptr;
    hamilt::OperatorLCAO<TK, TR>* V_local = nullptr;
    hamilt::OperatorLCAO<TK, TR>* V_hartree = nullptr;
    hamilt::OperatorLCAO<TK, TR>* V_exx_XC = nullptr;
    hamilt::OperatorLCAO<TK, TR>* V_dft_XC = nullptr;

#ifdef __EXX
    Exx_LRI<double>* Vxc_fromRI_d = nullptr;
    Exx_LRI<std::complex<double>>* Vxc_fromRI_c = nullptr;
    ModuleSymmetry::Symmetry_rotation symrot_exx;
    bool exx_spacegroup_symmetry = false;
#endif

    double etxc = 0.0;
    double vtxc = 0.0;
    bool only_exx_type = false;
    const int cal_E_type = 1;   // cal_type = 2 just support XC-functional without exx

    /****** these parameters are passed in from outside, don't need delete ******/
    // GK and GG are used for multi-k grid integration and gamma only algorithms respectively
    Gint_k* GK = nullptr;
    Gint_Gamma* GG = nullptr;
    Charge* charge = nullptr;

    // update after ion step
    const UnitCell* ucell = nullptr;
    const ModulePW::PW_Basis* rho_basis = nullptr;
    const ModuleBase::matrix* vloc = nullptr;
    const ModuleBase::ComplexMatrix* sf = nullptr;
    const LCAO_Orbitals* orb = nullptr;
    const TwoCenterBundle* two_center_bundle = nullptr;
};

}
#endif
