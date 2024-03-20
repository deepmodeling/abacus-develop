//==========================================================
// Author: Jingang Han
// DATE : 2024-03-11
//==========================================================
#ifndef RDMFT_H
#define RDMFT_H

// #include "module_rdmft/rdmft_tools.h"

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
// #include "module_elecstate/module_dm/cal_dm_psi.h"
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
#include <vector>
#include <iomanip>



namespace rdmft
{


template <typename TK, typename TR>
class RDMFT
{
  public:
    RDMFT();
    ~RDMFT();
    
    
    /****** these parameters are passed in from outside, don't need delete ******/

    Parallel_Orbitals* ParaV = nullptr;
    Parallel_2D para_Eij;
    
    // GK and GG are used for multi-k grid integration and gamma only algorithms respectively
    Gint_k* GK = nullptr;
    Gint_Gamma* GG = nullptr;
    Charge* charge = nullptr;

    // update after ion step
    UnitCell* ucell = nullptr;
    K_Vectors* kv = nullptr;
    LCAO_Matrix* LM = nullptr;
    ModulePW::PW_Basis* rho_basis = nullptr;
    ModuleBase::matrix* vloc = nullptr;
    ModuleBase::ComplexMatrix* sf = nullptr;
    Local_Orbital_Charge* loc = nullptr;  // would be delete in the future

    /****** these parameters are passed in from outside, don't need delete ******/


    int nk_total = 0;
    std::string XC_func_rdmft;
    double alpha_power = 0.656; // 0.656 for soilds, 0.525 for dissociation of H2, 0.55~0.58 for HEG

    // natrual occupation numbers and wavefunction
    ModuleBase::matrix occ_number;
    psi::Psi<TK> wfc;
    ModuleBase::matrix wg;
    ModuleBase::matrix wk_fun_occNum;

    // store the gradients of Etotal with respect to the natural occupation numbers and wfc respectively
    ModuleBase::matrix occNum_wfcHamiltWfc;
    psi::Psi<TK> occNum_HamiltWfc;

    // E_RDMFT[4] stores ETV, Ehartree, Exc, Etotal respectively
    double E_RDMFT[4] = {0.0};
    // std::vector<double> E_RDMFT(4);

    hamilt::HContainer<TR>* HR_TV = nullptr;
    hamilt::HContainer<TR>* HR_hartree = nullptr;
    hamilt::HContainer<TR>* HR_XC = nullptr;

    std::vector<TK> HK_TV;
    std::vector<TK> HK_hartree;
    std::vector<TK> HK_XC;

    ModuleBase::matrix Etotal_n_k;
    ModuleBase::matrix wfcHwfc_TV;
    ModuleBase::matrix wfcHwfc_hartree;
    ModuleBase::matrix wfcHwfc_XC;

    psi::Psi<TK> H_wfc_TV;
    psi::Psi<TK> H_wfc_hartree;
    psi::Psi<TK> H_wfc_XC;

    // just for temperate. in the future when realize psiDotPsi() without pzgemm_/pdgemm_,we don't need it
    std::vector<TK> Eij_TV;
    std::vector<TK> Eij_hartree;
    std::vector<TK> Eij_XC;

    hamilt::OperatorLCAO<TK, TR>* V_ekinetic_potential = nullptr;
    hamilt::OperatorLCAO<TK, TR>* V_nonlocal = nullptr;
    hamilt::OperatorLCAO<TK, TR>* V_local = nullptr;
    hamilt::OperatorLCAO<TK, TR>* V_hartree = nullptr;
    hamilt::OperatorLCAO<TK, TR>* V_XC = nullptr;

    Exx_LRI<double>* Vxc_fromRI_d = nullptr;
    Exx_LRI<std::complex<double>>* Vxc_fromRI_c = nullptr;

    void init(Gint_Gamma& GG_in, Gint_k& GK_in, Parallel_Orbitals& ParaV_in, UnitCell& ucell_in,
                        K_Vectors& kv_in, Charge& charge_in, std::string XC_func_rdmft_in = "HF", double alpha_power_in = 0.656);

    // update in ion-step and get V_TV
    void update_ion(UnitCell& ucell_in, LCAO_Matrix& LM_in, ModulePW::PW_Basis& rho_basis_in,
                        ModuleBase::matrix& vloc_in, ModuleBase::ComplexMatrix& sf_in, Local_Orbital_Charge& loc_in);

    // update in elec-step
    // Or we can use rdmft_solver.wfc/occ_number directly when optimizing, so that the update_elec() function does not require parameters.
    void update_elec(const ModuleBase::matrix& occ_number_in, const psi::Psi<TK>& wfc_in);

    // do all calculation after update occNum&wfc, get Etotal and the gradient of energy with respect to the occNum&wfc
    double Run(ModuleBase::matrix& E_gradient_occNum, psi::Psi<TK>&E_gradient_wfc);



  protected:

    void get_V_TV();

    void get_V_hartree();

    // get the special density matrix DM_XC(nk*nbasis_local*nbasis_local)
    void get_DM_XC(std::vector< std::vector<TK> >& DM_XC);

    // construct V_XC based on different XC_functional( i.e. RDMFT class member XC_func_rdmft)
    void get_V_XC();

    double cal_rdmft();

    void cal_Energy();



  private:
    








};






}




#endif