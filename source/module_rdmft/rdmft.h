//==========================================================
// Author: Jingang Han
// DATE : 2024-03-11
//==========================================================
#ifndef RDMFT_H
#define RDMFT_H

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
    
    Parallel_Orbitals* ParaV = nullptr;
    Parallel_2D para_Eij;
    
    // used for k-dependent grid integration.
    Gint_k* GK = nullptr;
    // used for gamma only algorithms.
    Gint_Gamma* GG = nullptr;

    K_Vectors kv;

    ModuleBase::matrix occ_number;
    psi::Psi<TK> wfc;
    ModuleBase::matrix wg;  // //wg is global matrix, wg.nr = nk_total, wg.nc = GlobalV::NBANDS
    ModuleBase::matrix wk_fun_occNum;

    // store the gradients of Etotal with respect to the natural occupation numbers and wfc respectively
    ModuleBase::matrix occNum_wfcHamiltWfc;
    psi::Psi<TK> occNum_HamiltWfc;

    hamilt::HContainer<TR> HR_TV;   // (GlobalC::ucell, ParaV)
    hamilt::HContainer<TR> HR_hartree;
    hamilt::HContainer<TR> HR_XC;

    std::vector<TK> HK_TV;  // ( ParaV->get_row_size()*ParaV->get_col_size() )
    std::vector<TK> HK_hartree;
    std::vector<TK> HK_XC;

    ModuleBase::matrix Etotal_n_k;  // (wg.nr, wg.nc, true)
    ModuleBase::matrix wfcHwfc_TV;
    ModuleBase::matrix wfcHwfc_hartree;
    ModuleBase::matrix wfcHwfc_XC;

    psi::Psi<TK> H_wfc_TV;  // (nk_total, nbands_local, nbasis_local)
    psi::Psi<TK> H_wfc_hartree;
    psi::Psi<TK> H_wfc_XC;

    // just for temperate. in the future when realize psiDotPsi() without pzgemm_/pdgemm_,we don't need it
    // const int nrow_bands = para_Eij.get_row_size();
    // const int ncol_bands = para_Eij.get_col_size();
    std::vector<TK> Eij_TV;     // (nrow_bands*ncol_bands)
    std::vector<TK> Eij_hartree;
    std::vector<TK> Eij_XC;

    hamilt::OperatorLCAO<TK, TR>* V_ekinetic_potential;
    hamilt::OperatorLCAO<TK, TR>* V_nonlocal;
    hamilt::OperatorLCAO<TK, TR>* V_local;
    hamilt::OperatorLCAO<TK, TR>* V_hartree;
    hamilt::OperatorLCAO<TK, TR>* V_XC;

    Exx_LRI<double> Vxc_fromRI_d;   // (GlobalC::exx_info.info_ri)
    Exx_LRI<std::complex<double>> Vxc_fromRI_c;

    void init(Gint_Gamma* GG_in, Gint_k* GK_in, Parallel_Orbitals* ParaV_in, UnitCell& cell);

  private:
    








};












}




#endif