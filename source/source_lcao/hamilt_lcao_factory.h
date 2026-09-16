#ifndef HAMILT_LCAO_FACTORY_H
#define HAMILT_LCAO_FACTORY_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_hamilt/hs_matrix_k.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/setup_deepks.h"
#include "source_pw/module_pwdft/dftu_base.h"

#include <string>
#include <vector>

namespace hamilt
{

/**
 * @brief products of the LCAO operator-chain factory.
 *
 * The factory builds the operator chain on top of the already-allocated
 * hsk/hR/sR buffers (owned by the caller); it only produces the chain head
 * and, for MLALGO builds, the DeePKS V_delta(R) handle.
 */
template <typename TK, typename TR>
struct LcaoOpsBundle
{
    OperatorLCAO<TK, TR>* ops = nullptr;   ///< built operator-chain head
    HContainer<TR>* v_delta_R = nullptr;   ///< DeePKS V_delta(R), MLALGO only
};

/**
 * @brief build the operator chain for the gamma-only case (TK == double).
 *
 * Appends overlap/kinetic/nonlocal/veff nodes (and optional DeePKS/DFTU)
 * onto the buffers hsk/hR/sR. hR is gamma-fixed here.
 *
 * @param kv k-point list (kvec_d/isk read from it)
 * @param hsk target H(k)/S(k) matrix buffer, already allocated by caller
 * @param hR target H(R) container, already allocated by caller
 * @param sR target S(R) container, already allocated by caller
 * @return LcaoOpsBundle with the chain head and DeePKS V_delta(R) handle
 */
template <typename TK, typename TR>
LcaoOpsBundle<TK, TR> build_gamma_ops(const UnitCell& ucell,
                                      const Grid_Driver& grid_d,
                                      const Parallel_Orbitals* paraV,
                                      elecstate::Potential* pot_in,
                                      const TwoCenterBundle& two_center_bundle,
                                      const LCAO_Orbitals& orb,
                                      elecstate::DensityMatrix<TK, double>* DM_in,
                                      Plus_U_Base* p_dftu,
                                      Setup_DeePKS<TK>& deepks,
                                      const Input_para& inp,
                                      const std::vector<std::string>& pot_register_in,
                                      const K_Vectors* kv,
                                      HS_Matrix_K<TK>* hsk,
                                      HContainer<TR>* hR,
                                      HContainer<TR>* sR);

/**
 * @brief build the operator chain for the multi-k case (TK == complex<double>).
 *
 * Appends veff/overlap/kinetic/nonlocal nodes (and optional DeePKS/TDDFT/
 * DFTU/spin-constrain) onto the buffers hsk/hR/sR.
 *
 * @param kv k-point list (kvec_d/isk read from it)
 * @param hsk target H(k)/S(k) matrix buffer, already allocated by caller
 * @param hR target H(R) container, already allocated by caller
 * @param sR target S(R) container, already allocated by caller
 * @return LcaoOpsBundle with the chain head and DeePKS V_delta(R) handle
 */
template <typename TK, typename TR>
LcaoOpsBundle<TK, TR> build_multik_ops(const UnitCell& ucell,
                                       const Grid_Driver& grid_d,
                                       const Parallel_Orbitals* paraV,
                                       elecstate::Potential* pot_in,
                                       const TwoCenterBundle& two_center_bundle,
                                       const LCAO_Orbitals& orb,
                                       elecstate::DensityMatrix<TK, double>* DM_in,
                                       Plus_U_Base* p_dftu,
                                       Setup_DeePKS<TK>& deepks,
                                       const Input_para& inp,
                                       const std::vector<std::string>& pot_register_in,
                                       const K_Vectors* kv,
                                       HS_Matrix_K<TK>* hsk,
                                       HContainer<TR>* hR,
                                       HContainer<TR>* sR);

} // namespace hamilt

#endif
