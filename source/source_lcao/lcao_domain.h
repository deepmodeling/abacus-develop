#ifndef LCAO_DOMAIN_H
#define LCAO_DOMAIN_H

#include "source_base/global_function.h"
#include "source_base/vector3.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_lcao/lcao_hs_arrays.h"
#include "source_lcao/force_stress_arrays.h"
#include "source_lcao/module_deepks/lcao_deepks.h"
#include "source_basis/module_ao/parallel_orbitals.h"

namespace LCAO_domain
{

void init_basis_lcao(Parallel_Orbitals& pv,
                     const double& onsite_radius,
                     const double& lcao_ecut,
                     const double& lcao_dk,
                     const double& lcao_dr,
                     const double& lcao_rmax,
                     UnitCell& ucell,
                     TwoCenterBundle& two_center_bundle,
                     LCAO_Orbitals& orb);

void build_Nonlocal_mu_new(const Parallel_Orbitals& pv,
                           ForceStressArrays& fsr, // mohan 2024-06-16
                           double* HlocR,
                           const bool& calc_deri,
                           const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           const TwoCenterIntegrator& intor_orb_beta,
                           const Grid_Driver* GridD);

/**
 * @brief set the elements of force-related matrices in LCAO method
 */
void set_force(const Parallel_Orbitals& pv,
               const int& iw1_all,
               const int& iw2_all,
               const double& vx,
               const double& vy,
               const double& vz,
               const char& dtype,
               double* dsloc_x,
               double* dsloc_y,
               double* dsloc_z,
               double* dhloc_fixed_x,
               double* dhloc_fixed_y,
               double* dhloc_fixed_z);

/**
 * @brief read-only environment for building S/T matrix elements.
 *
 * Everything here is fixed for the duration of one build_ST_new call:
 * the basis, the parallel layout, the unit cell and the spin/polarization
 * configuration. Passed by const reference into single_overlap /
 * single_deriv so those functions no longer read global INPUT state.
 */
struct ST_env
{
    const LCAO_Orbitals& orb;
    const TwoCenterBundle& two_center_bundle;
    const Parallel_Orbitals& pv;
    const UnitCell& ucell;
    const int nspin;
    const int npol;
    const bool cal_stress;
    const bool gamma_only_local;
};

/**
 * @brief one S/T matrix element <phi_1 | O | phi_2>.
 *
 * All inputs that vary per matrix element inside the build_ST_new loops:
 * the operator type, the global orbital indices, the angular quantum
 * numbers of both centres and their displacement.
 */
struct ST_elem
{
    const char dtype;
    const int iw1_all;
    const int iw2_all;
    const int m1;
    const int m2;
    const int t1;
    const int l1;
    const int n1;
    const int t2;
    const int l2;
    const int n2;
    const ModuleBase::Vector3<double> dtau;
    const int jj;
    const int jj0;
    const int kk;
    const int kk0;
};

/**
 * @brief read-only environment for building the nonlocal <psi|beta><beta|psi>
 * matrix elements.
 *
 * Everything here is fixed for the duration of one build_Nonlocal_mu_new
 * call: the parallel layout, the unit cell and the spin/polarization
 * configuration, plus whether derivatives (forces) are requested. Passed by
 * const reference into the per-element accumulation helpers so they no longer
 * read global INPUT state.
 */
struct NL_env
{
    const Parallel_Orbitals& pv;
    const UnitCell& ucell;
    const int nspin;
    const int npol;
    const bool gamma_only_local;
    const bool calc_deri;
};

/**
 * @brief identity of one nonlocal matrix element <psi1|beta><beta|psi2>.
 *
 * The per-matrix-element inputs inside the build_Nonlocal_mu_new loops: the
 * global orbital indices of both centres, the projector atom type and the
 * flat sparse slot this element accumulates into.
 */
struct NL_elem
{
    const int iw1_all;
    const int iw2_all;
    const int t0;
    const int nnr;
};

/**
 * @brief set each element without derivatives
 */
void single_overlap(const ST_env& env,
                    const ST_elem& e,
                    int& nnr,       // output value
                    int& total_nnr, // output value
                    double* olm,    // output value
                    double* HSloc); // output value

/**
 * @brief set each element of T matrices
 */
void single_deriv(const ST_env& env,
                  const ST_elem& e,
                  ForceStressArrays& fsr,
                  int& nnr,       // output value
                  int& total_nnr, // output value
                  double* olm);   // output value

/**
 * @brief set the elements of S and T matrices
 */
void build_ST_new(ForceStressArrays& fsr,
                  const char& dtype,
                  const bool& cal_deri,
                  const bool& cal_stress,
                  const UnitCell& ucell,
                  const LCAO_Orbitals& orb,
                  const Parallel_Orbitals& pv,
                  const TwoCenterBundle& two_center_bundle,
                  const Grid_Driver* GridD,
                  double* SHlocR,
                  bool cal_syns = false,
                  double dmax = 0.0);

/**
 * @brief set zeros for HSR matrices
 */
void zeros_HSR(const char& mtype, LCAO_HS_Arrays& HS_arrays);

#ifdef __MLALGO
template <typename T>
void DeePKS_init(const UnitCell& ucell,
                 Parallel_Orbitals& pv,
                 const int& nks,
                 const LCAO_Orbitals& orb,
                 LCAO_Deepks<T>& ld,
                 std::ofstream& ofs);
#endif

template <typename T>
void set_mat2d(const int& global_ir, const int& global_ic, const T& v, const Parallel_Orbitals& pv, T* mat);

} // namespace LCAO_domain

#endif
