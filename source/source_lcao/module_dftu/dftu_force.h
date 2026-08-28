/// @file dftu_force.h
/// @brief Free-function helpers for DFT+U force and stress, extracted from
///        Plus_U. The top-level force_stress takes a Plus_U& because it needs
///        to call Plus_U::pot_onsite_real/complex; the four inner
///        functions are fully decoupled and take their dependencies as
///        explicit parameters (mirroring the folding helpers in the same
///        DFTU_LCAO namespace).
#ifndef DFTU_FORCE_H
#define DFTU_FORCE_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_base/matrix.h"
#include "source_base/vector3.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/force_stress_arrays.h"

#include <complex>
#include <string>
#include <vector>

#ifdef __LCAO

class Plus_U;

namespace DFTU_LCAO {

/// @brief Top-level entry: drives force/stress from DFT+U.
/// Takes Plus_U& because it calls dftu.pot_onsite_real/complex,
/// which are still members of Plus_U (defined in dftu_tools.cpp).
void force_stress(Plus_U& dftu,
                  const UnitCell& ucell,
                  const Grid_Driver& gd,
                  std::vector<std::vector<double>>* dmk_d,
                  std::vector<std::vector<std::complex<double>>>* dmk_c,
                  const Parallel_Orbitals& pv,
                  ForceStressArrays& fsr,
                  ModuleBase::matrix& force_dftu,
                  ModuleBase::matrix& stress_dftu,
                  const K_Vectors& kv,
                  const int npol);

/// @brief Force contribution at a k-point (multik path).
void cal_force_k(int nlocal,
                 int npol,
                 const std::string& ks_solver,
                 const std::vector<double>& orb_cutoff,
                 const std::vector<int>& orbital_corr,
                 const std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>& iatlnmipol2iwt,
                 const UnitCell& ucell,
                 const Grid_Driver& gd,
                 ForceStressArrays& fsr,
                 const Parallel_Orbitals& pv,
                 const int ik,
                 const std::complex<double>* rho_VU,
                 ModuleBase::matrix& force_dftu,
                 const ModuleBase::Vector3<double>& kvec_d);

/// @brief Stress contribution at a k-point (multik path).
void cal_stress_k(int nlocal,
                  int npol,
                  const std::string& ks_solver,
                  const std::vector<double>& orb_cutoff,
                  const UnitCell& ucell,
                  const Grid_Driver& gd,
                  ForceStressArrays& fsr,
                  const Parallel_Orbitals& pv,
                  const int ik,
                  const std::complex<double>* rho_VU,
                  ModuleBase::matrix& stress_dftu,
                  const ModuleBase::Vector3<double>& kvec_d);

/// @brief Force contribution at gamma point.
void cal_force_gamma(int nlocal,
                     int npol,
                     const std::vector<int>& orbital_corr,
                     const std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>& iatlnmipol2iwt,
                     const UnitCell& ucell,
                     const double* rho_VU,
                     const Parallel_Orbitals& pv,
                     double* dsloc_x,
                     double* dsloc_y,
                     double* dsloc_z,
                     ModuleBase::matrix& force_dftu);

/// @brief Stress contribution at gamma point.
void cal_stress_gamma(int nlocal,
                      int npol,
                      const std::string& ks_solver,
                      const std::vector<double>& orb_cutoff,
                      const UnitCell& ucell,
                      const Parallel_Orbitals& pv,
                      const Grid_Driver* gd,
                      double* dsloc_x,
                      double* dsloc_y,
                      double* dsloc_z,
                      double* dh_r,
                      const double* rho_VU,
                      ModuleBase::matrix& stress_dftu);

} // namespace DFTU_LCAO

#endif // __LCAO

#endif // DFTU_FORCE_H
