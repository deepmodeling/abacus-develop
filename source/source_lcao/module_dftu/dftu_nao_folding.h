/// @file dftu_folding.h
/// @brief Free-function helpers for folding S/dS matrices, extracted from
///        Plus_U. Each function takes the data it needs (orb_cutoff,
///        ks_solver, npol, gamma_only_local, nspin) as direct parameters;
///        no Plus_U reference is required, so the helpers are fully decoupled
///        from the class and unit-testable.
#ifndef DFTU_FOLDING_H
#define DFTU_FOLDING_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/force_stress_arrays.h"
#include "source_hamilt/hamilt.h"

#include <complex>
#include <string>
#include <vector>

#ifdef __LCAO

namespace DFTU_LCAO {

/// @brief Judge whether atom pair (T1,I1) and (T2,I2,tau2) are adjacent
///        by direct orbital cutoff overlap or three-body bridging via a
///        common nonlocal projector center T0.
/// @param orb_cutoff orbital cutoff radii per atom type
/// @return true if the pair should be processed
bool is_adjacent_pair(const std::vector<double>& orb_cutoff,
                      const UnitCell& ucell,
                      const Grid_Driver& gd,
                      int T1,
                      int T2,
                      const ModuleBase::Vector3<double>& tau1,
                      const ModuleBase::Vector3<double>& tau2);

/// @brief Get the linear index of local matrix element (mu, nu) based on
///        ks_solver (column-major or row-major).
int get_linear_index(const std::string& ks_solver,
                     int mu,
                     int nu,
                     const Parallel_Orbitals& pv);

/// @brief Fold the dSR matrix for gamma-only calculations.
///        npol is the spin-polarization factor; orb_cutoff and ks_solver
///        are forwarded to is_adjacent_pair and get_linear_index.
void fold_dSR_gamma(int npol,
                    const std::string& ks_solver,
                    const std::vector<double>& orb_cutoff,
                    const UnitCell& ucell,
                    const Parallel_Orbitals& pv,
                    const Grid_Driver* gd,
                    double* dsloc_x,
                    double* dsloc_y,
                    double* dsloc_z,
                    double* dh_r,
                    int dim1,
                    int dim2,
                    double* dSR_gamma);

// dim1 = 0 : S, for Hamiltonian
// dim1 = 1-3 : dS, for force
// dim1 = 4-6 : dS * dR, for stress
void folding_matrix_k(int npol,
                      const std::string& ks_solver,
                      const std::vector<double>& orb_cutoff,
                      const UnitCell& ucell,
                      const Grid_Driver& gd,
                      ForceStressArrays& fsr,
                      const Parallel_Orbitals& pv,
                      int ik,
                      int dim1,
                      int dim2,
                      std::complex<double>* mat_k,
                      const ModuleBase::Vector3<double>& kvec_d);

/**
 * @brief new function of folding_S_matrix
 * only for Hamiltonian now, for force and stress will be developed later
 * use HContainer as input and output in mat_k
*/
void folding_matrix_k_new(const std::string& ks_solver,
                          bool gamma_only_local,
                          int nspin,
                          int ik,
                          hamilt::Hamilt<std::complex<double>>* p_ham);

} // namespace DFTU_LCAO

#endif // __LCAO

#endif // DFTU_FOLDING_H
