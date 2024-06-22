#ifndef GINT_FORCE_GPU_H
#define GINT_FORCE_GPU_H

#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
namespace GintKernel
{
void gint_fvl_gamma_gpu(hamilt::HContainer<double>* dm,
                          const double* vlocal,
                          double* force_in,
                          double* stress_in,
                          double dr,
                          double* rcut,
                          const int isforce,
                          const int isstress,
                          const Grid_Technique& gridt,
                          const UnitCell& ucell);

/**
 * @brief GPU task generator for forces.
 *
 * This function generates GPU tasks for force calculations.
 */

void gtask_force(const Grid_Technique& gridt,
                 const UnitCell& ucell,
                 const int grid_index_ij,
                 const int max_atom_per_bcell,
                 const int max_atom,
                 const int nczp,
                 const double vfactor,
                 const double* rcut,
                 const double* vlocal_global_value,
                 double* psi_input_double,
                 int* psi_input_int,
                 int* atom_num_per_bcell,
                 int* start_idx_per_bcell,
                 int* iat_per_z,
                 int& atom_per_z,
                 std::vector<bool>& gpu_mat_cal_flag);

void alloc_mult_force(const Grid_Technique& gridt,
                      const UnitCell& ucell,
                      const int grid_index_ij,
                      const int max_atom,
                      double* const psi_g,
                      double* const psi_dm_g,
                      double* const dm_matrix_g,
                      int& max_m,
                      int& max_n,
                      int& atom_pair_num,
                      int* mat_m,
                      int* mat_n,
                      int* mat_k,
                      int* mat_lda,
                      int* mat_ldb,
                      int* mat_ldc,
                      double** mat_A,
                      double** mat_B,
                      double** mat_C,
                      const std::vector<bool>& gpu_mat_cal_flag);

} // namespace GintKernel
#endif
