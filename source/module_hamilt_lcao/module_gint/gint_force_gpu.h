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
 *
 * @param gridt Reference to Grid_Technique .
 * @param ucell Reference to UnitCell .
 * @param grid_index_ij Index of the grid.
 * @param psi_size_max Maximum size of psi.
 * @param max_size Maximum size of atoms on a grid.
 * @param nczp Size parameter,stand for the current z-axis grids.
 * @param vfactor Scaling factor,stand for the Local potential.
 * @param rcut distance for each atom orbits
 * @param vlocal_global_value Global values of local potential.
 * @param iat_per_nbz save the number of the iat on per nbz grids.
 * @param atom_pair_num Number of atom pairs,stand for the max number of mat_n.
 * @param gpu_mat_cal_flag Establish whether to perform calculations between
 * atoms and grid points.
 * @param para Grid parameter in task generator,
 */

void gpu_task_generator_force(const Grid_Technique& gridt,
                              const UnitCell& ucell,
                              const int grid_index_ij,
                              const int psiSizeMax,
                              const int max_atom,
                              const int nczp,
                              const double vfactor,
                              const double* rcut,
                              const double* vlocal_global_value,
                              double* psi_input_double,
                              int* psi_input_int,
                              int* phi_num_per_bcell,
                              int* iat_per_z,
                              int& atom_pair_num,
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
                      std::vector<bool>& gpu_mat_cal_flag);

} // namespace GintKernel
#endif
