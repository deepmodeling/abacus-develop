#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_MODULE_GINT_GINT_FORCE_GPU_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_MODULE_GINT_GINT_FORCE_GPU_H

#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
namespace GintKernel
{
void gint_fvl_gpu(const hamilt::HContainer<double>* dm,
                        const double* vlocal,
                        double* force_in,
                        double* stress_in,
                        double dr,
                        const double* rcut,
                        const int isforce,
                        const int isstress,
                        const Grid_Technique& gridt,
                        const UnitCell& ucell);

void gtask_force(const Grid_Technique& gridt,
                 const UnitCell& ucell,
                 const int grid_index_ij,
                 const int nczp,
                 const double vfactor,
                 const double* vlocal_global_value,
                 int& atoms_per_z,
                 int* atoms_num_info,
                 int* iat_on_nbz,
                 uint8_t* atoms_type,
                 double* dr_part,
                 double* vldr3);

void alloc_mult_force(const hamilt::HContainer<double>* dm,
                      const Grid_Technique& gridt,
                      const UnitCell& ucell,
                      const int grid_index_ij,
                      const int max_atom,
                      const int *atoms_num_info,
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
                      double** mat_C);

} // namespace GintKernel
#endif
