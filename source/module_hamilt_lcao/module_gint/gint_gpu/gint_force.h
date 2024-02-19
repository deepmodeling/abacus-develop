#ifndef GINT_FORCE_H
#define GINT_FORCE_H
#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"

void gint_gamma_force_gpu(hamilt::HContainer<double> *DM,
                        const double vfactor,
                        const double *vlocal,
                        double *force,
                        double *stress,
                        const int nczp,
                        const double *ylmcoef_now,
                        const Grid_Technique &GridT);

void gpu_task_generator_force(const Grid_Technique &GridT, 
                            const int i, const int j,
                            const int psi_size_max, const int max_size,
                            const int nczp,
                            const double vfactor,
                            const double *vlocal_global_value,
                            int *iat_per_nbz,
                            double *psi_input_double, int *psi_input_int,
                            int *num_psir,
                            const int lgd,
                            double *psir_ylm_g,
                            double *psir_zeros_g,
                            double *dm_matrix_g,
                            int *mat_m,
                            int *mat_n,
                            int *mat_k,
                            int *mat_lda,
                            int *mat_ldb,
                            int *mat_ldc,
                            double **mat_A,
                            double **mat_B,
                            double **mat_C,
                            int &max_m,
                            int &max_n,
                            int &atom_pair_num,
                            double *rho_g,
                            double **vec_l,
                            double **vec_r,
                            double **dot_product,
                            int *vec_len,
                            int &dot_count 
                            );

#endif
