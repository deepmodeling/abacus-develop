#ifndef GINT_FORCE_H
#define GINT_FORCE_H

#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"

/**
 * @brief Calculate forces using GPU.
 *
 * This function calculates forces and stress for a given set of parameters.
 *
 * @param DM A pointer to hamilt::HContainer<double>.
 * @param vfactor Scaling factor for forces.
 * @param vlocal Local potential values.
 * @param force Output array for forces.
 * @param stress Output array for stress.
 * @param nczp Size parameter.
 * @param ylmcoef_now Coefficients for spherical harmonics.
 * @param GridT Reference to Grid_Technique object.
 */
void gint_gamma_force_gpu(hamilt::HContainer<double> *DM, const double vfactor,
                          const double *vlocal, double *force, double *stress,
                          const int nczp, const double *ylmcoef_now,
                          const Grid_Technique &GridT);

/**
 * @brief GPU task generator for forces.
 *
 * This function generates GPU tasks for force calculations.
 *
 * @param GridT Reference to Grid_Technique object.
 * @param i Value of i.
 * @param j Value of j.
 * @param psi_size_max Maximum size of psi.
 * @param max_size Maximum size.
 * @param nczp Size parameter.
 * @param vfactor Scaling factor.
 * @param vlocal_global_value Global values of local potential.
 * @param iat_per_nbz Array of iat values per nbz.
 * @param psi_input_double Double array for input psi values.
 * @param psi_input_int Integer array for input psi values.
 * @param num_psir Array for num_psir values.
 * @param lgd Value of lgd.
 * @param psir_ylm_g GPU array for psir_ylm.
 * @param psir_zeros_g GPU array for psir_zeros.
 * @param dm_matrix_g GPU array for dm_matrix.
 * @param mat_m Array for mat_m values.
 * @param mat_n Array for mat_n values.
 * @param mat_k Array for mat_k values.
 * @param mat_lda Array for mat_lda values.
 * @param mat_ldb Array for mat_ldb values.
 * @param mat_ldc Array for mat_ldc values.
 * @param mat_A Double pointer for mat_A.
 * @param mat_B Double pointer for mat_B.
 * @param mat_C Double pointer for mat_C.
 * @param max_m Maximum value of m.
 * @param max_n Maximum value of n.
 * @param atom_pair_num Number of atom pairs.
 * @param rho_g GPU array for rho.
 * @param vec_l Double pointer for vec_l.
 * @param vec_r Double pointer for vec_r.
 * @param dot_product Double pointer for dot_product.
 * @param vec_len Array for vec_len values.
 * @param dot_count Reference to dot_count.
 */
void gpu_task_generator_force(
    const Grid_Technique &GridT, const int i, const int j,
    const int psi_size_max, const int max_size, const int nczp,
    const double vfactor, const double *vlocal_global_value, int *iat_per_nbz,
    double *psi_input_double, int *psi_input_int, int *num_psir, const int lgd,
    double *psir_ylm_g, double *psir_zeros_g, double *dm_matrix_g, int *mat_m,
    int *mat_n, int *mat_k, int *mat_lda, int *mat_ldb, int *mat_ldc,
    double **mat_A, double **mat_B, double **mat_C, int &max_m, int &max_n,
    int &atom_pair_num, double *rho_g, double **vec_l, double **vec_r,
    double **dot_product, int *vec_len, int &dot_count);

#endif
