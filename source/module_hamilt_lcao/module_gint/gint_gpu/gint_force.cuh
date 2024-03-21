#ifndef GINT_FORCE_CUH
#define GINT_FORCE_CUH

#include <cuda_runtime.h>
namespace lcaoCudaKernel
{
    
/**
 * @brief GPU kernel to calculate the force.
 *
 * This kernel calculates the force based on provided input parameters.
 *
 * @param ylmcoef         Coefficients for Ylm.
 * @param delta_r_g       Delta r value.
 * @param bxyz_g          Bxyz values.
 * @param nwmax_g         Maximum nw value.
 * @param input_double    Array of input double values.
 * @param input_int       Array of input int values.
 * @param num_psir        Array representing the number of psir.
 * @param psi_size_max    Maximum size of psi.
 * @param ucell_atom_nwl  Array representing the unit cell atom nwl.
 * @param atom_iw2_new    Array representing the atom iw2 new.
 * @param atom_iw2_ylm    Array representing the atom iw2 ylm.
 * @param atom_iw2_l      Array representing the atom iw2 l.
 * @param atom_nw         Array representing the atom nw.
 * @param nr_max          Maximum nr value.
 * @param psi_u           Array representing psi_u values.
 * @param psir_ylm_right  Array representing psir ylm right values.
 * @param dpsir_ylm_left_x Array representing dpsir ylm left x values.
 * @param dpsir_ylm_left_y Array representing dpsir ylm left y values.
 * @param dpsir_ylm_left_z Array representing dpsir ylm left z values.
 * @param ddpsir_ylm_left_xx Array representing ddpsir ylm left xx values.
 * @param ddpsir_ylm_left_xy Array representing ddpsir ylm left xy values.
 * @param ddpsir_ylm_left_xz Array representing ddpsir ylm left xz values.
 * @param ddpsir_ylm_left_yy Array representing ddpsir ylm left yy values.
 * @param ddpsir_ylm_left_yz Array representing ddpsir ylm left yz values.
 * @param ddpsir_ylm_left_zz Array representing ddpsir ylm left zz values.
 */
__global__ void get_psi_force(double *ylmcoef, double delta_r_g, double bxyz_g,
                              double nwmax_g, double *input_double, int *input_int,
                              int *num_psir, int psi_size_max, int *ucell_atom_nwl,
                              bool *atom_iw2_new, int *atom_iw2_ylm, int *atom_iw2_l,
                              int *atom_nw, int nr_max, double *psi_u,
                              double *psir_ylm_right, double *dpsir_ylm_left_x,
                              double *dpsir_ylm_left_y, double *dpsir_ylm_left_z,
                              double *ddpsir_ylm_left_xx, double *ddpsir_ylm_left_xy,
                              double *ddpsir_ylm_left_xz, double *ddpsir_ylm_left_yy,
                              double *ddpsir_ylm_left_yz, double *ddpsir_ylm_left_zz);

/**
 * @brief GPU kernel to calculate the dot product for stress.
 *
 * This kernel calculates the dot product for stress based on provided input parameters.
 *
 * @param n               Array of integers.
 * @param x_array_g       Array of x values.
 * @param incx            Increment for x array.
 * @param y_array_g       Array of y values.
 * @param incy            Increment for y array.
 * @param results_g       Array of results.
 * @param batchcount      Batch count.
 */
__global__ void psir_dot_stress(int *n, double **x_array_g, int incx,
                                double **y_array_g, int incy,
                                double **results_g, int batchcount);

/**
 * @brief GPU kernel to calculate the dot product for stress.
 *
 * This kernel calculates the dot product for stress based on provided input parameters.
 *
 * @param ddpsir_ylm_left_xx Array representing ddpsir ylm left xx values.
 * @param ddpsir_ylm_left_xy Array representing ddpsir ylm left xy values.
 * @param ddpsir_ylm_left_xz Array representing ddpsir ylm left xz values.
 * @param ddpsir_ylm_left_yy Array representing ddpsir ylm left yy values.
 * @param ddpsir_ylm_left_yz Array representing ddpsir ylm left yz values.
 * @param ddpsir_ylm_left_zz Array representing ddpsir ylm left zz values.
 * @param psir_ylm_dm      Array representing psir ylm dm values.
 * @param stress_dot       Array representing stress dot values.
 * @param elements_num     Number of elements.
 */
__global__ void dot_product_stress(double *ddpsir_ylm_left_xx,
                                   double *ddpsir_ylm_left_xy,
                                   double *ddpsir_ylm_left_xz,
                                   double *ddpsir_ylm_left_yy,
                                   double *ddpsir_ylm_left_yz,
                                   double *ddpsir_ylm_left_zz,
                                   double *psir_ylm_dm, double *stress_dot,
                                   int elements_num);

/**
 * @brief GPU kernel to calculate the dot product for force.
 *
 * This kernel calculates the dot product for force based on provided input parameters.
 *
 * @param dpsir_ylm_left_x Array representing dpsir ylm left x values.
 * @param dpsir_ylm_left_y Array representing dpsir ylm left y values.
 * @param dpsir_ylm_left_z Array representing dpsir ylm left z values.
 * @param psir_ylm_dm      Array representing psir ylm dm values.
 * @param force_dot        Array representing force dot values.
 * @param iat              Array representing iat values.
 * @param nwmax            Maximum nw value.
 * @param max_size         Maximum size value.
 * @param elements_num     Number of elements.
 */
__global__ void dot_product_force(double *dpsir_ylm_left_x,
                                  double *dpsir_ylm_left_y,
                                  double *dpsir_ylm_left_z, double *psir_ylm_dm,
                                  double *force_dot, int *iat, int nwmax,
                                  int max_size, int elements_num);

} // namespace lcaoCudaKernel
#endif // GINT_VL_CUH