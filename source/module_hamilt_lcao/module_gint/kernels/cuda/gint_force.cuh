#ifndef GINT_FORCE_CUH
#define GINT_FORCE_CUH

#include <cuda_runtime.h>
namespace GintKernel
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
 * @param psir_r  Array representing psir ylm right values.
 * @param psir_lx Array representing dpsir ylm left x values.
 * @param psir_ly Array representing dpsir ylm left y values.
 * @param psir_lz Array representing dpsir ylm left z values.
 * @param psir_lxx Array representing ddpsir ylm left xx values.
 * @param psir_lxy Array representing ddpsir ylm left xy values.
 * @param psir_lxz Array representing ddpsir ylm left xz values.
 * @param psir_lyy Array representing ddpsir ylm left yy values.
 * @param psir_lyz Array representing ddpsir ylm left yz values.
 * @param psir_lzz Array representing ddpsir ylm left zz values.
 */
__global__ void get_psi_force(double* ylmcoef,
                              double delta_r_g,
                              int bxyz_g,
                              double nwmax_g,
                              double* input_double,
                              int* input_int,
                              int* num_psir,
                              int psi_size_max,
                              int* ucell_atom_nwl,
                              bool* atom_iw2_new,
                              int* atom_iw2_ylm,
                              int* atom_iw2_l,
                              int* atom_nw,
                              int nr_max,
                              double* psi_u,
                              double* psir_r,
                              double* psir_lx,
                              double* psir_ly,
                              double* psir_lz,
                              double* psir_lxx,
                              double* psir_lxy,
                              double* psir_lxz,
                              double* psir_lyy,
                              double* psir_lyz,
                              double* psir_lzz);



/**
 * @brief GPU kernel to calculate the dot product for stress.
 *
 * This kernel calculates the dot product for stress based on provided input
 * parameters.
 *
 * @param psir_lxx Array representing ddpsir ylm left xx values.
 * @param psir_lxy Array representing ddpsir ylm left xy values.
 * @param psir_lxz Array representing ddpsir ylm left xz values.
 * @param psir_lyy Array representing ddpsir ylm left yy values.
 * @param psir_lyz Array representing ddpsir ylm left yz values.
 * @param psir_lzz Array representing ddpsir ylm left zz values.
 * @param psir_ylm_dm      Array representing psir ylm dm values.
 * @param stress_dot       Array representing stress dot values.
 * @param elements_num     Number of elements.
 */
__global__ void dot_product_stress(double* psir_lxx,
                                   double* psir_lxy,
                                   double* psir_lxz,
                                   double* psir_lyy,
                                   double* psir_lyz,
                                   double* psir_lzz,
                                   double* psir_ylm_dm,
                                   double* stress_dot,
                                   int elements_num);

/**
 * @brief GPU kernel to calculate the dot product for force.
 *
 * This kernel calculates the dot product for force based on provided input
 * parameters.
 *
 * @param psir_lx Array representing dpsir ylm left x values.
 * @param psir_ly Array representing dpsir ylm left y values.
 * @param psir_lz Array representing dpsir ylm left z values.
 * @param psir_ylm_dm      Array representing psir ylm dm values.
 * @param force_dot        Array representing force dot values.
 * @param iat              Array representing iat values.
 * @param nwmax            Maximum nw value.
 * @param max_size         Maximum size value.
 * @param elements_num     Number of elements.
 */
__global__ void dot_product_force(double* psir_lx,
                                  double* psir_ly,
                                  double* psir_lz,
                                  double* psir_ylm_dm,
                                  double* force_dot,
                                  int* iat,
                                  int nwmax,
                                  int max_size,
                                  int elements_num);

} // namespace GintKernel
#endif // GINT_VL_CUH
