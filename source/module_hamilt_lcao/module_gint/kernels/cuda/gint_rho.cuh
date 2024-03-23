#ifndef GINT_RHO_CUH
#define GINT_RHO_CUH

#include <cuda_runtime.h>
namespace lcaoCudaKernel{
    
/**
 * @brief CUDA kernel to calculate psir.
 *
 * This kernel calculates the wave function psir using the provided input parameters.
 *
 * @param ylmcoef pointer to the array of Ylm coefficients.
 * @param delta_r_g value of delta_r_g.
 * @param bxyz_g number of meshcells in a bigcell.
 * @param nwmax_g maximum nw.
 * @param input_double the array of `double` type datas used to calculate psir.
 * @param input_int the array of `int` type datas used to calculate psir.
 * @param num_psir array records the number of atoms on each bigcell.
 * @param psi_size_max maximum number of atoms on bigcell.
 * @param ucell_atom_nwl array record nw of each type of atom.
 * @param atom_iw2_new 
 * @param atom_iw2_ylm 
 * @param atom_nw Pointer to the array of atom_nw values.
 * @param nr_max 
 * @param psi_u 
 * @param psir_ylm
 */
__global__ void get_psi(double *ylmcoef,
                        double delta_r_g,
                        double bxyz_g,
                        double nwmax_g,
                        double *input_double,
                        int *input_int,
                        int *num_psir,
                        int psi_size_max,
                        int *ucell_atom_nwl,
                        bool *atom_iw2_new,
                        int *atom_iw2_ylm,
                        int *atom_nw,
                        int nr_max,
                        double *psi_u,
                        double *psir_ylm);

/**
 * @brief Kernel function to calculate batch vector dot products.
 *
 * @param n             array of vector length.
 * @param vec_l_g       array of pointers to left vec.
 * @param incl          stride between consecutive elements in the `vec_l_g`.
 * @param vec_r_g       array of pointers to right vec.
 * @param incr          stride between consecutive elements in the `vec_r_g`.
 * @param results_g     array to store the dot product results.
 * @param batchcount    Number of dot products to compute.
 */
__global__ void psir_dot(int * n,
                        double **vec_l_g,
                        int incl,
                        double **vec_r_g,
                        int incr,
                        double **results_g,
                        int batchcount);

} // namespace lcaoCudaKernel
#endif // GINT_RHO_CUH