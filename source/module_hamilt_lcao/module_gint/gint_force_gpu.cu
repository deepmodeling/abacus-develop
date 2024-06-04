#include <omp.h>

#include <fstream>
#include <sstream>

#include "gint_force.h"
#include "kernels/cuda/cuda_tools.cuh"
#include "kernels/cuda/gint_force.cuh"
#include "module_base/ylm.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
namespace GintKernel
{

// Function to calculate forces using GPU-accelerated gamma point Gint
/**
 * @brief Calculate forces and stresses for the `gint_fvl_gamma_gpu` function.
 *
 * This function calculates forces and stresses based on given parameters.
 *
 * @param dm A pointer to the HContainer<double> object.
 * @param vfactor The scaling factor for some calculation.
 * @param vlocal A pointer to an array of doubles.
 * @param force A pointer to an array to store the calculated forces.
 * @param stress A pointer to an array to store the calculated stresses.
 * @param nczp An integer representing a parameter.
 * @param ylmcoef_now A pointer to an array of doubles representing Ylm
 * coefficients.
 * @param gridt A reference to a Grid_Technique object.
 */
/**
 * Function to calculate forces using GPU-accelerated gamma point Gint
 * @brief Calculate forces and stresses for the `gint_fvl_gamma_gpu` function.
 *
 * This function calculates forces and stresses based on given parameters.
 *
 * @param dm Pointer to the HContainer<double> object.
 * @param vfactor The scaling factor for the  gird calculation.
 * @param vlocal One-dimensional array that holds the local potential of each
 * gird.
 * @param force One-dimensional array that holds the force of each gird.
 * @param stress One-dimensional array that holds the stress of each gird.
 * @param nczp The number of grid layers in the C direction.
 * @param dr distance cut in calculate
 * @param rcut distance for each atom orbits
 * @param gridt The Grid_Technique object containing grid information.
 *
 * @note The grid integration on the GPU is mainly divided into the following
 * steps:
 * 1. Use the CPU to divide the grid integration into subtasks.
 * 2. Copy the subtask information to the GPU.
 * 3. Calculate the matrix elements on the GPU.
 * 4. Perform matrix multiplication on the GPU.
 * 5. stress dot on the GPU.
 * 6. force dot on the GPU.
 * 7. Copy the results back to the host.
 */
void gint_fvl_gamma_gpu(hamilt::HContainer<double>* dm,
                          const double vfactor,
                          const double* vlocal,
                          std::vector<double>& force,
                          std::vector<double>& stress,
                          const int nczp,
                          double dr,
                          double* rcut,
                          const int isforce,
                          const int isstress,
                          const Grid_Technique& gridt,
                          const UnitCell& ucell)
{
    const int nbz = gridt.nbzp;
    const int lgd = gridt.lgd;
    const int max_size = gridt.max_atom;
    const int nwmax = ucell.nwmax;
    const int bxyz = gridt.bxyz;
    const int atom_num_grid = nbz * bxyz * max_size;
    const int cuda_threads = 256;
    const int cuda_block
        = std::min(64, (gridt.psir_size + cuda_threads - 1) / cuda_threads);
    int iter_num = 0;
    DensityMat denstiy_mat;
    frc_strs_iat_gbl f_s_iat_dev;
    grid_para para;
    frc_strs_iat f_s_iat;

    calculateInit(denstiy_mat,
                  f_s_iat_dev,
                  dm,
                  gridt,
                  ucell,
                  lgd,
                  cuda_block,
                  atom_num_grid);
    /*cuda stream allocate */
    for (int i = 0; i < gridt.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(gridt.streams[i]));
    }

    /*compute the psi*/
    #pragma omp parallel for num_threads(gridt.nstreams) private(para,f_s_iat) collapse(2)
    for (int i = 0; i < gridt.nbx; i++)
    {
        for (int j = 0; j < gridt.nby; j++)
        {

            int max_m = 0;
            int max_n = 0;
            int atom_pair_num = 0;
            dim3 grid_psi(nbz, 32);
            dim3 block_psi(64);
            dim3 grid_dot_force(cuda_block);
            dim3 block_dot_force(cuda_threads);
            dim3 grid_dot(cuda_block);
            dim3 block_dot(cuda_threads);
            
            int pipeline_index = omp_get_thread_num();
            para_init(para, iter_num, nbz, pipeline_index,gridt);
            cal_init(f_s_iat,
                               pipeline_index,
                               cuda_block,
                               atom_num_grid,
                               max_size,
                               f_s_iat_dev);
            /*gpu task compute in CPU */
            gpu_task_generator_force(gridt,
                                     ucell,
                                     i,
                                     j,
                                     gridt.psi_size_max_z,
                                     max_size,
                                     nczp,
                                     vfactor,
                                     rcut,
                                     vlocal,
                                     f_s_iat.iat_host,
                                     lgd,
                                     denstiy_mat.density_mat_d,
                                     max_m,
                                     max_n,
                                     atom_pair_num,
                                     para);
            /*variables memcpy to gpu host*/
            para_mem_copy(para, 
                                 gridt, 
                                 nbz, 
                                 pipeline_index,
                                 atom_num_grid);
            cal_mem_cpy(f_s_iat,
                                 gridt,
                                 atom_num_grid,
                                 cuda_block,
                                 pipeline_index);
            checkCuda(cudaStreamSynchronize(gridt.streams[pipeline_index]));
            /* cuda stream compute and Multiplication of multinomial matrices */
            
            get_psi_force<<<grid_psi,
                            block_psi,
                            0,
                            gridt.streams[pipeline_index]>>>(
                gridt.ylmcoef_g,
                dr,
                gridt.bxyz,
                ucell.nwmax,
                para.input_double_g,
                para.input_int_g,
                para.num_psir_g,
                gridt.psi_size_max_z,
                gridt.atom_nwl_g,
                gridt.atom_new_g,
                gridt.atom_ylm_g,
                gridt.atom_l_g,
                gridt.atom_nw_g,
                gridt.nr_max,
                gridt.psi_u_g,
                para.psir_r_device,
                para.psir_lx_device,
                para.psir_ly_device,
                para.psir_lz_device,
                para.psir_lxx_device,
                para.psir_lxy_device,
                para.psir_lxz_device,
                para.psir_lyy_device,
                para.psir_lyz_device,
                para.psir_lzz_device);
            checkCudaLastError();
            gridt.fastest_matrix_mul(max_m,
                                     max_n,
                                     para.A_m_device,
                                     para.B_n_device,
                                     para.K_device,
                                     para.matrix_A_device,
                                     para.lda_device,
                                     para.matrix_B_device,
                                     para.ldb_device,
                                     para.matrix_C_device,
                                     para.ldc_device,
                                     atom_pair_num,
                                     gridt.streams[pipeline_index],
                                     nullptr);
            /* force compute in GPU */
            if (isforce){
            dot_product_force<<<grid_dot_force,
                                block_dot_force,
                                0,
                                gridt.streams[pipeline_index]>>>(
                para.psir_lx_device,
                para.psir_ly_device,
                para.psir_lz_device,
                para.psir_dm_device,
                f_s_iat.force_device,
                f_s_iat.iat_device,
                nwmax,
                max_size,
                gridt.psir_size / nwmax);
            }
            /*stress compute in GPU*/
            if (isstress){
            dot_product_stress<<<grid_dot,
                                 block_dot,
                                 0,
                                 gridt.streams[pipeline_index]>>>(
                para.psir_lxx_device,
                para.psir_lxy_device,
                para.psir_lxz_device,
                para.psir_lyy_device,
                para.psir_lyz_device,
                para.psir_lzz_device,
                para.psir_dm_device,
                f_s_iat.stress_device,
                gridt.psir_size);
            }
            /* stress compute in CPU*/
            if (isstress){
                cal_stress_add(f_s_iat, stress, cuda_block);
            }
            if (isforce){
                cal_force_add(f_s_iat, force, atom_num_grid);
            }
            iter_num++;
            delete[] f_s_iat.stress_host;
            delete[] f_s_iat.force_host;
            delete[] f_s_iat.iat_host;
        }
    }
    delete[] denstiy_mat.density_mat_h;
    /*free variables in CPU host*/
    for (int i = 0; i < gridt.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(gridt.streams[i]));
    }
}

} // namespace GintKernel
