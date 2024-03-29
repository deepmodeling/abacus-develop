#include <fstream>
#include <sstream>
#include <omp.h>

#include "kernels/cuda/gint_force.cuh"
#include "gint_force.h"
#include "module_base/ylm.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "kernels/cuda/cuda_tools.cuh"

namespace GintKernel{

// Function to calculate forces using GPU-accelerated gamma point Gint
/**
 * @brief Calculate forces and stresses for the `gint_gamma_force_gpu` function.
 *
 * This function calculates forces and stresses based on given parameters.
 *
 * @param dm A pointer to the HContainer<double> object.
 * @param vfactor The scaling factor for some calculation.
 * @param vlocal A pointer to an array of doubles.
 * @param force A pointer to an array to store the calculated forces.
 * @param stress A pointer to an array to store the calculated stresses.
 * @param nczp An integer representing a parameter.
 * @param ylmcoef_now A pointer to an array of doubles representing Ylm coefficients.
 * @param gridt A reference to a Grid_Technique object.
 */
/**
 * Function to calculate forces using GPU-accelerated gamma point Gint
 * @brief Calculate forces and stresses for the `gint_gamma_force_gpu` function.
 *
 * This function calculates forces and stresses based on given parameters.
 *
 * @param dm Pointer to the HContainer<double> object.
 * @param vfactor The scaling factor for the  gird calculation.
 * @param vlocal One-dimensional array that holds the local potential of each gird.
 * @param force One-dimensional array that holds the force of each gird.
 * @param stress One-dimensional array that holds the stress of each gird.
 * @param nczp The number of grid layers in the C direction.
 * @param ylmcoef_now Pointer to an array of doubles representing Ylm coefficients.
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
void gint_gamma_force_gpu(hamilt::HContainer<double> *dm, 
                            const double vfactor,
                            const double *vlocal, 
                            double *force, 
                            double *stress,
                            const int nczp, 
                            const double *ylmcoef_now,
                            const Grid_Technique &gridt, 
                            const LCAO_Orbitals &ORB,
                            const UnitCell &ucell)
{
    const int nbz = gridt.nbzp;
    const int lgd = gridt.lgd;
    const int max_size = gridt.max_atom;
    const int nwmax = ucell.nwmax;
    const int bxyz = gridt.bxyz;
    const int atomNumOnGrids=nbz*bxyz*max_size;
    const int cudaThreads=256;
    const int cudaBlocks=std::min(64,(gridt.psir_size+cudaThreads-1)/cudaThreads);
    int iter_num = 0;
    DensityMat           densityMat;
    ForceStressIatGlobal forceStressIatG;
    SGridParameter       para;
    ForceStressIat       forceStressIat;
    
    CalculateInit(densityMat,forceStressIatG,dm,gridt,ucell,
                lgd,cudaBlocks,atomNumOnGrids);
    /*cuda stream allocate */
    for (int i = 0; i < gridt.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(gridt.streams[i]));
    }
    
    /*compute the psi*/
    for (int i = 0; i < gridt.nbx; i++)
    {
        for (int j = 0; j < gridt.nby; j++)
        {
            
            int max_m = 0;
            int max_n = 0;
            int atom_pair_num = 0;
            dim3 GridPsi(nbz, 8);
            dim3 blockPsi(64);
            dim3 girdDotForce(cudaBlocks);
            dim3 blockDotForce(cudaThreads);
            dim3 gridDot(cudaBlocks);
            dim3 blockDot(cudaThreads);

            CalculateGridInit(para,iter_num,nbz,gridt);
            ForceStressIatInit(forceStressIat,para.streamNum,cudaBlocks,atomNumOnGrids,
                                max_size,forceStressIatG.stressGlobal,
                                forceStressIatG.forceGlobal,forceStressIatG.iatGlobal);
            checkCuda(cudaStreamSynchronize(gridt.streams[para.streamNum]));

            /*gpu task compute in CPU */
            gpu_task_generator_force(gridt,ORB,ucell,i,j,
                                    gridt.psi_size_max_z,
                                    max_size,nczp,vfactor,vlocal,
                                    forceStressIat.iatHost,lgd,densityMat.densityMatDev,
                                    max_m,max_n,atom_pair_num,
                                    para);
            /*variables memcpy to gpu host*/
            CalculateGridMemCpy(para,gridt,nbz,atomNumOnGrids);
            ForceStressIatMemCpy(forceStressIat,gridt,atomNumOnGrids,cudaBlocks,para.streamNum);

            /* cuda stream compute and Multiplication of multinomial matrices */
            get_psi_force<<<GridPsi, blockPsi, 0, gridt.streams[para.streamNum]>>>(
                gridt.ylmcoef_g, ORB.dr_uniform, gridt.bxyz, ucell.nwmax,para.psi_input_double_g, 
                para.input_int_g, para.num_psirDevice, gridt.psi_size_max_z,
                gridt.ucell_atom_nwl_g, gridt.atom_iw2_new_g,gridt.atom_iw2_ylm_g,gridt.atom_iw2_l_g, 
                gridt.atom_nw_g, gridt.nr_max, gridt.psi_u_g, para.psir_r_device,
                para.psir_lx_device,para.psir_ly_device,para.psir_lz_device, 
                para.psir_lxx_device,para.psir_lxy_device, para.psir_lxz_device, 
                para.psir_lyy_device, para.psir_lyz_device, para.psir_lzz_device);
            checkCudaLastError();
            gridt.fastest_matrix_mul(max_m, max_n, 
                                        para.A_m_device, para.B_n_device,
                                        para.K_device, para.matrix_ADev, 
                                        para.lda_device,para.matrix_BDev, 
                                        para.ldb_device, para.matrix_CDev,
                                        para.ldc_device, atom_pair_num,
                                        gridt.streams[para.streamNum], nullptr);

            /* force compute in GPU */
            dot_product_force<<<girdDotForce, blockDotForce, 0,gridt.streams[para.streamNum]>>>(
                para.psir_lx_device, para.psir_ly_device, 
                para.psir_lz_device, para.psir_dm_device, 
                forceStressIat.forceDev,forceStressIat.iatDev, nwmax, 
                max_size, 
                gridt.psir_size / nwmax);
            /* force compute in CPU*/
            ForceCalculate(forceStressIat,force,atomNumOnGrids);

            /*stress compute in GPU*/
            dot_product_stress<<<gridDot, blockDot, 0, gridt.streams[para.streamNum]>>>(
                        para.psir_lxx_device, para.psir_lxy_device, 
                        para.psir_lxz_device, para.psir_lyy_device,
                        para.psir_lyz_device, para.psir_lzz_device, 
                        para.psir_dm_device, forceStressIat.stressDev,
                        gridt.psir_size);
            /* stress compute in CPU*/
            StressCalculate(forceStressIat,stress,cudaBlocks);
            iter_num++;
        }
    }
    /*free variables in CPU host*/
    for (int i = 0; i < gridt.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(gridt.streams[i]));
    }
}

} // namespace GintKernel
