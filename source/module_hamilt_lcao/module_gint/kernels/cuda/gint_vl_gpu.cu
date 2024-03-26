#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_hamilt_lcao/module_gint/kernels/cuda/gint_vl.h"
#include "module_hamilt_lcao/module_gint/kernels/cuda/gint_vl.cuh"
#include "module_base/ylm.h"
#include "vbatch_matrix_multiple/vbatch_matrix_mul.cuh"
#include "vbatch_matrix_multiple/cuda_tools.cuh"
#include <omp.h>

namespace GintKernel{

/**
 * Computes the gamma component of the VL (Vlocal) integral on the GPU.
 *
 * @param hRGint Pointer to the HContainer<double> object to store the computed integrals.
 * @param lgd Dimension information for the computation results.
 * @param max_size The maximum number of neighboring atoms for a grid point.
 * @param vfactor Related to volume. The scaling factor for the Vlocal integrals.
 * @param vlocal Pointer to the Vlocal array.
 * @param ylmcoef_now Pointer to the Ylm coefficients array.
 * @param nczp The number of grid layers in the C direction.
 * @param nbxx The total number of grid points.
 * @param dr The grid spacing.
 * @param rcut Pointer to the cutoff radius array.
 * @param gridt The Grid_Technique object containing grid information.
 * @param ucell The UnitCell object containing unit cell information.
 * 
 * @note The grid integration on the GPU is mainly divided into the following steps:
 * 1. Use the CPU to divide the grid integration into subtasks.
 * 2. Copy the subtask information to the GPU.
 * 3. Calculate the matrix elements on the GPU.
 * 4. Perform matrix multiplication on the GPU.
 * 5. Copy the results back to the host.
 */
void gint_gamma_vl_gpu(hamilt::HContainer<double> *hRGint,
                       const int lgd,
                       const int max_size,
                       const double vfactor,
                       const double *vlocal,
                       const double *ylmcoef_now,
                       const int nczp,
                       const int nbxx,
                       const double dr,
                       const double *rcut,
                       const Grid_Technique &gridt,
                       const UnitCell &ucell)
{
    const int nbz = gridt.nbzp;
    checkCuda(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
    {
        int iter_num = 0;
        for (int iat1 = 0; iat1 < ucell.nat; iat1++)
        {
            for (int iat2 = 0; iat2 < ucell.nat; iat2++)
            {
                int stream_num = iter_num % gridt.nstreams;
                int it1 = ucell.iat2it[iat1];
                int lo1 = gridt.trace_lo[ucell.itiaiw2iwt(
                    it1, ucell.iat2ia[iat1], 0)];

                int it2 = ucell.iat2it[iat2];
                int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(
                    it2, ucell.iat2ia[iat2], 0)];

                if (lo1 <= lo2) {
                    hamilt::AtomPair<double> *tmp_ap = hRGint->find_pair(iat1, iat2);
                    if (tmp_ap == nullptr)
                    {
                        continue;
                    }
                    int atom_pair_nw = ucell.atoms[it1].nw * ucell.atoms[it2].nw;
                    if (gridt.GridVlocal_v2_g[iat1 * ucell.nat + iat2] == nullptr)
                    {
                        checkCuda(cudaMallocAsync((void **)&gridt.GridVlocal_v2_g[iat1 * ucell.nat + iat2],
                                            atom_pair_nw * sizeof(double), gridt.streams[stream_num]));
                    }
                    checkCuda(cudaMemsetAsync(gridt.GridVlocal_v2_g[iat1 * ucell.nat + iat2],
                                            0,
                                            atom_pair_nw * sizeof(double), gridt.streams[stream_num]));
                    iter_num++;
                }
            }
        }
    }
    for (int i = 0; i < gridt.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(gridt.streams[i]));
    }

    #pragma omp parallel for num_threads(gridt.nstreams) collapse(2)
    for (int i = 0; i < gridt.nbx; i++)
    {
        for (int j = 0; j < gridt.nby; j++)
        {
            int stream_num = omp_get_thread_num();
            checkCuda(cudaStreamSynchronize(gridt.streams[stream_num]));
            double *psi_input_double = &gridt.psi_input_double_global[gridt.psi_size_max * stream_num * 5];
            int *psi_input_int = &gridt.psi_input_int_global[gridt.psi_size_max * stream_num * 2];
            int *num_psir = &gridt.num_psir_global[nbz * stream_num];
            int *atom_pair_A_m = &gridt.atom_pair_left_info_global[gridt.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_B_n = &gridt.atom_pair_right_info_global[gridt.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_k = &gridt.atom_pair_k_info_global[gridt.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_lda = &gridt.atom_pair_lda_global[gridt.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_ldb = &gridt.atom_pair_ldb_global[gridt.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_ldc = &gridt.atom_pair_ldc_global[gridt.atom_pair_size_over_nbz * stream_num];

            double *psi_input_double_g = &gridt.psi_input_double_global_g[gridt.psi_size_max * stream_num * 5];
            int *psi_input_int_g = &gridt.psi_input_int_global_g[gridt.psi_size_max * stream_num * 2];
            int *num_psir_g = &gridt.num_psir_global_g[nbz * stream_num];
            double *psir_ylm_left_g = &gridt.psir_ylm_left_global_g[gridt.psir_size * stream_num];
            double *psir_ylm_right_g = &gridt.psir_ylm_right_global_g[gridt.psir_size * stream_num];

            int *atom_pair_A_m_g = &gridt.atom_pair_left_info_global_g[gridt.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_B_n_g = &gridt.atom_pair_right_info_global_g[gridt.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_k_g = &gridt.atom_pair_k_info_global_g[gridt.atom_pair_size_over_nbz * stream_num];

            int *atom_pair_lda_g = &gridt.atom_pair_lda_global_g[gridt.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_ldb_g = &gridt.atom_pair_ldb_global_g[gridt.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_ldc_g = &gridt.atom_pair_ldc_global_g[gridt.atom_pair_size_over_nbz * stream_num];

            double **atom_pair_mat_A_array = &gridt.atom_pair_left_global[gridt.atom_pair_size_over_nbz * stream_num];
            double **atom_pair_mat_B_array = &gridt.atom_pair_right_global[gridt.atom_pair_size_over_nbz * stream_num];
            double **atom_pair_mat_C_array = &gridt.atom_pair_output_global[gridt.atom_pair_size_over_nbz * stream_num];

            double **atom_pair_mat_A_array_g = &gridt.atom_pair_left_global_g[gridt.atom_pair_size_over_nbz * stream_num];
            double **atom_pair_mat_B_array_g = &gridt.atom_pair_right_global_g[gridt.atom_pair_size_over_nbz * stream_num];
            double **atom_pair_mat_C_array_g = &gridt.atom_pair_output_global_g[gridt.atom_pair_size_over_nbz * stream_num];
            int atom_pair_num = 0;
            int max_m = 0;
            int max_n = 0;

            gtask_vlocal(gridt, 
                         rcut, ucell,
                         i, j,
                         max_size, nczp,
                         vfactor,
                         vlocal,
                         psir_ylm_left_g,
                         psir_ylm_right_g,
                         psi_input_double,
                         psi_input_int,
                         num_psir,
                         atom_pair_A_m,
                         atom_pair_B_n,
                         atom_pair_lda,
                         atom_pair_ldb,
                         atom_pair_ldc,
                         atom_pair_mat_A_array,
                         atom_pair_mat_B_array,
                         atom_pair_mat_C_array,
                         atom_pair_num,
                         max_m,
                         max_n);

            for(int z = 0; z < gridt.atom_pair_size_over_nbz; z++){
                atom_pair_k[z] = gridt.bxyz;
            }

            checkCuda(cudaMemcpyAsync(psi_input_double_g, psi_input_double, gridt.psi_size_max * 5 * sizeof(double), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(psi_input_int_g, psi_input_int, gridt.psi_size_max * 2 * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(num_psir_g, num_psir, nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));


            checkCuda(cudaMemcpyAsync(atom_pair_A_m_g, atom_pair_A_m, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_B_n_g, atom_pair_B_n, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_k_g, atom_pair_k, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_lda_g, atom_pair_lda, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_ldb_g, atom_pair_ldb, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_ldc_g, atom_pair_ldc, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_mat_A_array_g, atom_pair_mat_A_array, gridt.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_mat_B_array_g, atom_pair_mat_B_array, gridt.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, gridt.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_mat_C_array_g, atom_pair_mat_C_array, gridt.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, gridt.streams[stream_num]));

            checkCuda(cudaMemsetAsync(psir_ylm_left_g, 0, gridt.psir_size * sizeof(double), gridt.streams[stream_num]));
            checkCuda(cudaMemsetAsync(psir_ylm_right_g, 0, gridt.psir_size * sizeof(double), gridt.streams[stream_num]));

            dim3 grid_psi(nbz, 8);
            dim3 block_psi(64);

            get_psi_and_vldr3<<<grid_psi, block_psi, 0, gridt.streams[stream_num]>>>(gridt.ylmcoef_g,
                                                                                dr,
                                                                                gridt.bxyz,
                                                                                ucell.nwmax,
                                                                                psi_input_double_g,
                                                                                psi_input_int_g,
                                                                                num_psir_g,
                                                                                gridt.psi_size_max_per_z,
                                                                                gridt.ucell_atom_nwl_g,
                                                                                gridt.atom_iw2_new_g,
                                                                                gridt.atom_iw2_ylm_g,
                                                                                gridt.atom_nw_g,
                                                                                gridt.nr_max,
                                                                                gridt.psi_u_g,
                                                                                psir_ylm_left_g,
                                                                                psir_ylm_right_g);
            checkCudaLastError();
            gridt.fastest_matrix_mul(max_m, max_n,
                                     atom_pair_A_m_g, atom_pair_B_n_g, atom_pair_k_g,
                                     atom_pair_mat_A_array_g, atom_pair_lda_g,
                                     atom_pair_mat_B_array_g, atom_pair_ldb_g,
                                     atom_pair_mat_C_array_g, atom_pair_ldc_g,
                                     atom_pair_num, gridt.streams[stream_num], nullptr);
            // checkCuda(cudaStreamSynchronize(gridt.streams[stream_num]));
        }
    }
    for (int i = 0; i < gridt.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(gridt.streams[i]));
    }
    {
        int iter_num = 0;
        for (int iat1 = 0; iat1 < ucell.nat; iat1++)
        {
            for (int iat2 = 0; iat2 < ucell.nat; iat2++)
            {
                int stream_num = iter_num % gridt.nstreams;
                int it1 = ucell.iat2it[iat1];
                int lo1 = gridt.trace_lo[ucell.itiaiw2iwt(
                    it1, ucell.iat2ia[iat1], 0)];

                int it2 = ucell.iat2it[iat2];
                int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(
                    it2, ucell.iat2ia[iat2], 0)];
                if (lo1 <= lo2) {
                    int atom_pair_nw = ucell.atoms[it1].nw * ucell.atoms[it2].nw;
                    hamilt::AtomPair<double> *tmp_ap = hRGint->find_pair(iat1, iat2);
                    if (tmp_ap == nullptr)
                    {
                        continue;
                    }
                    checkCuda(cudaMemcpyAsync(tmp_ap->get_pointer(0),
                                            gridt.GridVlocal_v2_g[iat1 * ucell.nat + iat2],
                                            atom_pair_nw * sizeof(double),
                                            cudaMemcpyDeviceToHost, gridt.streams[stream_num]));
                    iter_num++;
                }
            }
        }
    }
    for (int i = 0; i < gridt.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(gridt.streams[i]));
    }
}

} // namespace GintKernel