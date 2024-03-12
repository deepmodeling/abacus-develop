#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_vl.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_vl.cuh"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "vbatch_matrix_multiple/vbatch_matrix_mul.cuh"
#include "vbatch_matrix_multiple/cuda_tools.cuh"
#include <omp.h>
void gint_gamma_vl_gpu(hamilt::HContainer<double> *hRGint,
                       const int lgd,
                       const int max_size,
                       const double vfactor,
                       const double *vlocal,
                       const double *ylmcoef_now,
                       const int nczp,
                       const int NLOCAL_now,
                       const int nbxx,
                       const Grid_Technique &GridT)
{

    const int nbz = GridT.nbzp;
    checkCuda(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
    {
        int iter_num = 0;
        for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
        {
            for (int iat2 = 0; iat2 < GlobalC::ucell.nat; iat2++)
            {
                int stream_num = iter_num % GridT.nstreams;
                int it1 = GlobalC::ucell.iat2it[iat1];
                int lo1 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
                    it1, GlobalC::ucell.iat2ia[iat1], 0)];

                int it2 = GlobalC::ucell.iat2it[iat2];
                int lo2 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
                    it2, GlobalC::ucell.iat2ia[iat2], 0)];

                if (lo1 <= lo2) {
                    hamilt::AtomPair<double> *tmp_ap = hRGint->find_pair(iat1, iat2);
                    if (tmp_ap == nullptr)
                    {
                        continue;
                    }
                    int atom_pair_nw = GlobalC::ucell.atoms[it1].nw * GlobalC::ucell.atoms[it2].nw;
                    if (GridT.GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2] == nullptr)
                    {
                        checkCuda(cudaMallocAsync((void **)&GridT.GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2],
                                            atom_pair_nw * sizeof(double), GridT.streams[stream_num]));
                    }
                    checkCuda(cudaMemsetAsync(GridT.GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2],
                                            0,
                                            atom_pair_nw * sizeof(double), GridT.streams[stream_num]));
                    iter_num++;
                }
            }
        }
    }
    for (int i = 0; i < GridT.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(GridT.streams[i]));
    }

    #pragma omp parallel for num_threads(GridT.nstreams) collapse(2)
    for (int i = 0; i < GridT.nbx; i++)
    {
        for (int j = 0; j < GridT.nby; j++)
        {
            int stream_num = omp_get_thread_num();
            checkCuda(cudaStreamSynchronize(GridT.streams[stream_num]));
            double *psi_input_double = &GridT.psi_input_double_global[GridT.psi_size_max * stream_num * 5];
            int *psi_input_int = &GridT.psi_input_int_global[GridT.psi_size_max * stream_num * 2];
            int *num_psir = &GridT.num_psir_global[nbz * stream_num];
            int *atom_pair_A_m = &GridT.atom_pair_left_info_global[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_B_n = &GridT.atom_pair_right_info_global[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_k = &GridT.atom_pair_k_info_global[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_lda = &GridT.atom_pair_lda_global[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_ldb = &GridT.atom_pair_ldb_global[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_ldc = &GridT.atom_pair_ldc_global[GridT.atom_pair_size_over_nbz * stream_num];

            double *psi_input_double_g = &GridT.psi_input_double_global_g[GridT.psi_size_max * stream_num * 5];
            int *psi_input_int_g = &GridT.psi_input_int_global_g[GridT.psi_size_max * stream_num * 2];
            int *num_psir_g = &GridT.num_psir_global_g[nbz * stream_num];
            double *psir_ylm_left_g = &GridT.psir_ylm_left_global_g[GridT.psir_size * stream_num];
            double *psir_ylm_right_g = &GridT.psir_ylm_right_global_g[GridT.psir_size * stream_num];

            int *atom_pair_A_m_g = &GridT.atom_pair_left_info_global_g[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_B_n_g = &GridT.atom_pair_right_info_global_g[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_k_g = &GridT.atom_pair_k_info_global_g[GridT.atom_pair_size_over_nbz * stream_num];

            int *atom_pair_lda_g = &GridT.atom_pair_lda_global_g[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_ldb_g = &GridT.atom_pair_ldb_global_g[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_ldc_g = &GridT.atom_pair_ldc_global_g[GridT.atom_pair_size_over_nbz * stream_num];

            double **atom_pair_mat_A_array = &GridT.atom_pair_left_global[GridT.atom_pair_size_over_nbz * stream_num];
            double **atom_pair_mat_B_array = &GridT.atom_pair_right_global[GridT.atom_pair_size_over_nbz * stream_num];
            double **atom_pair_mat_C_array = &GridT.atom_pair_output_global[GridT.atom_pair_size_over_nbz * stream_num];

            double **atom_pair_mat_A_array_g = &GridT.atom_pair_left_global_g[GridT.atom_pair_size_over_nbz * stream_num];
            double **atom_pair_mat_B_array_g = &GridT.atom_pair_right_global_g[GridT.atom_pair_size_over_nbz * stream_num];
            double **atom_pair_mat_C_array_g = &GridT.atom_pair_output_global_g[GridT.atom_pair_size_over_nbz * stream_num];
            int atom_pair_num = 0;
            int max_m = 0;
            int max_n = 0;

            gpu_task_generate_vlocal(GridT, i, j,
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

            for(int z = 0; z < GridT.atom_pair_size_over_nbz; z++){
                atom_pair_k[z] = GridT.bxyz;
            }

            checkCuda(cudaMemcpyAsync(psi_input_double_g, psi_input_double, GridT.psi_size_max * 5 * sizeof(double), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(psi_input_int_g, psi_input_int, GridT.psi_size_max * 2 * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(num_psir_g, num_psir, nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));


            checkCuda(cudaMemcpyAsync(atom_pair_A_m_g, atom_pair_A_m, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_B_n_g, atom_pair_B_n, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_k_g, atom_pair_k, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_lda_g, atom_pair_lda, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_ldb_g, atom_pair_ldb, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_ldc_g, atom_pair_ldc, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_mat_A_array_g, atom_pair_mat_A_array, GridT.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_mat_B_array_g, atom_pair_mat_B_array, GridT.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_mat_C_array_g, atom_pair_mat_C_array, GridT.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            checkCuda(cudaMemsetAsync(psir_ylm_left_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(psir_ylm_right_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));

            dim3 grid_psi(nbz, 8);
            dim3 block_psi(64);

            get_psi_and_vldr3<<<grid_psi, block_psi, 0, GridT.streams[stream_num]>>>(GridT.ylmcoef_g,
                                                                                GlobalC::ORB.dr_uniform,
                                                                                GridT.bxyz,
                                                                                GlobalC::ucell.nwmax,
                                                                                psi_input_double_g,
                                                                                psi_input_int_g,
                                                                                num_psir_g,
                                                                                GridT.psi_size_max_per_z,
                                                                                GridT.ucell_atom_nwl_g,
                                                                                GridT.atom_iw2_new_g,
                                                                                GridT.atom_iw2_ylm_g,
                                                                                GridT.atom_nw_g,
                                                                                GridT.nr_max,
                                                                                GridT.psi_u_g,
                                                                                psir_ylm_left_g,
                                                                                psir_ylm_right_g);
            checkCudaLastError();
            GridT.fastest_matrix_mul(max_m, max_n,
                                     atom_pair_A_m_g, atom_pair_B_n_g, atom_pair_k_g,
                                     atom_pair_mat_A_array_g, atom_pair_lda_g,
                                     atom_pair_mat_B_array_g, atom_pair_ldb_g,
                                     atom_pair_mat_C_array_g, atom_pair_ldc_g,
                                     atom_pair_num, GridT.streams[stream_num], nullptr);
            // checkCuda(cudaStreamSynchronize(GridT.streams[stream_num]));
        }
    }
    for (int i = 0; i < GridT.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(GridT.streams[i]));
    }
    {
        int iter_num = 0;
        for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
        {
            for (int iat2 = 0; iat2 < GlobalC::ucell.nat; iat2++)
            {
                int stream_num = iter_num % GridT.nstreams;
                int it1 = GlobalC::ucell.iat2it[iat1];
                int lo1 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
                    it1, GlobalC::ucell.iat2ia[iat1], 0)];

                int it2 = GlobalC::ucell.iat2it[iat2];
                int lo2 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(
                    it2, GlobalC::ucell.iat2ia[iat2], 0)];
                if (lo1 <= lo2) {
                    int atom_pair_nw = GlobalC::ucell.atoms[it1].nw * GlobalC::ucell.atoms[it2].nw;
                    hamilt::AtomPair<double> *tmp_ap = hRGint->find_pair(iat1, iat2);
                    if (tmp_ap == nullptr)
                    {
                        continue;
                    }
                    checkCuda(cudaMemcpyAsync(tmp_ap->get_pointer(0),
                                            GridT.GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2],
                                            atom_pair_nw * sizeof(double),
                                            cudaMemcpyDeviceToHost, GridT.streams[stream_num]));
                    iter_num++;
                }
            }
        }
    }
    for (int i = 0; i < GridT.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(GridT.streams[i]));
    }
}