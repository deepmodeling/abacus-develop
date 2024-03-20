#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_rho.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_rho.cuh"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "vbatch_matrix_multiple/vbatch_matrix_mul.cuh"
#include "vbatch_matrix_multiple/cuda_tools.cuh"

namespace lcaoCudaKernel{
void gint_gamma_rho_gpu(hamilt::HContainer<double> *dm,
                        double *rho,
                        const int nczp,
                        const double *ylmcoef_now,
                        const LCAO_Orbitals &ORB,
                        const Grid_Technique &gridt,
                        const UnitCell &ucell)
{
    const int nbz = gridt.nbzp;
    const int lgd = gridt.lgd;
    const int max_size = gridt.max_atom;

    double *dm_matrix_h = new double[lgd * lgd];
    ModuleBase::GlobalFunc::ZEROS(dm_matrix_h, lgd*lgd);
    checkCuda(cudaMemset(gridt.rho_g, 0, gridt.ncxyz * sizeof(double)));
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        for (int iat2 = 0; iat2 < ucell.nat; iat2++)
        {
            int it1 = ucell.iat2it[iat1];
            int it2 = ucell.iat2it[iat2];
            int lo1 = gridt.trace_lo[ucell.itiaiw2iwt(it1, ucell.iat2ia[iat1], 0)];
            int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(it2, ucell.iat2ia[iat2], 0)];

                hamilt::AtomPair<double> *tmp_ap = dm->find_pair(iat1, iat2);
                int orb_index = 0;
                if (tmp_ap == NULL)
                {
                    continue; 
                }
                for (int orb_i = 0; orb_i < tmp_ap->get_row_size(); orb_i++)
                {
                    for (int orb_j = 0; orb_j < tmp_ap->get_col_size(); orb_j++)
                    {
                        
                        dm_matrix_h[(lo1 + orb_i) * lgd + (lo2 + orb_j)]=tmp_ap->get_pointer(0)[orb_index];
                        orb_index++;
                    }
                }
        }
    }
    double *dm_matrix_g;
    checkCuda(cudaMalloc((void **)&dm_matrix_g,lgd * lgd * sizeof(double)));
    checkCuda(cudaMemcpy(dm_matrix_g,dm_matrix_h,lgd * lgd *sizeof(double),cudaMemcpyHostToDevice));


    for (int i = 0; i < gridt.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(gridt.streams[i]));
    }

    int iter_num = 0;
    for (int i = 0; i < gridt.nbx; i++)
    {
        for (int j = 0; j < gridt.nby; j++)
        {
            int stream_num = iter_num % gridt.nstreams;
            double *psi_input_double = &gridt.psi_input_double_global[gridt.psi_size_max * stream_num * 5];
            int *psi_input_int = &gridt.psi_input_int_global[gridt.psi_size_max * stream_num * 2];
            int *num_psir = &gridt.num_psir_global[nbz * stream_num];
            double *atom_pair_alpha = &gridt.atom_pair_alpha_global[gridt.atom_pair_size_over_nbz * stream_num];
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

            double *atom_pair_alpha_g = &gridt.atom_pair_alpha_global_g[gridt.atom_pair_size_over_nbz * stream_num];
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

            double *rho_g = gridt.rho_g;
            int dot_count = 0;
            int *vec_len = &gridt.vec_len[gridt.num_mcell * stream_num];
            double **vec_l = &gridt.vec_l[gridt.num_mcell * stream_num];
            double **vec_r = &gridt.vec_r[gridt.num_mcell * stream_num];
            double **dot_product = &gridt.dot_product[gridt.num_mcell * stream_num];

            int *vec_len_g = &gridt.vec_len_g[gridt.num_mcell * stream_num];
            double **vec_l_g = &gridt.vec_l_g[gridt.num_mcell * stream_num];
            double **vec_r_g = &gridt.vec_r_g[gridt.num_mcell * stream_num];
            double **dot_product_g = &gridt.dot_product_g[gridt.num_mcell * stream_num];

            int max_m = 0;
            int max_n = 0;
            int atom_pair_num = 0;

            checkCuda(cudaStreamSynchronize(gridt.streams[stream_num]));
            // TODO
            gpu_task_generator_rho(gridt, i, j,
                        gridt.psi_size_max_per_z,
                        max_size, nczp,
                        ucell, ORB,
                        psi_input_double,
                        psi_input_int,
                        num_psir,
                        lgd,
                        psir_ylm_left_g,
                        psir_ylm_right_g,
                        dm_matrix_g,
                        atom_pair_alpha,
                        atom_pair_A_m,
                        atom_pair_B_n,
                        atom_pair_k,
                        atom_pair_lda,
                        atom_pair_ldb,
                        atom_pair_ldc,
                        atom_pair_mat_A_array,
                        atom_pair_mat_B_array,
                        atom_pair_mat_C_array,
                        max_m,
                        max_n,
                        atom_pair_num,
                        rho_g,
                        vec_l,
                        vec_r,
                        dot_product,
                        vec_len,
                        dot_count
                        );

            checkCuda(cudaMemcpyAsync(psi_input_double_g, psi_input_double, gridt.psi_size_max * 5 * sizeof(double), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(psi_input_int_g, psi_input_int, gridt.psi_size_max * 2 * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(num_psir_g, num_psir, nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_alpha_g, atom_pair_alpha, gridt.atom_pair_size_over_nbz * sizeof(double), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_A_m_g, atom_pair_A_m, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_B_n_g, atom_pair_B_n, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_k_g, atom_pair_k, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_lda_g, atom_pair_lda, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_ldb_g, atom_pair_ldb, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_ldc_g, atom_pair_ldc, gridt.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_mat_A_array_g, atom_pair_mat_A_array, gridt.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_mat_B_array_g, atom_pair_mat_B_array, gridt.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_mat_C_array_g, atom_pair_mat_C_array, gridt.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, gridt.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(vec_len_g, vec_len, gridt.num_mcell * sizeof(int), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(vec_l_g, vec_l, gridt.num_mcell * sizeof(double *), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(vec_r_g, vec_r, gridt.num_mcell * sizeof(double *), cudaMemcpyHostToDevice, gridt.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(dot_product_g, dot_product, gridt.num_mcell * sizeof(double *), cudaMemcpyHostToDevice, gridt.streams[stream_num]));

            checkCuda(cudaMemsetAsync(psir_ylm_left_g, 0, gridt.psir_size * sizeof(double), gridt.streams[stream_num]));
            checkCuda(cudaMemsetAsync(psir_ylm_right_g, 0, gridt.psir_size * sizeof(double), gridt.streams[stream_num]));

            dim3 grid_psi(nbz, 8);
            dim3 block_psi(64);

            get_psi<<<grid_psi, block_psi, 0, gridt.streams[stream_num]>>>(gridt.ylmcoef_g,
                                                                            ORB.dr_uniform,
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
                                                                            psir_ylm_left_g);
            checkCudaLastError();
            gridt.fastest_matrix_mul(max_m, max_n,
                                     atom_pair_A_m_g, atom_pair_B_n_g, atom_pair_k_g,
                                     atom_pair_mat_A_array_g, atom_pair_lda_g,
                                     atom_pair_mat_B_array_g, atom_pair_ldb_g,
                                     atom_pair_mat_C_array_g, atom_pair_ldc_g,
                                     atom_pair_num, gridt.streams[stream_num], atom_pair_alpha_g);
          
            dim3 grid_dot(128);
            dim3 block_dot(64);
            int incx = 1;
            int incy = 1;
            psir_dot<<<grid_dot, block_dot, 0, gridt.streams[stream_num]>>>(vec_len_g,
                                                                            vec_l_g,
                                                                            incx,
                                                                            vec_r_g,
                                                                            incy,
                                                                            dot_product_g,
                                                                            dot_count);

            iter_num++;
        }   
    }
    for (int i = 0; i < gridt.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(gridt.streams[i]));
    }
    checkCuda(cudaMemcpy(rho, gridt.rho_g, nczp * gridt.ncx * gridt.ncy * sizeof(double), cudaMemcpyDeviceToHost));
    checkCuda(cudaFree(dm_matrix_g));
    delete[] dm_matrix_h;
}
}