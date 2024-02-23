#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_rho.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_rho.cuh"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "vbatch_matrix_multiple/vbatch_matrix_mul.cuh"
#include "vbatch_matrix_multiple/cuda_tools.cuh"


void gint_gamma_rho_gpu(hamilt::HContainer<double> *DM,
                        double *rho,
                        const int nczp,
                        const double *ylmcoef_now,
                        const Grid_Technique &GridT)
{
    const int nbz = GridT.nbzp;
    const int lgd = GridT.lgd;
    const int max_size = GridT.max_atom;

    double *dm_matrix_h = new double[lgd * lgd];
    ModuleBase::GlobalFunc::ZEROS(dm_matrix_h, lgd*lgd);
    checkCuda(cudaMemset(GridT.rho_g, 0, GridT.ncxyz * sizeof(double)));
    for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
    {
        for (int iat2 = 0; iat2 < GlobalC::ucell.nat; iat2++)
        {
            int it1 = GlobalC::ucell.iat2it[iat1];
            int it2 = GlobalC::ucell.iat2it[iat2];
            int lo1 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(it1, GlobalC::ucell.iat2ia[iat1], 0)];
            int lo2 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(it2, GlobalC::ucell.iat2ia[iat2], 0)];

                hamilt::AtomPair<double> *tmp_ap = DM->find_pair(iat1, iat2);
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


    for (int i = 0; i < GridT.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(GridT.streams[i]));
    }

    int iter_num = 0;
    for (int i = 0; i < GridT.nbx; i++)
    {
        for (int j = 0; j < GridT.nby; j++)
        {
            int stream_num = iter_num % GridT.nstreams;
            double *psi_input_double = &GridT.psi_input_double_global[GridT.psi_size_max * stream_num * 5];
            int *psi_input_int = &GridT.psi_input_int_global[GridT.psi_size_max * stream_num * 2];
            int *num_psir = &GridT.num_psir_global[nbz * stream_num];
            double *atom_pair_alpha = &GridT.atom_pair_alpha_global[GridT.atom_pair_size_over_nbz * stream_num];
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

            double *atom_pair_alpha_g = &GridT.atom_pair_alpha_global_g[GridT.atom_pair_size_over_nbz * stream_num];
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

            double *rho_g = GridT.rho_g;
            int dot_count = 0;
            int *vec_len = &GridT.vec_len[GridT.num_mcell * stream_num];
            double **vec_l = &GridT.vec_l[GridT.num_mcell * stream_num];
            double **vec_r = &GridT.vec_r[GridT.num_mcell * stream_num];
            double **dot_product = &GridT.dot_product[GridT.num_mcell * stream_num];

            int *vec_len_g = &GridT.vec_len_g[GridT.num_mcell * stream_num];
            double **vec_l_g = &GridT.vec_l_g[GridT.num_mcell * stream_num];
            double **vec_r_g = &GridT.vec_r_g[GridT.num_mcell * stream_num];
            double **dot_product_g = &GridT.dot_product_g[GridT.num_mcell * stream_num];

            int max_m = 0;
            int max_n = 0;
            int atom_pair_num = 0;

            checkCuda(cudaStreamSynchronize(GridT.streams[stream_num]));
            // TODO
            gpu_task_generator_rho(GridT, i, j,
                        GridT.psi_size_max_per_z,
                        max_size, nczp,
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

            checkCuda(cudaMemcpyAsync(psi_input_double_g, psi_input_double, GridT.psi_size_max * 5 * sizeof(double), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(psi_input_int_g, psi_input_int, GridT.psi_size_max * 2 * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(num_psir_g, num_psir, nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_alpha_g, atom_pair_alpha, GridT.atom_pair_size_over_nbz * sizeof(double), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_A_m_g, atom_pair_A_m, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_B_n_g, atom_pair_B_n, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_k_g, atom_pair_k, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_lda_g, atom_pair_lda, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_ldb_g, atom_pair_ldb, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_ldc_g, atom_pair_ldc, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_mat_A_array_g, atom_pair_mat_A_array, GridT.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_mat_B_array_g, atom_pair_mat_B_array, GridT.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_mat_C_array_g, atom_pair_mat_C_array, GridT.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(vec_len_g, vec_len, GridT.num_mcell * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(vec_l_g, vec_l, GridT.num_mcell * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(vec_r_g, vec_r, GridT.num_mcell * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(dot_product_g, dot_product, GridT.num_mcell * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            checkCuda(cudaMemsetAsync(psir_ylm_left_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(psir_ylm_right_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));

            dim3 grid_psi(nbz, 8);
            dim3 block_psi(64);

            get_psi<<<grid_psi, block_psi, 0, GridT.streams[stream_num]>>>(GridT.ylmcoef_g,
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
                                                                            psir_ylm_left_g);
            checkCudaLastError();
            GridT.fastest_matrix_mul(max_m, max_n,
                                     atom_pair_A_m_g, atom_pair_B_n_g, atom_pair_k_g,
                                     atom_pair_mat_A_array_g, atom_pair_lda_g,
                                     atom_pair_mat_B_array_g, atom_pair_ldb_g,
                                     atom_pair_mat_C_array_g, atom_pair_ldc_g,
                                     atom_pair_num, GridT.streams[stream_num], atom_pair_alpha_g);
          

            // new add streamsynchronize
            // checkCuda(cudaStreamSynchronize(GridT.streams[stream_num]));

            dim3 grid_dot(128);
            dim3 block_dot(64);
            int incx = 1;
            int incy = 1;
            psir_dot<<<grid_dot, block_dot, 0, GridT.streams[stream_num]>>>(vec_len_g,
                                                                            vec_l_g,
                                                                            incx,
                                                                            vec_r_g,
                                                                            incy,
                                                                            dot_product_g,
                                                                            dot_count);

            iter_num++;
        }   
    }
    for (int i = 0; i < GridT.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(GridT.streams[i]));
    }
    checkCuda(cudaMemcpy(rho, GridT.rho_g, nczp * GridT.ncx * GridT.ncy * sizeof(double), cudaMemcpyDeviceToHost));
    checkCuda(cudaFree(dm_matrix_g));
    delete[] dm_matrix_h;
}