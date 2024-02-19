#include <fstream>
#include <sstream>

#include "gint_force.cuh"
#include "gint_force.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "omp.h"
#include "vbatch_matrix_multiple/cuda_tools.cuh"

// Function to calculate forces using GPU-accelerated gamma point Gint
void gint_gamma_force_gpu(hamilt::HContainer<double> *DM,
                        const double vfactor,
                        const double *vlocal,
                        double *force,
                        double *stress,
                        const int nczp,
                        const double *ylmcoef_now,
                        const Grid_Technique &GridT)
{
    const int nbz = GridT.nbzp;
    const int lgd = GridT.lgd;
    const int max_size = GridT.max_atom;
    const int nat =GlobalC::ucell.nat;
    const int nwmax =GlobalC::ucell.nwmax;
    const int bxyz=GridT.bxyz;
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

    const int threadsPerBlock=256;
    const int blocksPerGrid = std::min(64, (GridT.psir_size + threadsPerBlock - 1) / threadsPerBlock);
    const int threadsPerBlock_force=256;
    const int blocksPerGrid_force =std::min(64,(nbz*bxyz*max_size+threadsPerBlock_force-1)/threadsPerBlock_force);
    double *stress_dot=new double[6*blocksPerGrid];

    for (int i=0;i<6*blocksPerGrid;i++)
    {
        stress_dot[i]=0.0;
    }
    double *stress_dot_global_g;
    checkCuda(cudaMalloc((void **)&stress_dot_global_g, 6*blocksPerGrid * GridT.nstreams * sizeof(double)));
    checkCuda(cudaMemset(stress_dot_global_g, 0, 6*blocksPerGrid * GridT.nstreams * sizeof(double)));


    double *force_dot_global_g;
    checkCuda(cudaMalloc((void **)&force_dot_global_g,3*nbz*bxyz*max_size* GridT.nstreams * sizeof(double)));
    checkCuda(cudaMemset(force_dot_global_g,0,3*nbz*bxyz*max_size*GridT.nstreams*sizeof(double)));

    

    int *iat_global_g;
    checkCuda(cudaMalloc((void **)&iat_global_g,nbz*bxyz*max_size*GridT.nstreams*sizeof(int)));
    checkCuda(cudaMemset(iat_global_g,0,nbz*bxyz*max_size*GridT.nstreams*sizeof(int)));

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
            int *atom_pair_A_m = &GridT.atom_pair_left_info_global[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_B_n = &GridT.atom_pair_right_info_global[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_k = &GridT.atom_pair_k_info_global[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_lda = &GridT.atom_pair_lda_global[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_ldb = &GridT.atom_pair_ldb_global[GridT.atom_pair_size_over_nbz * stream_num];
            int *atom_pair_ldc = &GridT.atom_pair_ldc_global[GridT.atom_pair_size_over_nbz * stream_num];

            double *psi_input_double_g = &GridT.psi_input_double_global_g[GridT.psi_size_max * stream_num * 5];
            int *psi_input_int_g = &GridT.psi_input_int_global_g[GridT.psi_size_max * stream_num * 2];
            int *num_psir_g = &GridT.num_psir_global_g[nbz * stream_num];
            double *psir_ylm_dm_g = &GridT.psir_ylm_dm_global_g[GridT.psir_size * stream_num];
            double *psir_ylm_right_g = &GridT.psir_ylm_right_global_g[GridT.psir_size * stream_num];
            double *dpsir_ylm_left_x_g = &GridT.dpsir_ylm_left_x_global_g[GridT.psir_size * stream_num];
            double *dpsir_ylm_left_y_g = &GridT.dpsir_ylm_left_y_global_g[GridT.psir_size * stream_num];
            double *dpsir_ylm_left_z_g = &GridT.dpsir_ylm_left_z_global_g[GridT.psir_size * stream_num];
            double *ddpsir_ylm_left_xx_g = &GridT.ddpsir_ylm_left_xx_global_g[GridT.psir_size * stream_num];
            double *ddpsir_ylm_left_xy_g = &GridT.ddpsir_ylm_left_xy_global_g[GridT.psir_size * stream_num];
            double *ddpsir_ylm_left_xz_g = &GridT.ddpsir_ylm_left_xz_global_g[GridT.psir_size * stream_num];
            double *ddpsir_ylm_left_yy_g = &GridT.ddpsir_ylm_left_yy_global_g[GridT.psir_size * stream_num];
            double *ddpsir_ylm_left_yz_g = &GridT.ddpsir_ylm_left_yz_global_g[GridT.psir_size * stream_num];
            double *ddpsir_ylm_left_zz_g = &GridT.ddpsir_ylm_left_zz_global_g[GridT.psir_size * stream_num];

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

            double *stress_dot_g =&stress_dot_global_g[6*blocksPerGrid*stream_num];
            double *force_dot_g =&force_dot_global_g[3*nbz*bxyz*max_size*stream_num];

            int *iat_g=&iat_global_g[nbz*bxyz*max_size*stream_num];
            int *iat=new int[nbz*bxyz*max_size];
            for (int index=0;index<nbz*bxyz*max_size;index++)
            {
                iat[index]=-max_size-1;
            }
            double *force_h=new double[3*nbz*bxyz*max_size];
            for (int index=0;index<3*nbz*bxyz*max_size;index++)
            {
                force_h[index]=0.0;
            }

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
            gpu_task_generator_force(GridT, i, j,
                        GridT.psi_size_max_per_z,
                        max_size, nczp,
                        vfactor,
                        vlocal,
                        iat,
                        psi_input_double,
                        psi_input_int,
                        num_psir,
                        lgd,
                        psir_ylm_right_g,
                        psir_ylm_dm_g,
                        dm_matrix_g,
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
            checkCuda(cudaMemcpyAsync(iat_g,iat,nbz*bxyz*max_size*sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_A_m_g, atom_pair_A_m, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_B_n_g, atom_pair_B_n, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_k_g, atom_pair_k, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_lda_g, atom_pair_lda, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_ldb_g, atom_pair_ldb, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_ldc_g, atom_pair_ldc, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_mat_A_array_g, atom_pair_mat_A_array, GridT.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_mat_B_array_g, atom_pair_mat_B_array, GridT.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_mat_C_array_g, atom_pair_mat_C_array, GridT.atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            // checkCuda(cudaMemcpyAsync(vec_len_g, vec_len, GridT.num_mcell * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            // checkCuda(cudaMemcpyAsync(vec_l_g, vec_l, GridT.num_mcell * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            // checkCuda(cudaMemcpyAsync(vec_r_g, vec_r, GridT.num_mcell * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            // checkCuda(cudaMemcpyAsync(dot_product_g, dot_product, GridT.num_mcell * sizeof(double *), cudaMemcpyHostToDevice, GridT.streams[stream_num]));

            checkCuda(cudaMemsetAsync(psir_ylm_dm_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(psir_ylm_right_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(dpsir_ylm_left_x_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(dpsir_ylm_left_y_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(dpsir_ylm_left_z_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(ddpsir_ylm_left_xx_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(ddpsir_ylm_left_xy_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(ddpsir_ylm_left_xz_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(ddpsir_ylm_left_yy_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(ddpsir_ylm_left_yz_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(ddpsir_ylm_left_zz_g, 0, GridT.psir_size * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(stress_dot_g, 0, 6*blocksPerGrid * sizeof(double), GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(force_dot_g,0, 3*nbz*bxyz*max_size*sizeof(double), GridT.streams[stream_num]));
            dim3 grid_psi(nbz, 8);
            dim3 block_psi(64);

            get_psi_force<<<grid_psi, block_psi, 0, GridT.streams[stream_num]>>>(GridT.ylmcoef_g,
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
                                                                            GridT.atom_iw2_l_g,
                                                                            GridT.atom_nw_g,
                                                                            GridT.nr_max,
                                                                            GridT.psi_u_g,
                                                                            psir_ylm_right_g,
                                                                            dpsir_ylm_left_x_g,
                                                                            dpsir_ylm_left_y_g,
                                                                            dpsir_ylm_left_z_g,
                                                                            ddpsir_ylm_left_xx_g,
                                                                            ddpsir_ylm_left_xy_g,
                                                                            ddpsir_ylm_left_xz_g,
                                                                            ddpsir_ylm_left_yy_g,
                                                                            ddpsir_ylm_left_yz_g,
                                                                            ddpsir_ylm_left_zz_g);
            checkCudaLastError();
            GridT.fastest_matrix_mul(max_m, max_n,
                                     atom_pair_A_m_g, atom_pair_B_n_g, atom_pair_k_g,
                                     atom_pair_mat_A_array_g, atom_pair_lda_g,
                                     atom_pair_mat_B_array_g, atom_pair_ldb_g,
                                     atom_pair_mat_C_array_g, atom_pair_ldc_g,
                                     atom_pair_num, GridT.streams[stream_num],nullptr);

            dim3 grid_dot_force(blocksPerGrid_force);
            dim3 block_dot_force(threadsPerBlock_force);
            dim3 grid_dot(blocksPerGrid);
            dim3 block_dot(threadsPerBlock);
            dot_product_force<<<grid_dot_force,block_dot_force,0,GridT.streams[stream_num]>>>(dpsir_ylm_left_x_g,
                                                                                         dpsir_ylm_left_y_g,
                                                                                         dpsir_ylm_left_z_g,
                                                                                         psir_ylm_dm_g,
                                                                                         force_dot_g,
                                                                                         iat_g,
                                                                                         nwmax,
                                                                                         max_size,
                                                                                         GridT.psir_size/nwmax);

            checkCuda(cudaMemcpy(force_h, force_dot_g ,3*nbz*bxyz*max_size*sizeof(double), cudaMemcpyDeviceToHost));
            for (int index1=0;index1<nbz*bxyz*max_size;index1++)
            {
                int iat1=iat[index1];
                if (iat1>=0)
                {
                for (int index2=0;index2<3;index2++)
                    {    
                        force[iat1*3+index2]+=force_h[index1*3+index2];
                    }
                }
            }
            dot_product_stress<<<grid_dot, block_dot, 0, GridT.streams[stream_num]>>>(ddpsir_ylm_left_xx_g,
                                                                                                ddpsir_ylm_left_xy_g,
                                                                                                ddpsir_ylm_left_xz_g,
                                                                                                ddpsir_ylm_left_yy_g,
                                                                                                ddpsir_ylm_left_yz_g,
                                                                                                ddpsir_ylm_left_zz_g,
                                                                                                psir_ylm_dm_g,
                                                                                                stress_dot_g,
                                                                                                GridT.psir_size);
            
            checkCuda(cudaMemcpy(stress_dot, stress_dot_g ,6*blocksPerGrid * sizeof(double), cudaMemcpyDeviceToHost));
            for (int i=0;i<6;i++)
            {
                for (int index=0;index<blocksPerGrid;index++)
                {
                    stress[i]+=stress_dot[i*blocksPerGrid+index];
                }
            }

            iter_num++;
        }   
    }
    checkCuda(cudaFree(stress_dot_global_g));
    for (int i = 0; i < GridT.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(GridT.streams[i]));
    }
}