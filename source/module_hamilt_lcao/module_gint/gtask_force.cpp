#include <omp.h>

#include "gint_force_gpu.h"
#include "module_base/ylm.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
namespace GintKernel
{

/**
 * @brief Description of the function.
 *
 * Detailed description of the function.
 *
 * @param gridt The Grid_Technique object.
 * @param i The integer parameter i.
 * @param j The integer parameter j.
 * @param psiSizeMax The maximum size of psi.
 * @param max_size The maximum size.
 * @param nczp The nczp parameter.
 * @param vfactor The vfactor parameter.
 * @param vlocal_global_value The array of vlocal_global_value.
 * @param iat_per_nbz The array of iat_per_nbz.
 * @param input_dou The double array of psi_input.
 * @param psiInputInt The integer array of psi_input.
 * @param num_psir The array of num_psir.
 * @param lgd The lgd parameter.
 * @param psir_ylm_g The double array of psir_ylm_g.
 * @param psir_zeros_g The double array of psir_zeros_g.
 * @param dm_matrix_g The double array of dm_matrix_g.
 * @param mat_m The array of mat_m.
 * @param mat_n The array of mat_n.
 * @param mat_k The array of mat_k.
 * @param mat_lda The array of mat_lda.
 * @param mat_ldb The array of mat_ldb.
 * @param mat_ldc The array of mat_ldc.
 * @param mat_A The pointer to mat_A.
 * @param mat_B The pointer to mat_B.
 * @param mat_C The pointer to mat_C.
 * @param max_m The reference to max_m.
 * @param max_n The reference to max_n.
 * @param atom_pair_num The reference to atom_pair_num.
 */
void gpu_task_generator_force(const Grid_Technique& gridt,
                              const UnitCell& ucell,
                              const int grid_index_ij,
                              const int psiSizeMax,
                              const int max_size,
                              const int nczp,
                              const double vfactor,
                              const double* rcut,
                              const double* vlocal_global_value,
                              std::vector<int>& iat_per_nbz,
                              int& atom_pair_num,
                              std::vector<bool>& gpu_mat_cal_flag,
                              grid_para& para)
{
    
    const int nwmax = ucell.nwmax;
    // psir generate
    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        int num_get_psi = 0;
        int grid_index = grid_index_ij + z_index;
        int num_psi_pos = psiSizeMax * z_index;
        int calc_flag_index = max_size * z_index;
        int bcell_start_index = gridt.bcell_start[grid_index];
        int na_grid = gridt.how_many_atoms[grid_index];

        for (int id = 0; id < na_grid; id++)
        {
            int ib = 0;
            int mcell_index = bcell_start_index + id;
            int imcell = gridt.which_bigcell[mcell_index];
            int iat = gridt.which_atom[mcell_index];
            int it_temp = ucell.iat2it[iat];
            int start_ind_grid = gridt.start_ind[grid_index];

            for (int bx_index = 0; bx_index < gridt.bx; bx_index++)
            {
                for (int by_index = 0; by_index < gridt.by; by_index++)
                {
                    for (int bz_index = 0; bz_index < gridt.bz; bz_index++)
                    {
                        double dr_temp[3];
                        dr_temp[0] = gridt.meshcell_pos[ib][0]
                                     + gridt.meshball_positions[imcell][0]
                                     - gridt.tau_in_bigcell[iat][0];
                        dr_temp[1] = gridt.meshcell_pos[ib][1]
                                     + gridt.meshball_positions[imcell][1]
                                     - gridt.tau_in_bigcell[iat][1];
                        dr_temp[2] = gridt.meshcell_pos[ib][2]
                                     + gridt.meshball_positions[imcell][2]
                                     - gridt.tau_in_bigcell[iat][2];
                        /* compute distance in and allocate the paramter in
                         * z_index */
                        double distance = sqrt(dr_temp[0] * dr_temp[0]
                                               + dr_temp[1] * dr_temp[1]
                                               + dr_temp[2] * dr_temp[2]);
                        if (distance <= rcut[it_temp])
                        {
                            gpu_mat_cal_flag[calc_flag_index + id] = true;
                            int pos_temp_double = num_psi_pos + num_get_psi;
                            int pos_temp_int = pos_temp_double * 2;
                            pos_temp_double *= 5;
                            if (distance < 1.0E-9)
                            {
                                distance += 1.0E-9;
                            }
                            para.input_dou[pos_temp_double] = dr_temp[0];
                            para.input_dou[pos_temp_double + 1] = dr_temp[1];
                            para.input_dou[pos_temp_double + 2] = dr_temp[2];
                            para.input_dou[pos_temp_double + 3] = distance;
                            int vindex_global = bx_index * gridt.ncy * nczp
                                                + by_index * nczp + bz_index
                                                + start_ind_grid;
                            para.input_dou[pos_temp_double + 4]
                                = vlocal_global_value[vindex_global] * vfactor;

                            para.input_int[pos_temp_int] = it_temp;
                            para.input_int[pos_temp_int + 1]
                                = (z_index * gridt.bxyz + ib) * max_size * nwmax
                                  + id * nwmax;
                            iat_per_nbz[z_index * gridt.bxyz * max_size
                                        + ib * max_size + id]
                                = iat;
                            num_get_psi++;
                        }
                        ib++;
                    }
                }
            }
        }
        para.num_psir[z_index] = num_get_psi;
    }
}

/* allocate the Multiplication of multinomial matrices */
void alloc_mult_force(const Grid_Technique& gridt,
                                    const UnitCell& ucell,
                                    const int grid_index_ij,
                                    const int max_size,
                                    const int lgd,
                                    double* dm_matrix_g,
                                    int& max_m,
                                    int& max_n,
                                    int& atom_pair_num,
                                    std::vector<bool>& gpu_mat_cal_flag,
                                    grid_para& para)
{
    int tid = 0;
    max_m = 0;
    max_n = 0;
    const int nwmax=ucell.nwmax;
    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        int grid_index = grid_index_ij + z_index;
        int calc_flag_index = max_size * z_index;
        int bcell_start_index = gridt.bcell_start[grid_index];
        int bcell_start_psir = z_index * gridt.bxyz * max_size * nwmax;

        for (int atom1 = 0; atom1 < gridt.how_many_atoms[grid_index]; atom1++)
        {
            if (!gpu_mat_cal_flag[calc_flag_index + atom1])
            {
                continue;
            }
            const int mcell_index1 = bcell_start_index + atom1;
            int iat1 = gridt.which_atom[mcell_index1];
            int it1 = ucell.iat2it[iat1];
            int lo1
                = gridt.trace_lo[ucell.itiaiw2iwt(it1, ucell.iat2ia[iat1], 0)];
            int nw1 = ucell.atoms[it1].nw;

            for (int atom2 = 0; atom2 < gridt.how_many_atoms[grid_index];atom2++)
            {
                if (!gpu_mat_cal_flag[calc_flag_index + atom2])
                {
                    continue;
                }
                const int mcell_index2 = bcell_start_index + atom2;
                int iat2 = gridt.which_atom[mcell_index2];
                int it2 = ucell.iat2it[iat2];
                int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(it2,
                                                            ucell.iat2ia[iat2],
                                                            0)];
                int nw2 = ucell.atoms[it2].nw;

                int mat_A_idx = bcell_start_psir + atom2 * nwmax;
                int mat_B_idx = lgd * lo1 + lo2;
                int mat_C_idx = bcell_start_psir + atom1 * nwmax;
                para.atom_pair_A_m[tid] = gridt.bxyz;
                para.atom_pair_B_n[tid] = nw1;
                para.atom_pair_K[tid] = nw2;
                para.atom_pair_lda[tid] = nwmax * max_size;
                para.atom_pair_ldb[tid] = lgd;
                para.atom_pair_ldc[tid] = nwmax * max_size;
                para.matrix_A[tid] = para.psir_r_device + mat_A_idx;
                para.matrix_B[tid] = dm_matrix_g + mat_B_idx;
                para.matrix_C[tid] = para.psir_dm_device + mat_C_idx;

                if (para.atom_pair_A_m[tid] > max_m)
                {
                    max_m = para.atom_pair_A_m[tid];
                }

                if (para.atom_pair_B_n[tid] > max_n)
                {
                    max_n = para.atom_pair_B_n[tid];
                }

                tid++;
            }
        }
    }
    atom_pair_num = tid;
}


void allocateDm(std::vector<double> &matrixHost,
                hamilt::HContainer<double>* dm,
                const Grid_Technique& gridt,
                const UnitCell& ucell)
{
    matrixHost = std::vector<double>(gridt.lgd * gridt.lgd, 0);
    for (int iatRow = 0; iatRow < ucell.nat; iatRow++)
    {
        for (int iatColumn = 0; iatColumn < ucell.nat; iatColumn++)
        {
            int indexTypeRow = ucell.iat2it[iatRow];
            int indexTypeColumn = ucell.iat2it[iatColumn];
            int localOrbitRow
                = gridt.trace_lo[ucell.itiaiw2iwt(indexTypeRow,
                                                  ucell.iat2ia[iatRow],
                                                  0)];
            int localOrbitColumn
                = gridt.trace_lo[ucell.itiaiw2iwt(indexTypeColumn,
                                                  ucell.iat2ia[iatColumn],
                                                  0)];
            hamilt::AtomPair<double>* tmpAtomPair
                = dm->find_pair(iatRow, iatColumn);
            int orbitIndex = 0;
            if (tmpAtomPair == NULL)
            {
                continue;
            }
            for (int orbitRow = 0; orbitRow < tmpAtomPair->get_row_size();
                 orbitRow++)
            {
                for (int orbitColumn = 0;
                     orbitColumn < tmpAtomPair->get_col_size();
                     orbitColumn++)
                {
                    matrixHost[(localOrbitRow + orbitRow) * gridt.lgd
                               + (localOrbitColumn + orbitColumn)]
                        = tmpAtomPair->get_pointer(0)[orbitIndex];
                    orbitIndex++;
                }
            }
        }
    }
    return;
}
void calculateInit(DensityMat& denstiy_mat,
                   frc_strs_iat_gbl& f_s_iat_dev,
                   hamilt::HContainer<double>* dm,
                   const Grid_Technique& gridt,
                   const UnitCell& ucell,
                   const int lgd,
                   const int cuda_block,
                   const int nat,
                   const int atom_num_grid)
{
    denstiy_mat.density_mat_h = std::vector<double>(lgd * lgd, 0);
    allocateDm(denstiy_mat.density_mat_h, dm, gridt, ucell);

    checkCuda(cudaMalloc((void**)&denstiy_mat.density_mat_d,
                         lgd * lgd * sizeof(double)));
    checkCuda(cudaMemcpy(denstiy_mat.density_mat_d,
                         denstiy_mat.density_mat_h.data(),
                         lgd * lgd * sizeof(double),
                         cudaMemcpyHostToDevice));

    checkCuda(cudaMalloc((void**)&f_s_iat_dev.stress_global,
                         6  * gridt.nstreams * sizeof(double)));
    checkCuda(cudaMemset(f_s_iat_dev.stress_global,
                         0,
                         6  * gridt.nstreams * sizeof(double)));

    checkCuda(cudaMalloc((void**)&f_s_iat_dev.force_global,
                         3 * nat * gridt.nstreams * sizeof(double)));
    checkCuda(cudaMemset(f_s_iat_dev.force_global,
                         0,
                         3 * nat * gridt.nstreams * sizeof(double)));

    checkCuda(cudaMalloc((void**)&f_s_iat_dev.iat_global,
                         atom_num_grid * gridt.nstreams * sizeof(int)));
    checkCuda(cudaMemset(f_s_iat_dev.iat_global,
                         0,
                         atom_num_grid * gridt.nstreams * sizeof(int)));
}

/**
 * @brief grid parameter Init
 *
 * GridParameter init
 *
 * @param para double *,contained the destiyMatHost
 * @param iter_num int , used for calcute the stream
 * @param nbz int,stand for the number of Z-axis
 * @param gridt Grid_Technique,stored the major method in the the gint.
 */
void para_init(grid_para& para,
                       const int iter_num,
                       const int nbz,
                       const int pipeline_index,
                       const Grid_Technique& gridt)
{

    // pipeline_index stand for nstreams
    
    //input_dou and input _int used for the Spherical Harmonics
    para.input_dou
        = &gridt.psi_dbl_gbl[gridt.psi_size_max * pipeline_index * 5];
    para.input_int
        = &gridt.psi_int_gbl[gridt.psi_size_max * pipeline_index * 2];
    para.num_psir = &gridt.num_psir_gbl[nbz * pipeline_index];
    //one dimension,record the length and the leading dimension of three matrix
    para.atom_pair_A_m
        = &gridt.l_info_global[gridt.atom_pair_nbz * pipeline_index];
    para.atom_pair_B_n
        = &gridt.r_info_global[gridt.atom_pair_nbz * pipeline_index];
    para.atom_pair_K
        = &gridt.k_info_global[gridt.atom_pair_nbz * pipeline_index];
    para.atom_pair_lda
        = &gridt.lda_info_global[gridt.atom_pair_nbz * pipeline_index];
    para.atom_pair_ldb
        = &gridt.ldb_info_global[gridt.atom_pair_nbz * pipeline_index];
    para.atom_pair_ldc
        = &gridt.ldc_info_global[gridt.atom_pair_nbz * pipeline_index];
    //input_double_g and input_int_g used for the Spherical Harmonics on GPU
    para.input_double_g
        = &gridt.psi_dbl_gbl_g[gridt.psi_size_max * pipeline_index * 5];
    para.input_int_g
        = &gridt.psi_int_gbl_g[gridt.psi_size_max * pipeline_index * 2];
    para.num_psir_g = &gridt.num_psir_gbl_g[nbz * pipeline_index];
    para.psir_dm_device = &gridt.dm_global_g[gridt.psir_size * pipeline_index];
    para.psir_r_device
        = &gridt.right_global_g[gridt.psir_size * pipeline_index];
    //psi function ,record the force in x y z,and the stress in six dimension
    para.psir_lx_device = &gridt.d_left_x_g[gridt.psir_size * pipeline_index];
    para.psir_ly_device = &gridt.d_left_y_g[gridt.psir_size * pipeline_index];
    para.psir_lz_device = &gridt.d_left_z_g[gridt.psir_size * pipeline_index];
    para.psir_lxx_device
        = &gridt.dd_left_xx_g[gridt.psir_size * pipeline_index];
    para.psir_lxy_device
        = &gridt.dd_left_xy_g[gridt.psir_size * pipeline_index];
    para.psir_lxz_device
        = &gridt.dd_left_xz_g[gridt.psir_size * pipeline_index];
    para.psir_lyy_device
        = &gridt.dd_left_yy_g[gridt.psir_size * pipeline_index];
    para.psir_lyz_device
        = &gridt.dd_left_yz_g[gridt.psir_size * pipeline_index];
    para.psir_lzz_device
        = &gridt.dd_left_zz_g[gridt.psir_size * pipeline_index];
    //one dimension,record the length and the leading dimension of three matrix on GPU
    para.A_m_device
        = &gridt.l_info_global_g[gridt.atom_pair_nbz * pipeline_index];
    para.B_n_device
        = &gridt.r_info_global_g[gridt.atom_pair_nbz * pipeline_index];
    para.K_device
        = &gridt.k_info_global_g[gridt.atom_pair_nbz * pipeline_index];
    para.lda_device
        = &gridt.lda_info_gbl_g[gridt.atom_pair_nbz * pipeline_index];
    para.ldb_device
        = &gridt.ldb_info_gbl_g[gridt.atom_pair_nbz * pipeline_index];
    para.ldc_device
        = &gridt.ldc_info_gbl_g[gridt.atom_pair_nbz * pipeline_index];
    //two dimension,record number to compute
    para.matrix_A = &gridt.ap_left_gbl[gridt.atom_pair_nbz * pipeline_index];
    para.matrix_B = &gridt.ap_right_gbl[gridt.atom_pair_nbz * pipeline_index];
    para.matrix_C = &gridt.ap_output_gbl[gridt.atom_pair_nbz * pipeline_index];
    para.matrix_A_device
        = &gridt.ap_left_gbl_g[gridt.atom_pair_nbz * pipeline_index];
    para.matrix_B_device
        = &gridt.ap_right_gbl_g[gridt.atom_pair_nbz * pipeline_index];
    para.matrix_C_device
        = &gridt.ap_output_gbl_g[gridt.atom_pair_nbz * pipeline_index];
}
/**
 * @brief frc_strs_iat on host and device Init
 *
 * GridParameter init
 *
 * @param frc_strs_iat frc_strs_iat,contains the Force Stree Iat on Host
 * @param pipeline_index int , record the stream in GPU
 * @param cuda_block in stress compute,used for Block nums
 * @param atom_num_grid in force calculate,used for Block nums
 * @param max_size Maximum size of atoms on a grid.
 * @param frc_strs_iat_gbl frc_strs_iat_gbl,contains the Force Stree Iat on Host
 */
void cal_init(frc_strs_iat& f_s_iat,
                const int pipeline_index,
                const int cuda_block,
                const int atom_num_grid,
                const int nat,
                const int max_size,
                frc_strs_iat_gbl& f_s_iat_dev)
{
    const int iat_min = -max_size - 1;
    f_s_iat.stress_device
        = &f_s_iat_dev.stress_global[6 * pipeline_index];
    f_s_iat.force_device
        = &f_s_iat_dev.force_global[3 * nat * pipeline_index];
    f_s_iat.iat_device
        = &f_s_iat_dev.iat_global[atom_num_grid * pipeline_index];
    f_s_iat.iat_host = vector<int>(atom_num_grid, iat_min);
}





/**
 * @brief GridParameter memCpy,from Host to Device
 *
 * parameter init,which contains the gpu task and multi matrix multiplication
 *
 * @param para Grid parameter in task generator.
 * @param f_s_iat frc_strs_iat,contains the Force Stree Iat.
 * @param gridt Grid_Technique,stored the major method in the the gint.
 * @param nbz int,stand for the number of Z-axis
 * @param atom_num_grid in force calculate,used for Block nums
 */
void mem_copy(grid_para& para,
                    frc_strs_iat& f_s_iat,
                    const Grid_Technique& gridt,
                    const int nbz,
                    const int pipeline_index,
                    const int atom_num_grid)
{
    checkCuda(cudaMemcpyAsync(para.input_double_g,
                              para.input_dou,
                              gridt.psi_size_max * 5 * sizeof(double),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.input_int_g,
                              para.input_int,
                              gridt.psi_size_max * 2 * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.num_psir_g,
                              para.num_psir,
                              nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.A_m_device,
                              para.atom_pair_A_m,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.B_n_device,
                              para.atom_pair_B_n,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.K_device,
                              para.atom_pair_K,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.lda_device,
                              para.atom_pair_lda,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.ldb_device,
                              para.atom_pair_ldb,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.ldc_device,
                              para.atom_pair_ldc,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.matrix_A_device,
                              para.matrix_A,
                              gridt.atom_pair_nbz * sizeof(double*),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.matrix_B_device,
                              para.matrix_B,
                              gridt.atom_pair_nbz * sizeof(double*),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(para.matrix_C_device,
                              para.matrix_C,
                              gridt.atom_pair_nbz * sizeof(double*),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_dm_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_r_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_lx_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_ly_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_lz_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_lxx_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_lxy_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_lxz_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_lyy_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_lyz_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemsetAsync(para.psir_lzz_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[pipeline_index]));
    checkCuda(cudaMemcpyAsync(f_s_iat.iat_device,
                              f_s_iat.iat_host.data(),
                              atom_num_grid * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[pipeline_index]));
}

} // namespace GintKernel
