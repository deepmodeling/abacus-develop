#include "gint_force.cuh"
#include "interp.cuh"
#include "module_hamilt_lcao/module_gint/gint_force.h"
#include "module_hamilt_lcao/module_gint/kernels/cuda/cuda_tools.cuh"
#include "module_hamilt_lcao/module_gint/kernels/cuda/gint_force.cuh"
#include "module_hamilt_lcao/module_gint/kernels/cuda/sph.cuh"
// CUDA kernel to calculate psi and force
namespace GintKernel
{

/*!
 * \file
 * \brief CUDA kernel to calculate psi and force
 *
 * CUDA kernel that performs calculations on psi and force.
 *
 * \param ylmcoef Pointer to the Ylm coefficients
 * \param delta_r_g Delta r value
 * \param bxyz_g Bxyz value
 * \param nwmax_g Nwmax value
 * \param input_double Array of double input values
 * \param input_int Array of int input values
 * \param num_psir Array containing the number of psi for each block
 * \param psi_size_max Maximum size of psi
 * \param ucell_atom_nwl Array containing Ucell atom nwl values
 * \param atom_iw2_new Array indicating whether atom_iw2 is new
 * \param atom_iw2_ylm Array of atom_iw2 Ylm values
 * \param atom_iw2_l Array of atom_iw2 l values
 * \param atom_nw Array of atom_nw values
 * \param nr_max Maximum nr value
 * \param psi_u Array for psi_u values,recording the Spherical Harmonics from psi
 * \param psir_r Array for psir_r values,recored the distance from psi
 * \param psir_lx Array for psir_lx values,recored the force left in x
 * \param psir_ly Array for psir_ly values,recored the force left in y
 * \param psir_lz Array for psir_lz values,recored the force left in z
 * \param psir_lxx Array for psir_lxx values,recored the stress left in xx
 * \param psir_lxy Array for psir_lxy values,recored the stress left in xy
 * \param psir_lxz Array for psir_lxz values,recored the stress left in xz
 * \param psir_lyy Array for psir_lyy values,recored the stress left in yy
 * \param psir_lyz Array for psir_lyz values,recored the stress left in yz
 * \param psir_lzz Array for psir_lzz values,recored the stress left in zz
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
                              double* psir_lzz)
{
    // Get the size of psi for the current block
    int size = num_psir[blockIdx.x];
    int start_index = psi_size_max * blockIdx.x;
    int end_index = start_index + size;
    start_index += threadIdx.x + blockDim.x * blockIdx.y;
    // Loop over the psi indices for the current block
    for (int index = start_index; index < end_index;
         index += blockDim.x * gridDim.y)
    {
        // Extract information from input arrays
        double dr[3];
        int index_double = index * 5;
        dr[0] = input_double[index_double];
        dr[1] = input_double[index_double + 1];
        dr[2] = input_double[index_double + 2];
        double distance = input_double[index_double + 3];
        distance = distance * distance;
        double vlbr3_value = input_double[index_double + 4];
        // begin calculation
        double ylma[49]; // Attention!!! At present, we only use L=5 at
                         // most. So (L+1) * (L+1)=36
        double grly[49][3];
        int index_int = index * 2;
        int it = input_int[index_int];
        int dist_tmp = input_int[index_int + 1];

        int nwl = ucell_atom_nwl[it];
        spherical_harmonics_d(dr, distance, grly, nwl, ylma, ylmcoef);

        interpolate_f(distance,
                      delta_r_g,
                      it,
                      nwmax_g,
                      nr_max,
                      atom_nw,
                      atom_iw2_new,
                      psi_u,
                      atom_iw2_l,
                      atom_iw2_ylm,
                      psir_r,
                      dist_tmp,
                      ylma,
                      vlbr3_value,
                      psir_lx,
                      dr,
                      grly,
                      psir_ly,
                      psir_lz,
                      psir_lxx,
                      psir_lxy,
                      psir_lxz,
                      psir_lyy,
                      psir_lyz,
                      psir_lzz);
    }
}


/**
 * \brief Compute dot product of stress components and partial derivatives.
 *
 * This CUDA kernel computes the dot product of stress components and partial
 * derivatives based on the input arrays.
 *
 * \param psir_lxx Array of psir_lxx values.
 * \param psir_lxy Array of psir_lxy values.
 * \param psir_lxz Array of psir_lxz values.
 * \param psir_lyy Array of psir_lyy values.
 * \param psir_lyz Array of psir_lyz values.
 * \param psir_lzz Array of psir_lzz values.
 * \param psir_ylm_dm Array of psir_ylm_dm values.
 * \param stress_dot Output array for the dot product of stress components.
 * \param elements_num Number of elements in the input arrays.
 */

__global__ void dot_product_stress(double* psir_lxx,
                                   double* psir_lxy,
                                   double* psir_lxz,
                                   double* psir_lyy,
                                   double* psir_lyz,
                                   double* psir_lzz,
                                   double* psir_ylm_dm,
                                   double* stress_dot,
                                   int elements_num)
{

    __shared__ double cache[256][6]; // == threadsPerBlock
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int cacheIndex = threadIdx.x;
    double tmp[6] = {0.0};
    while (tid < elements_num)
    {
        tmp[0] += psir_lxx[tid] * psir_ylm_dm[tid] * 2;
        tmp[1] += psir_lxy[tid] * psir_ylm_dm[tid] * 2;
        tmp[2] += psir_lxz[tid] * psir_ylm_dm[tid] * 2;
        tmp[3] += psir_lyy[tid] * psir_ylm_dm[tid] * 2;
        tmp[4] += psir_lyz[tid] * psir_ylm_dm[tid] * 2;
        tmp[5] += psir_lzz[tid] * psir_ylm_dm[tid] * 2;
        tid += blockDim.x * gridDim.x;
    }

    for (int i = 0; i < 6; i++)
    {
        cache[cacheIndex][i] = tmp[i];
    }
    __syncthreads();

    int i = blockDim.x / 2;
    while (i != 0)
    {
        if (cacheIndex < i)
        {
            for (int index = 0; index < 6; index++)
            {
                cache[cacheIndex][index] += cache[cacheIndex + i][index];
            }
        }
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0){
        for (int index = 0; index < 6; index++)
        {
            stress_dot[blockIdx.x + gridDim.x * index] = cache[0][index];
        }
    }
}

/**
 * @brief Calculate the dot product force.
 *
 * This function calculates the dot product force based on the provided
 * parameters.
 *
 * @param psir_lx Pointer to the array of psir_lx values.
 * @param psir_ly Pointer to the array of psir_ly values.
 * @param psir_lz Pointer to the array of psir_lz values.
 * @param psir_ylm_dm Pointer to the array of psir_ylm_dm values.
 * @param force_dot Pointer to the array where the calculated force will be
 * stored.
 * @param iat Pointer to the array of iat values.
 * @param nwmax Maximum value for nwmax.
 * @param max_size Maximum size for arrays.
 * @param elements_num Number of elements to process.
 */

__global__ void dot_product_force(double* psir_lx,
                                  double* psir_ly,
                                  double* psir_lz,
                                  double* psir_ylm_dm,
                                  double* force_dot,
                                  int* iat,
                                  int nwmax,
                                  int max_size,
                                  int elements_num)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < elements_num)
    {
        int iat_on_nbz = iat[tid];
        if (iat_on_nbz <= -1)
        {
            tid += blockDim.x * gridDim.x;
            continue;
        }

        int iat_index = tid * 3;
        int dist = tid * nwmax;
        double tmp[3] = {0.0};

        for (int i = 0; i < nwmax; i++)
        {
            tmp[0] += psir_lx[dist + i] * psir_ylm_dm[dist + i] * 2;
            tmp[1] += psir_ly[dist + i] * psir_ylm_dm[dist + i] * 2;
            tmp[2] += psir_lz[dist + i] * psir_ylm_dm[dist + i] * 2;
        }

        for (int i = 0; i < 3; i++)
        {
            force_dot[iat_index + i] = tmp[i];
        }
        tid += blockDim.x * gridDim.x;
    }
}
void calculateInit(DensityMat& denstiy_mat,
                   ForceStressIatGlobal& f_s_iat_dev,
                   hamilt::HContainer<double>* dm,
                   const Grid_Technique& gridt,
                   const UnitCell& ucell,
                   int lgd,
                   int cuda_block,
                   int atom_num_grid)
{
    denstiy_mat.density_mat_h = new double[lgd * lgd];
    allocateDm(denstiy_mat.density_mat_h, dm, gridt, ucell);

    checkCuda(cudaMalloc((void**)&denstiy_mat.density_mat_d,
                         lgd * lgd * sizeof(double)));
    checkCuda(cudaMemcpy(denstiy_mat.density_mat_d,
                         denstiy_mat.density_mat_h,
                         lgd * lgd * sizeof(double),
                         cudaMemcpyHostToDevice));

    checkCuda(cudaMalloc((void**)&f_s_iat_dev.stress_global,
                         6 * cuda_block * gridt.nstreams * sizeof(double)));
    checkCuda(cudaMemset(f_s_iat_dev.stress_global,
                         0,
                         6 * cuda_block * gridt.nstreams * sizeof(double)));

    checkCuda(cudaMalloc((void**)&f_s_iat_dev.force_global,
                         3 * atom_num_grid * gridt.nstreams * sizeof(double)));
    checkCuda(cudaMemset(f_s_iat_dev.force_global,
                         0,
                         3 * atom_num_grid * gridt.nstreams * sizeof(double)));

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
void para_init(SGridParameter& para,
                       int iter_num,
                       int nbz,
                       const Grid_Technique& gridt)
{

    // stream_num stand for nstreams
    para.stream_num = iter_num % gridt.nstreams;
    //input_dou and input _int used for the Spherical Harmonics
    para.input_dou
        = &gridt.psi_dbl_gbl[gridt.psi_size_max * para.stream_num * 5];
    para.input_int
        = &gridt.psi_int_gbl[gridt.psi_size_max * para.stream_num * 2];
    para.num_psir = &gridt.num_psir_gbl[nbz * para.stream_num];
    //one dimension,record the length and the leading dimension of three matrix
    para.atom_pair_A_m
        = &gridt.l_info_global[gridt.atom_pair_nbz * para.stream_num];
    para.atom_pair_B_n
        = &gridt.r_info_global[gridt.atom_pair_nbz * para.stream_num];
    para.atom_pair_K
        = &gridt.k_info_global[gridt.atom_pair_nbz * para.stream_num];
    para.atom_pair_lda
        = &gridt.lda_info_global[gridt.atom_pair_nbz * para.stream_num];
    para.atom_pair_ldb
        = &gridt.ldb_info_global[gridt.atom_pair_nbz * para.stream_num];
    para.atom_pair_ldc
        = &gridt.ldc_info_global[gridt.atom_pair_nbz * para.stream_num];
    //input_double_g and input_int_g used for the Spherical Harmonics on GPU
    para.input_double_g
        = &gridt.psi_dbl_gbl_g[gridt.psi_size_max * para.stream_num * 5];
    para.input_int_g
        = &gridt.psi_int_gbl_g[gridt.psi_size_max * para.stream_num * 2];
    para.num_psir_g = &gridt.num_psir_gbl_g[nbz * para.stream_num];
    para.psir_dm_device = &gridt.dm_global_g[gridt.psir_size * para.stream_num];
    para.psir_r_device
        = &gridt.right_global_g[gridt.psir_size * para.stream_num];
    //psi function ,record the force in x y z,and the stress in six dimension
    para.psir_lx_device = &gridt.d_left_x_g[gridt.psir_size * para.stream_num];
    para.psir_ly_device = &gridt.d_left_y_g[gridt.psir_size * para.stream_num];
    para.psir_lz_device = &gridt.d_left_z_g[gridt.psir_size * para.stream_num];
    para.psir_lxx_device
        = &gridt.dd_left_xx_g[gridt.psir_size * para.stream_num];
    para.psir_lxy_device
        = &gridt.dd_left_xy_g[gridt.psir_size * para.stream_num];
    para.psir_lxz_device
        = &gridt.dd_left_xz_g[gridt.psir_size * para.stream_num];
    para.psir_lyy_device
        = &gridt.dd_left_yy_g[gridt.psir_size * para.stream_num];
    para.psir_lyz_device
        = &gridt.dd_left_yz_g[gridt.psir_size * para.stream_num];
    para.psir_lzz_device
        = &gridt.dd_left_zz_g[gridt.psir_size * para.stream_num];
    //one dimension,record the length and the leading dimension of three matrix on GPU
    para.A_m_device
        = &gridt.l_info_global_g[gridt.atom_pair_nbz * para.stream_num];
    para.B_n_device
        = &gridt.r_info_global_g[gridt.atom_pair_nbz * para.stream_num];
    para.K_device
        = &gridt.k_info_global_g[gridt.atom_pair_nbz * para.stream_num];
    para.lda_device
        = &gridt.lda_info_gbl_g[gridt.atom_pair_nbz * para.stream_num];
    para.ldb_device
        = &gridt.ldb_info_gbl_g[gridt.atom_pair_nbz * para.stream_num];
    para.ldc_device
        = &gridt.ldc_info_gbl_g[gridt.atom_pair_nbz * para.stream_num];
    //two dimension,record number to compute
    para.matrix_A = &gridt.ap_left_gbl[gridt.atom_pair_nbz * para.stream_num];
    para.matrix_B = &gridt.ap_right_gbl[gridt.atom_pair_nbz * para.stream_num];
    para.matrix_C = &gridt.ap_output_gbl[gridt.atom_pair_nbz * para.stream_num];
    para.matrix_A_device
        = &gridt.ap_left_gbl_g[gridt.atom_pair_nbz * para.stream_num];
    para.matrix_B_device
        = &gridt.ap_right_gbl_g[gridt.atom_pair_nbz * para.stream_num];
    para.matrix_C_device
        = &gridt.ap_output_gbl_g[gridt.atom_pair_nbz * para.stream_num];
}
/**
 * @brief ForceStressIat on host and device Init
 *
 * GridParameter init
 *
 * @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Host
 * @param stream_num int , record the stream in GPU
 * @param cuda_block in stress compute,used for Block nums
 * @param atom_num_grid in force calculate,used for Block nums
 * @param max_size Maximum size of atoms on a grid.
 * @param ForceStressIatGlobal ForceStressIatGlobal,contains the Force Stree Iat on Host
 */
void cal_init(ForceStressIat& f_s_iat,
                        const int stream_num,
                        const int cuda_block,
                        const int atom_num_grid,
                        const int max_size,
                        const ForceStressIatGlobal& f_s_iat_dev)
{
    const int iat_min = -max_size - 1;
    f_s_iat.stress_host = new double[6 * cuda_block];
    f_s_iat.stress_device
        = &f_s_iat_dev.stress_global[6 * cuda_block * stream_num];
    f_s_iat.force_device
        = &f_s_iat_dev.force_global[3 * atom_num_grid * stream_num];
    f_s_iat.iat_device
        = &f_s_iat_dev.iat_global[atom_num_grid * stream_num];
    f_s_iat.iat_host = new int[atom_num_grid];
    for (int index = 0; index < atom_num_grid; index++)
    {
        f_s_iat.iat_host[index] = iat_min;
    }
    f_s_iat.force_host = new double[3 * atom_num_grid];
    ModuleBase::GlobalFunc::ZEROS(f_s_iat.force_host,
                                  3 * atom_num_grid);
}

/**
 * @brief GridParameter memCpy,from Host to Device
 *
 * parameter init,which contains the gpu task and multi matrix multiplication
 *
 * @param para Grid parameter in task generator,
 * @param gridt Grid_Technique,stored the major method in the the gint.
 * @param nbz int,stand for the number of Z-axis
 * @param atom_num_grid in force calculate,used for Block nums
 */
void para_mem_copy(SGridParameter& para,
                         const Grid_Technique& gridt,
                         const int nbz,
                         const int atom_num_grid)
{
    checkCuda(cudaMemcpyAsync(para.input_double_g,
                              para.input_dou,
                              gridt.psi_size_max * 5 * sizeof(double),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.input_int_g,
                              para.input_int,
                              gridt.psi_size_max * 2 * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.num_psir_g,
                              para.num_psir,
                              nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.A_m_device,
                              para.atom_pair_A_m,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.B_n_device,
                              para.atom_pair_B_n,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.K_device,
                              para.atom_pair_K,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.lda_device,
                              para.atom_pair_lda,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.ldb_device,
                              para.atom_pair_ldb,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.ldc_device,
                              para.atom_pair_ldc,
                              gridt.atom_pair_nbz * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.matrix_A_device,
                              para.matrix_A,
                              gridt.atom_pair_nbz * sizeof(double*),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.matrix_B_device,
                              para.matrix_B,
                              gridt.atom_pair_nbz * sizeof(double*),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemcpyAsync(para.matrix_C_device,
                              para.matrix_C,
                              gridt.atom_pair_nbz * sizeof(double*),
                              cudaMemcpyHostToDevice,
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_dm_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_r_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_lx_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_ly_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_lz_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_lxx_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_lxy_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_lxz_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_lyy_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_lyz_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
    checkCuda(cudaMemsetAsync(para.psir_lzz_device,
                              0,
                              gridt.psir_size * sizeof(double),
                              gridt.streams[para.stream_num]));
}
/**
 * @brief Force Stress Force Iat memCpy,from Host to Device
 *
 *  @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 *  @param gridt Grid_Technique,stored the major method in the the gint.
 *  @param atom_num_grid in force calculate,used for Block nums
 *  @param cuda_block in stress compute,used for Block nums
 *  @param stream_num int , record the stream in GPU
 */
void cal_mem_cpy(ForceStressIat& f_s_iat,
                          const Grid_Technique& gridt,
                          const int atom_num_grid,
                          const int cuda_block,
                          const int stream_num)
{
    checkCuda(cudaMemcpyAsync(f_s_iat.iat_device,
                              f_s_iat.iat_host,
                              atom_num_grid * sizeof(int),
                              cudaMemcpyHostToDevice,
                              gridt.streams[stream_num]));
    checkCuda(cudaMemsetAsync(f_s_iat.stress_device,
                              0,
                              6 * cuda_block * sizeof(double),
                              gridt.streams[stream_num]));
    checkCuda(cudaMemsetAsync(f_s_iat.force_device,
                              0,
                              3 * atom_num_grid * sizeof(double),
                              gridt.streams[stream_num]));
}
/*
 * @brief Force Calculate on Host
 *
 * @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 * @param force stored the force for each atom on each directions
 * @param atom_num_grid in force calculate,used for Block nums
 */
void cal_force_add(ForceStressIat& f_s_iat,
                     double* force,
                    const int atom_num_grid)
{
    checkCuda(cudaMemcpy(f_s_iat.force_host,
                         f_s_iat.force_device,
                         3 * atom_num_grid * sizeof(double),
                         cudaMemcpyDeviceToHost));
    for (int index1 = 0; index1 < atom_num_grid; index1++)
    {
        int iat1 = f_s_iat.iat_host[index1];
        if (iat1 >= 0)
        {
            for (int index2 = 0; index2 < 3; index2++)
            {
                force[iat1 * 3 + index2]
                    += f_s_iat.force_host[index1 * 3 + index2];
            }
        }
    }
}
/**
 * @brief Stress Calculate on Host
 *
 * @param ForceStressIat ForceStressIat,contains the Force Stree Iat on Device
 * and Host
 * @param stress stored the stress for each directions
 * @param cuda_block in stress compute,used for Block nums
 */
void cal_stress_add(ForceStressIat& f_s_iat,
                     double* stress,
                     const int cuda_block)
{
    checkCuda(cudaMemcpy(f_s_iat.stress_host,
                         f_s_iat.stress_device,
                         6 * cuda_block * sizeof(double),
                         cudaMemcpyDeviceToHost));
    for (int i = 0; i < 6; i++)
    {
        for (int index = 0; index < cuda_block; index++)
        {
            // printf("the stress is %f\n",stress[i]);
            stress[i] += f_s_iat.stress_host[i * cuda_block + index];
        }
    }
}
} // namespace GintKernel
