#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_vl.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_vl.cuh"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "vbatch_matrix_multiple/vbatch_matrix_mul.cuh"
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(__DEBUG)
    if (result != cudaSuccess)
    {
        fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
#endif
    return result;
}

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

            gpu_task_generate_vlocal(GridT, i, j,
                        GridT.atom_pair_size_of_meshcell,
                        GridT.psi_size_max_per_z,
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
                        GridT.GridVlocal_v2_g,
                        atom_pair_mat_A_array,
                        atom_pair_mat_B_array,
                        atom_pair_mat_C_array,
                        atom_pair_num);

            checkCuda(cudaStreamSynchronize(GridT.streams[stream_num]));

            checkCuda(cudaMemcpyAsync(psi_input_double_g, psi_input_double, GridT.psi_size_max * 5 * sizeof(double), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(psi_input_int_g, psi_input_int, GridT.psi_size_max * 2 * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(num_psir_g, num_psir, nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));


            checkCuda(cudaMemcpyAsync(atom_pair_A_m_g, atom_pair_A_m, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_B_n_g, atom_pair_B_n, GridT.atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, GridT.streams[stream_num]));
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

            vbatch_gemm<double>(GlobalC::ucell.nwmax, GlobalC::ucell.nwmax,
                                atom_pair_A_m_g, atom_pair_B_n_g, GridT.bxyz,
                                atom_pair_mat_A_array_g, atom_pair_lda_g,
                                atom_pair_mat_B_array_g, atom_pair_ldb_g,
                                atom_pair_mat_C_array_g, atom_pair_ldc_g,
                                atom_pair_num, GridT.streams[stream_num]);
            /*
            dim3 grid_multiple(atom_pair_num);
            dim3 block_multiple(128);
            psi_multiple<<<grid_multiple, block_multiple, 0, GridT.streams[stream_num]>>>(atom_pair_A_m_g, atom_pair_B_n_g,
                                                                                    atom_pair_mat_A_array_g,
                                                                                    atom_pair_mat_B_array_g,
                                                                                    atom_pair_mat_C_array_g);*/
            iter_num++;
        }
    }
    for (int i = 0; i < GridT.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(GridT.streams[i]));
    }

    for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
    {
        for (int iat2 = 0; iat2 < GlobalC::ucell.nat; iat2++)
        {
            int stream_num = iter_num % GridT.nstreams;
            hamilt::AtomPair<double> *tmp_ap = hRGint->find_pair(iat1, iat2);
            if (tmp_ap == nullptr)
                continue;

             checkCuda(cudaMemcpyAsync(tmp_ap->get_pointer(0),
                                    GridT.GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2],
                                    tmp_ap->get_row_size() * tmp_ap->get_col_size() * sizeof(double),
                                    cudaMemcpyDeviceToHost, GridT.streams[stream_num]));
            checkCuda(cudaMemsetAsync(GridT.GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2],
                                    0,
                                    tmp_ap->get_row_size() * tmp_ap->get_col_size() * sizeof(double), GridT.streams[stream_num]));
      /*
            checkCuda(cudaMemcpy(tmp_ap->get_pointer(0),
                                    GridT.GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2],
                                    tmp_ap->get_row_size() * tmp_ap->get_col_size() * sizeof(double),
                                    cudaMemcpyDeviceToHost));
            checkCuda(cudaMemset(GridT.GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2],
                                    0,
                                    tmp_ap->get_row_size() * tmp_ap->get_col_size() * sizeof(double)));*/

        }
    }
    for (int i = 0; i < GridT.nstreams; i++)
    {
        checkCuda(cudaStreamSynchronize(GridT.streams[i]));
    }
}

void test_vbatch_gemm(const int max_m, const int max_n, const int k,
                 const int batchCount)
{

    int* h_m = new int[batchCount];
    int* h_n = new int[batchCount];
    int* h_global_lda = new int[batchCount];
    int* h_global_ldb = new int[batchCount];
    int* h_global_ldc = new int[batchCount];

    double** h_global_A_array = new double*[batchCount];
    double** h_global_B_array = new double*[batchCount];
    double** h_global_C_array = new double*[batchCount];

    double* h_global_A = new double[batchCount * max_m * k];
    double* h_global_B = new double[batchCount * max_n * k ];
    double* h_global_C = new double[batchCount * max_m * max_n];

    for (int i = 0; i < batchCount * max_m * k; ++i) {
        h_global_A[i] = i * 0.1;
    }
    for (int i = 0; i < batchCount * max_n * k; ++i) {
        h_global_B[i] = i * 0.2;
    }
    for (int i = 0; i < batchCount * max_m * max_n; ++i) {
        h_global_C[i] = 0.0;
    }
    // Allocate device memory
    int* d_m;
    int* d_n;
    int* d_global_lda;
    int* d_global_ldb;
    int* d_global_ldc;

    double** d_global_A_array;
    double** d_global_B_array;
    double** d_global_C_array;

    double* d_global_A;
    double* d_global_B;
    double* d_global_C;

    cudaMalloc(&d_m, batchCount * sizeof(int));
    cudaMalloc(&d_n, batchCount * sizeof(int));
    cudaMalloc(&d_global_lda, batchCount * sizeof(int));
    cudaMalloc(&d_global_ldb, batchCount * sizeof(int));
    cudaMalloc(&d_global_ldc, batchCount * sizeof(int));
    cudaMalloc(&d_global_A_array, batchCount * sizeof(double*));
    cudaMalloc(&d_global_B_array, batchCount * sizeof(double*));
    cudaMalloc(&d_global_C_array, batchCount * sizeof(double*));

    cudaMalloc(&d_global_A, batchCount * max_m * k * sizeof(double));
    cudaMalloc(&d_global_B, batchCount * max_n * k * sizeof(double));
    cudaMalloc(&d_global_C, batchCount * max_m * max_n * sizeof(double));


    for (int i = 0; i < batchCount; ++i) {
        h_m[i] = max_m;
        h_n[i] = max_n;
        h_global_lda[i] = k;
        h_global_ldb[i] = k;
        h_global_ldc[i] = max_n;

        h_global_A_array[i] = &d_global_A[i * max_m * k];
        h_global_B_array[i] = &d_global_B[i * max_n * k];
        h_global_C_array[i] = &d_global_C[0];  // test atom add
    }

    cudaMemcpy(d_m, h_m, batchCount * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_n, h_n, batchCount * sizeof(int), cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_global_lda, h_global_lda, batchCount * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_ldb, h_global_ldb, batchCount * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_ldc, h_global_ldc, batchCount * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_global_A_array, h_global_A_array, batchCount * sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_B_array, h_global_B_array, batchCount * sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_C_array, h_global_C_array, batchCount * sizeof(double*), cudaMemcpyHostToDevice);

    cudaMemcpy(d_global_A, h_global_A, batchCount * max_m * k * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_B, h_global_B, batchCount * max_n * k * sizeof(double), cudaMemcpyHostToDevice);


    cudaStream_t stream;
    cudaStreamCreate(&stream);

    vbatch_gemm<double>(max_n, max_m,  
                        d_n, d_m,  k,
                        d_global_B_array, d_global_ldb,
                        d_global_A_array, d_global_lda,
                        d_global_C_array, d_global_ldc,
                        batchCount, stream);

    cudaStreamSynchronize(stream);
    cudaMemcpy(h_global_C, d_global_C, batchCount * max_m * max_n * sizeof(double), cudaMemcpyDeviceToHost);
    {
        int index = 0;
        for (size_t i = 0; i < batchCount; i++)
        {
            for (size_t j = 0; j < max_m * k; j++)
            {
                std::cout << h_global_A[index++] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    {
        int index = 0;
        for (size_t i = 0; i < batchCount; i++)
        {
            for (size_t j = 0; j < max_n * k; j++)
            {
                std::cout << h_global_B[index++] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }


    {
        int index = 0;
        for (size_t i = 0; i < batchCount; i++)
        {
            for (size_t j = 0; j < max_m * max_n; j++)
            {
                std::cout << h_global_C[index++] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    delete[] h_global_A_array;
    delete[] h_global_B_array;
    delete[] h_global_C_array;
    
    delete[] h_m;
    delete[] h_n;
    
    delete[] h_global_lda;
    delete[] h_global_ldb;
    delete[] h_global_ldc;

    delete[] h_global_A;
    delete[] h_global_B;
    delete[] h_global_C;

    // Cleanup
    cudaFree(d_global_A_array);
    cudaFree(d_global_B_array);
    cudaFree(d_global_C_array);

    cudaFree(d_m);
    cudaFree(d_n);

    cudaFree(d_global_lda);
    cudaFree(d_global_ldb);
    cudaFree(d_global_ldc);

    cudaFree(d_global_A_array);
    cudaFree(d_global_B_array);
    cudaFree(d_global_C_array);

    cudaFree(d_global_A);
    cudaFree(d_global_B);
    cudaFree(d_global_C);

    cudaStreamDestroy(stream);
}