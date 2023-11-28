#include "omp.h"
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


    // test_vbatch_gemm(2, 3, 2, 1);

    const Numerical_Orbital_Lm *pointer;
    const int nbx = GridT.nbx;
    const int nby = GridT.nby;
    const int nbz = GridT.nbzp;
    const int nwmax = GlobalC::ucell.nwmax;
    const int ntype = GlobalC::ucell.ntype;
    gint_gamma_vl_upload_const(max_size, ylmcoef_now, GridT.bxyz);
    double max_cut = 0;
    for (int i = 0; i < ntype; i++)
    {
        if (GlobalC::ORB.Phi[i].getRcut() > max_cut)
        {
            max_cut = GlobalC::ORB.Phi[i].getRcut();
        }
    }

    int atom_nw_now[ntype];
    int ucell_atom_nwl_now[ntype];
    for (int i = 0; i < ntype; i++)
    {
        atom_nw_now[i] = GlobalC::ucell.atoms[i].nw;
        ucell_atom_nwl_now[i] = GlobalC::ucell.atoms[i].nwl;
    }

    int nr_max = static_cast<int>(1000 * max_cut) + 10;
    double psi_u_now[ntype * nwmax * nr_max * 2];
    memset(psi_u_now, 0, ntype * nwmax * nr_max * 2 * sizeof(double));
    bool atom_iw2_new_now[ntype * nwmax];
    memset(atom_iw2_new_now, 0, ntype * nwmax * sizeof(bool));
    int atom_iw2_ylm_now[ntype * nwmax];
    memset(atom_iw2_ylm_now, 0, ntype * nwmax * sizeof(int));

    Atom *atomx;
    for (int i = 0; i < ntype; i++)
    {
        atomx = &GlobalC::ucell.atoms[i];
        for (int j = 0; j < nwmax; j++)
        {
            if (j < atomx->nw)
            {
                atom_iw2_new_now[i * nwmax + j] = atomx->iw2_new[j];
                atom_iw2_ylm_now[i * nwmax + j] = atomx->iw2_ylm[j];
                pointer = &GlobalC::ORB.Phi[i].PhiLN(atomx->iw2l[j], atomx->iw2n[j]);
                for (int k = 0; k < nr_max; k++)
                {
                    int index_temp = (i * nwmax * nr_max + j * nr_max + k) * 2;
                    if (k < pointer->nr_uniform)
                    {
                        psi_u_now[index_temp] = pointer->psi_uniform[k];
                        psi_u_now[index_temp + 1] = pointer->dpsi_uniform[k];
                    }
                }
            }
        }
    }

    int *atom_nw_g;
    checkCuda(cudaMalloc((void **)&atom_nw_g, ntype * sizeof(int)));
    checkCuda(cudaMemcpy(atom_nw_g, atom_nw_now, ntype * sizeof(int), cudaMemcpyHostToDevice));

    int *ucell_atom_nwl_g;
    checkCuda(cudaMalloc((void **)&ucell_atom_nwl_g, ntype * sizeof(int)));
    checkCuda(cudaMemcpy(ucell_atom_nwl_g, ucell_atom_nwl_now, ntype * sizeof(int), cudaMemcpyHostToDevice));

    double *psi_u_g;
    checkCuda(cudaMalloc((void **)&psi_u_g, ntype * nwmax * nr_max * sizeof(double) * 2));
    checkCuda(cudaMemcpy(psi_u_g, psi_u_now, ntype * nwmax * nr_max * sizeof(double) * 2, cudaMemcpyHostToDevice));

    bool *atom_iw2_new_g;
    int *atom_iw2_ylm_g;
    checkCuda(cudaMalloc((void **)&atom_iw2_new_g, ntype * nwmax * sizeof(bool)));
    checkCuda(cudaMalloc((void **)&atom_iw2_ylm_g, ntype * nwmax * sizeof(int)));
    checkCuda(cudaMemcpy(atom_iw2_new_g, atom_iw2_new_now, ntype * nwmax * sizeof(bool), cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpy(atom_iw2_ylm_g, atom_iw2_ylm_now, ntype * nwmax * sizeof(int), cudaMemcpyHostToDevice));

    const int max_atom_pair_number = GlobalC::ucell.nat * GlobalC::ucell.nat;
    double *GridVlocal_v2_g[max_atom_pair_number];
    for (int i = 0; i < max_atom_pair_number; i++)
    {
        GridVlocal_v2_g[i] = nullptr;
    }

    const int nStreams = 4;
    const int psir_size = nbz * max_size * GridT.bxyz * nwmax;
    double *psir_ylm_left_global_g;
    double *psir_ylm_right_global_g;
    checkCuda(cudaMalloc((void **)&psir_ylm_left_global_g, psir_size * nStreams * sizeof(double)));
    checkCuda(cudaMalloc((void **)&psir_ylm_right_global_g, psir_size * nStreams * sizeof(double)));
    checkCuda(cudaMemset(psir_ylm_left_global_g, 0, psir_size * nStreams * sizeof(double)));
    checkCuda(cudaMemset(psir_ylm_right_global_g, 0, psir_size * nStreams * sizeof(double)));

    const int atom_pair_size_of_meshcell = max_size * max_size;
    const int atom_pair_size_over_nbz = atom_pair_size_of_meshcell * nbz;

    int *atom_pair_left_info_global;
    checkCuda(cudaMallocHost((void **)&atom_pair_left_info_global, atom_pair_size_over_nbz * nStreams * sizeof(int)));
    int *atom_pair_left_info_global_g;
    checkCuda(cudaMalloc((void **)&atom_pair_left_info_global_g, atom_pair_size_over_nbz * nStreams * sizeof(int)));

    int *atom_pair_right_info_global;
    checkCuda(cudaMallocHost((void **)&atom_pair_right_info_global, atom_pair_size_over_nbz * nStreams * sizeof(int)));
    int *atom_pair_right_info_global_g;
    checkCuda(cudaMalloc((void **)&atom_pair_right_info_global_g, atom_pair_size_over_nbz * nStreams * sizeof(int)));

    int *atom_pair_lda_global;
    checkCuda(cudaMallocHost((void **)&atom_pair_lda_global, atom_pair_size_over_nbz * nStreams * sizeof(int)));
    int *atom_pair_lda_global_g;
    checkCuda(cudaMalloc((void **)&atom_pair_lda_global_g, atom_pair_size_over_nbz * nStreams * sizeof(int)));

    int *atom_pair_ldb_global;
    checkCuda(cudaMallocHost((void **)&atom_pair_ldb_global, atom_pair_size_over_nbz * nStreams * sizeof(int)));
    int *atom_pair_ldb_global_g;
    checkCuda(cudaMalloc((void **)&atom_pair_ldb_global_g, atom_pair_size_over_nbz * nStreams * sizeof(int)));

    int *atom_pair_ldc_global;
    checkCuda(cudaMallocHost((void **)&atom_pair_ldc_global, atom_pair_size_over_nbz * nStreams * sizeof(int)));
    int *atom_pair_ldc_global_g;
    checkCuda(cudaMalloc((void **)&atom_pair_ldc_global_g, atom_pair_size_over_nbz * nStreams * sizeof(int)));

    double **atom_pair_left_global;
    double **atom_pair_right_global;
    double **atom_pair_output_global;
    checkCuda(cudaMallocHost((void **)&atom_pair_left_global, atom_pair_size_over_nbz * nStreams * sizeof(double *)));
    checkCuda(cudaMallocHost((void **)&atom_pair_right_global, atom_pair_size_over_nbz * nStreams * sizeof(double *)));
    checkCuda(cudaMallocHost((void **)&atom_pair_output_global, atom_pair_size_over_nbz * nStreams * sizeof(double *)));

    double **atom_pair_left_global_g;
    double **atom_pair_right_global_g;
    double **atom_pair_output_global_g;
    checkCuda(cudaMalloc((void **)&atom_pair_left_global_g, atom_pair_size_over_nbz * nStreams * sizeof(double *)));
    checkCuda(cudaMalloc((void **)&atom_pair_right_global_g, atom_pair_size_over_nbz * nStreams * sizeof(double *)));
    checkCuda(cudaMalloc((void **)&atom_pair_output_global_g, atom_pair_size_over_nbz * nStreams * sizeof(double *)));

    const int psi_size_max = max_size * GridT.bxyz * nbz;
    const int psi_size_max_per_z = max_size * GridT.bxyz;
    double *psi_input_double_global;
    checkCuda(cudaMallocHost((void **)&psi_input_double_global, psi_size_max * nStreams * 5 * sizeof(double)));
    double *psi_input_double_global_g;
    checkCuda(cudaMalloc((void **)&psi_input_double_global_g, psi_size_max * nStreams * 5 * sizeof(double)));

    int *psi_input_int_global;
    checkCuda(cudaMallocHost((void **)&psi_input_int_global, psi_size_max * nStreams * 2 * sizeof(int)));
    int *psi_input_int_global_g;
    checkCuda(cudaMalloc((void **)&psi_input_int_global_g, psi_size_max * nStreams * 2 * sizeof(int)));

    int *num_psir_global;
    checkCuda(cudaMallocHost((void **)&num_psir_global, nbz * nStreams * sizeof(int)));
    int *num_psir_global_g;
    checkCuda(cudaMalloc((void **)&num_psir_global_g, nbz * nStreams * sizeof(int)));

    checkCuda(cudaSetDevice(0));
    cudaStream_t stream[nStreams];
    for (int i = 0; i < nStreams; ++i)
    {
        checkCuda(cudaStreamCreate(&stream[i]));
    }
    {
        int iter_num = 0;
#pragma omp parallel for
        for (int i = 0; i < nbx; i++)
        {
#pragma omp parallel for
            for (int j = 0; j < nby; j++)
            {
                int stream_num = iter_num % nStreams;
                double *psi_input_double = &psi_input_double_global[psi_size_max * stream_num * 5];
                int *psi_input_int = &psi_input_int_global[psi_size_max * stream_num * 2];
                int *num_psir = &num_psir_global[nbz * stream_num];
                int *atom_pair_A_m = &atom_pair_left_info_global[atom_pair_size_over_nbz * stream_num];
                int *atom_pair_B_n = &atom_pair_right_info_global[atom_pair_size_over_nbz * stream_num];
                int *atom_pair_lda = &atom_pair_lda_global[atom_pair_size_over_nbz * stream_num];
                int *atom_pair_ldb = &atom_pair_ldb_global[atom_pair_size_over_nbz * stream_num];
                int *atom_pair_ldc = &atom_pair_ldc_global[atom_pair_size_over_nbz * stream_num];

                double *psi_input_double_g = &psi_input_double_global_g[psi_size_max * stream_num * 5];
                int *psi_input_int_g = &psi_input_int_global_g[psi_size_max * stream_num * 2];
                int *num_psir_g = &num_psir_global_g[nbz * stream_num];
                double *psir_ylm_left_g = &psir_ylm_left_global_g[psir_size * stream_num];
                double *psir_ylm_right_g = &psir_ylm_right_global_g[psir_size * stream_num];

                int *atom_pair_A_m_g = &atom_pair_left_info_global_g[atom_pair_size_over_nbz * stream_num];
                int *atom_pair_B_n_g = &atom_pair_right_info_global_g[atom_pair_size_over_nbz * stream_num];
                int *atom_pair_lda_g = &atom_pair_lda_global_g[atom_pair_size_over_nbz * stream_num];
                int *atom_pair_ldb_g = &atom_pair_ldb_global_g[atom_pair_size_over_nbz * stream_num];
                int *atom_pair_ldc_g = &atom_pair_ldc_global_g[atom_pair_size_over_nbz * stream_num];

                double **atom_pair_mat_A_array = &atom_pair_left_global[atom_pair_size_over_nbz * stream_num];
                double **atom_pair_mat_B_array = &atom_pair_right_global[atom_pair_size_over_nbz * stream_num];
                double **atom_pair_mat_C_array = &atom_pair_output_global[atom_pair_size_over_nbz * stream_num];

                double **atom_pair_mat_A_array_g = &atom_pair_left_global_g[atom_pair_size_over_nbz * stream_num];
                double **atom_pair_mat_B_array_g = &atom_pair_right_global_g[atom_pair_size_over_nbz * stream_num];
                double **atom_pair_mat_C_array_g = &atom_pair_output_global_g[atom_pair_size_over_nbz * stream_num];
                int atom_pair_num = 0;

                gpu_task_generate_vlocal(GridT, i, j,
                            atom_pair_size_of_meshcell,
                            psi_size_max_per_z,
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
                            GridVlocal_v2_g,
                            atom_pair_mat_A_array,
                            atom_pair_mat_B_array,
                            atom_pair_mat_C_array,
                            atom_pair_num);

                checkCuda(cudaStreamSynchronize(stream[stream_num]));

                checkCuda(cudaMemcpyAsync(psi_input_double_g, psi_input_double, psi_size_max * 5 * sizeof(double), cudaMemcpyHostToDevice, stream[stream_num]));
                checkCuda(cudaMemcpyAsync(psi_input_int_g, psi_input_int, psi_size_max * 2 * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));
                checkCuda(cudaMemcpyAsync(num_psir_g, num_psir, nbz * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));


                checkCuda(cudaMemcpyAsync(atom_pair_A_m_g, atom_pair_A_m, atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));
                checkCuda(cudaMemcpyAsync(atom_pair_B_n_g, atom_pair_B_n, atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));
                checkCuda(cudaMemcpyAsync(atom_pair_lda_g, atom_pair_lda, atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));
                checkCuda(cudaMemcpyAsync(atom_pair_ldb_g, atom_pair_ldb, atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));
                checkCuda(cudaMemcpyAsync(atom_pair_ldc_g, atom_pair_ldc, atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));

                checkCuda(cudaMemcpyAsync(atom_pair_mat_A_array_g, atom_pair_mat_A_array, atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, stream[stream_num]));
                checkCuda(cudaMemcpyAsync(atom_pair_mat_B_array_g, atom_pair_mat_B_array, atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, stream[stream_num]));
                checkCuda(cudaMemcpyAsync(atom_pair_mat_C_array_g, atom_pair_mat_C_array, atom_pair_size_over_nbz * sizeof(double *), cudaMemcpyHostToDevice, stream[stream_num]));

                checkCuda(cudaMemsetAsync(psir_ylm_left_g, 0, psir_size * sizeof(double), stream[stream_num]));
                checkCuda(cudaMemsetAsync(psir_ylm_right_g, 0, psir_size * sizeof(double), stream[stream_num]));

                dim3 grid_psi(nbz, 8);
                dim3 block_psi(64);
                dim3 grid_multiple(atom_pair_num);
                dim3 block_multiple(128);

                get_psi_and_vldr3<<<grid_psi, block_psi, 0, stream[stream_num]>>>(psi_input_double_g,
                                                                                  psi_input_int_g,
                                                                                  num_psir_g,
                                                                                  psi_size_max_per_z,
                                                                                  ucell_atom_nwl_g,
                                                                                  atom_iw2_new_g,
                                                                                  atom_iw2_ylm_g,
                                                                                  atom_nw_g,
                                                                                  nr_max,
                                                                                  psi_u_g,
                                                                                  psir_ylm_left_g,
                                                                                  psir_ylm_right_g);

                vbatch_gemm<double>(nwmax, nwmax,
                                    atom_pair_A_m_g, atom_pair_B_n_g, GridT.bxyz,
                                    atom_pair_mat_A_array_g, atom_pair_lda_g,
                                    atom_pair_mat_B_array_g, atom_pair_ldb_g,
                                    atom_pair_mat_C_array_g, atom_pair_ldc_g,
                                    atom_pair_num, stream[stream_num]);
                /*psi_multiple<<<grid_multiple, block_multiple, 0, stream[stream_num]>>>(atom_pair_A_m_g, atom_pair_B_n_g,
                                                                                        atom_pair_mat_A_array_g,
                                                                                        atom_pair_mat_B_array_g,
                                                                                        atom_pair_mat_C_array_g);*/
                iter_num++;
            }
        }
    }
    {
        int iter_num = 0;
        for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
        {
            for (int iat2 = 0; iat2 < GlobalC::ucell.nat; iat2++)
            {
                hamilt::AtomPair<double> *tmp_ap = hRGint->find_pair(iat1, iat2);
                if (tmp_ap == nullptr)
                    continue;
                int stream_num = iter_num % nStreams;

                checkCuda(cudaMemcpyAsync(tmp_ap->get_pointer(0),
                                     GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2],
                                     tmp_ap->get_row_size() * tmp_ap->get_col_size() * sizeof(double),
                                     cudaMemcpyDeviceToHost, stream[stream_num]));
                iter_num++;
            }
        }
    }
    // free
    for (int i = 0; i < nStreams; ++i)
        checkCuda(cudaStreamDestroy(stream[i]));

    checkCuda(cudaFree(ucell_atom_nwl_g));
    checkCuda(cudaFree(psi_u_g));
    checkCuda(cudaFree(atom_iw2_new_g));
    checkCuda(cudaFree(atom_iw2_ylm_g));
    checkCuda(cudaFree(atom_nw_g));

    checkCuda(cudaFreeHost(psi_input_double_global));
    checkCuda(cudaFreeHost(psi_input_int_global));
    checkCuda(cudaFreeHost(num_psir_global));

    checkCuda(cudaFree(psi_input_double_global_g));
    checkCuda(cudaFree(psi_input_int_global_g));
    checkCuda(cudaFree(num_psir_global_g));
    checkCuda(cudaFree(psir_ylm_left_global_g));
    checkCuda(cudaFree(psir_ylm_right_global_g));

    checkCuda(cudaFreeHost(atom_pair_left_info_global));
    checkCuda(cudaFree(atom_pair_left_info_global_g));

    checkCuda(cudaFreeHost(atom_pair_right_info_global));
    checkCuda(cudaFree(atom_pair_right_info_global_g));

    checkCuda(cudaFreeHost(atom_pair_lda_global));
    checkCuda(cudaFree(atom_pair_lda_global_g));

    checkCuda(cudaFreeHost(atom_pair_ldb_global));
    checkCuda(cudaFree(atom_pair_ldb_global_g));

    checkCuda(cudaFreeHost(atom_pair_ldc_global));
    checkCuda(cudaFree(atom_pair_ldc_global_g));


    checkCuda(cudaFreeHost(atom_pair_left_global));
    checkCuda(cudaFreeHost(atom_pair_right_global));
    checkCuda(cudaFreeHost(atom_pair_output_global));

    checkCuda(cudaFree(atom_pair_left_global_g));
    checkCuda(cudaFree(atom_pair_right_global_g));
    checkCuda(cudaFree(atom_pair_output_global_g));

    for (int i = 0; i < max_atom_pair_number; i++)
    {
        if (GridVlocal_v2_g[i] != nullptr)
        {
            checkCuda(cudaFree(GridVlocal_v2_g[i]));
        }
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