#include "omp.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_vl.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/gint_vl.cuh"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

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

    int *atom_nw_now = new int[ntype];
    int *ucell_atom_nwl_now = new int[ntype];
    for (int i = 0; i < ntype; i++)
    {
        atom_nw_now[i] = GlobalC::ucell.atoms[i].nw;
        ucell_atom_nwl_now[i] = GlobalC::ucell.atoms[i].nwl;
    }

    int nr_max = static_cast<int>(1000 * max_cut) + 10;
    double *psi_u_now = new double[ntype * nwmax * nr_max * 2];
    memset(psi_u_now, 0, ntype * nwmax * nr_max * 2 * sizeof(double));
    bool *atom_iw2_new_now = new bool[ntype * nwmax];
    memset(atom_iw2_new_now, 0, ntype * nwmax * sizeof(bool));
    int *atom_iw2_ylm_now = new int[ntype * nwmax];
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

    int *ucell_atom_nwl;
    checkCuda(cudaMalloc((void **)&ucell_atom_nwl, ntype * sizeof(int)));
    checkCuda(cudaMemcpy(ucell_atom_nwl, ucell_atom_nwl_now, ntype * sizeof(int), cudaMemcpyHostToDevice));

    double *psi_u;
    checkCuda(cudaMalloc((void **)&psi_u, ntype * nwmax * nr_max * sizeof(double) * 2));
    checkCuda(cudaMemcpy(psi_u, psi_u_now, ntype * nwmax * nr_max * sizeof(double) * 2, cudaMemcpyHostToDevice));

    bool *atom_iw2_new;
    int *atom_iw2_ylm;
    checkCuda(cudaMalloc((void **)&atom_iw2_new, ntype * nwmax * sizeof(bool)));
    checkCuda(cudaMalloc((void **)&atom_iw2_ylm, ntype * nwmax * sizeof(int)));
    checkCuda(cudaMemcpy(atom_iw2_new, atom_iw2_new_now, ntype * nwmax * sizeof(bool), cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpy(atom_iw2_ylm, atom_iw2_ylm_now, ntype * nwmax * sizeof(int), cudaMemcpyHostToDevice));

    const int max_atom_pair_number = GlobalC::ucell.nat * GlobalC::ucell.nat;
    double *GridVlocal_v2_g[max_atom_pair_number];
    for (int i = 0; i < max_atom_pair_number; i++)
    {
        GridVlocal_v2_g[i] = nullptr;
    }

    const int nStreams = 4;
    const int psir_size = nbz * max_size * GridT.bxyz * nwmax;
    double *psir_ylm_left_global;
    double *psir_ylm_right_global;
    checkCuda(cudaMalloc((void **)&psir_ylm_left_global, psir_size * nStreams * sizeof(double)));
    checkCuda(cudaMalloc((void **)&psir_ylm_right_global, psir_size * nStreams * sizeof(double)));
    checkCuda(cudaMemset(psir_ylm_left_global, 0, psir_size * nStreams * sizeof(double)));
    checkCuda(cudaMemset(psir_ylm_right_global, 0, psir_size * nStreams * sizeof(double)));

    const int atom_pair_size_of_meshcell_v2 = max_size * max_size;
    const int atom_pair_size_over_nbz_v2 = atom_pair_size_of_meshcell_v2 * nbz;
    const int atom_pair_size_over_nbz = atom_pair_size_of_meshcell_v2 * 2 * nbz;

    int *atom_pair_input_info_global;
    checkCuda(cudaMallocHost((void **)&atom_pair_input_info_global, atom_pair_size_over_nbz * nStreams * sizeof(int)));
    int *atom_pair_input_info_g_global;
    checkCuda(cudaMalloc((void **)&atom_pair_input_info_g_global, atom_pair_size_over_nbz * nStreams * sizeof(int)));

    double ** atom_pair_left_global;
    double ** atom_pair_right_global;
    double ** atom_pair_output_global;
    checkCuda(cudaMallocHost((void **)&atom_pair_left_global, atom_pair_size_over_nbz_v2 * nStreams * sizeof(double *)));
    checkCuda(cudaMallocHost((void **)&atom_pair_right_global, atom_pair_size_over_nbz_v2 * nStreams * sizeof(double *)));
    checkCuda(cudaMallocHost((void **)&atom_pair_output_global, atom_pair_size_over_nbz_v2 * nStreams * sizeof(double *)));

    double ** atom_pair_left_global_g;
    double ** atom_pair_right_global_g;
    double ** atom_pair_output_global_g;
    checkCuda(cudaMalloc((void **)&atom_pair_left_global_g, atom_pair_size_over_nbz_v2 * nStreams * sizeof(double *)));
    checkCuda(cudaMalloc((void **)&atom_pair_right_global_g, atom_pair_size_over_nbz_v2 * nStreams * sizeof(double *)));
    checkCuda(cudaMalloc((void **)&atom_pair_output_global_g, atom_pair_size_over_nbz_v2 * nStreams * sizeof(double *)));

    int *num_atom_pair_global;
    checkCuda(cudaMallocHost((void **)&num_atom_pair_global, nbz * nStreams * sizeof(int)));
    int *num_atom_pair_g_global;
    checkCuda(cudaMalloc((void **)&num_atom_pair_g_global, nbz * nStreams * sizeof(int)));

    const int psi_size_max = max_size * GridT.bxyz * nbz;
    const int psi_size_max_per_z = max_size * GridT.bxyz;
    double *psi_input_double_global;
    checkCuda(cudaMallocHost((void **)&psi_input_double_global, psi_size_max * nStreams * 5 * sizeof(double)));
    double *psi_input_double_g_global;
    checkCuda(cudaMalloc((void **)&psi_input_double_g_global, psi_size_max * nStreams * 5 * sizeof(double)));

    int *psi_input_int_global;
    checkCuda(cudaMallocHost((void **)&psi_input_int_global, psi_size_max * nStreams * 2 * sizeof(int)));
    int *psi_input_int_g_global;
    checkCuda(cudaMalloc((void **)&psi_input_int_g_global, psi_size_max * nStreams * 2 * sizeof(int)));

    int *num_psir_global;
    checkCuda(cudaMallocHost((void **)&num_psir_global, nbz * nStreams * sizeof(int)));
    int *num_psir_g_global;
    checkCuda(cudaMalloc((void **)&num_psir_g_global, nbz * nStreams * sizeof(int)));

    int devId = 0;
    cudaDeviceProp prop;
    checkCuda(cudaGetDeviceProperties(&prop, devId));
    // printf("Device : %s\n", prop.name);
    checkCuda(cudaSetDevice(devId));

    cudaStream_t stream[nStreams];
    for (int i = 0; i < nStreams; ++i)
    {
        checkCuda(cudaStreamCreate(&stream[i]));
    }

    int iter_num = 0;
    //int omp_thread_num = omp_get_num_threads();
    //omp_set_num_threads(nStreams);

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
            int *atom_pair_input_info = &atom_pair_input_info_global[atom_pair_size_over_nbz * stream_num];
            int *num_atom_pair = &num_atom_pair_global[nbz * stream_num];

            double *psi_input_double_g = &psi_input_double_g_global[psi_size_max * stream_num * 5];
            int *psi_input_int_g = &psi_input_int_g_global[psi_size_max * stream_num * 2];
            int *num_psir_g = &num_psir_g_global[nbz * stream_num];
            double *psir_ylm_left_g = &psir_ylm_left_global[psir_size * stream_num];
            double *psir_ylm_right_g = &psir_ylm_right_global[psir_size * stream_num];
            int *atom_pair_input_info_g = &atom_pair_input_info_g_global[atom_pair_size_over_nbz * stream_num];
            int *num_atom_pair_g = &num_atom_pair_g_global[nbz * stream_num];


            double ** atom_pair_left = &atom_pair_left_global[atom_pair_size_over_nbz_v2 * stream_num];
            double ** atom_pair_right = &atom_pair_right_global[atom_pair_size_over_nbz_v2 * stream_num];
            double ** atom_pair_output = &atom_pair_output_global[atom_pair_size_over_nbz_v2 * stream_num];

            double ** atom_pair_left_g = &atom_pair_left_global_g[atom_pair_size_over_nbz_v2 * stream_num];
            double ** atom_pair_right_g = &atom_pair_right_global_g[atom_pair_size_over_nbz_v2 * stream_num];
            double ** atom_pair_output_g = &atom_pair_output_global_g[atom_pair_size_over_nbz_v2 * stream_num];

            gpu_task_generate_vlocal(GridT, i, j,
                                     atom_pair_size_of_meshcell_v2,
                                     psi_size_max_per_z,
                                     max_size, nczp,
                                     vfactor,
                                     vlocal,
                                     psir_ylm_left_g,
                                     psir_ylm_right_g,
                                     psi_input_double,
                                     psi_input_int,
                                     num_psir,
                                     atom_pair_input_info,
                                     num_atom_pair, GridVlocal_v2_g,
                                     atom_pair_left,
                                     atom_pair_right,
                                     atom_pair_output);

            checkCuda(cudaStreamSynchronize(stream[stream_num]));

            checkCuda(cudaMemcpyAsync(psi_input_double_g, psi_input_double, psi_size_max * 5 * sizeof(double), cudaMemcpyHostToDevice, stream[stream_num]));
            checkCuda(cudaMemcpyAsync(psi_input_int_g, psi_input_int, psi_size_max * 2 * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));
            checkCuda(cudaMemcpyAsync(num_psir_g, num_psir, nbz * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_input_info_g, atom_pair_input_info, atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));
            checkCuda(cudaMemcpyAsync(num_atom_pair_g, num_atom_pair, nbz * sizeof(int), cudaMemcpyHostToDevice, stream[stream_num]));

            checkCuda(cudaMemcpyAsync(atom_pair_left_g, atom_pair_left, atom_pair_size_over_nbz_v2 * sizeof(double*), cudaMemcpyHostToDevice, stream[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_right_g, atom_pair_right, atom_pair_size_over_nbz_v2 * sizeof(double*), cudaMemcpyHostToDevice, stream[stream_num]));
            checkCuda(cudaMemcpyAsync(atom_pair_output_g, atom_pair_output, atom_pair_size_over_nbz_v2 * sizeof(double*), cudaMemcpyHostToDevice, stream[stream_num]));

            checkCuda(cudaMemsetAsync(psir_ylm_left_g, 0, psir_size * sizeof(double), stream[stream_num]));
            checkCuda(cudaMemsetAsync(psir_ylm_right_g, 0, psir_size * sizeof(double), stream[stream_num]));

            dim3 grid_psi(nbz, 8);
            dim3 block_psi(64);
            dim3 grid_multiple(nbz, 512);
            dim3 block_multiple(128);

            get_psi_and_vldr3<<<grid_psi, block_psi, 0, stream[stream_num]>>>(psi_input_double_g,
                                                                              psi_input_int_g,
                                                                              num_psir_g,
                                                                              psi_size_max_per_z,
                                                                              ucell_atom_nwl,
                                                                              atom_iw2_new,
                                                                              atom_iw2_ylm,
                                                                              atom_nw_g,
                                                                              nr_max,
                                                                              psi_u,
                                                                              psir_ylm_left_g,
                                                                              psir_ylm_right_g);
            psi_multiple<<<grid_multiple, block_multiple, 0, stream[stream_num]>>>(atom_pair_left_g,
                                                                                   atom_pair_right_g,
                                                                                   atom_pair_output_g,
                                                                                   atom_pair_input_info_g,
                                                                                   num_atom_pair_g,
                                                                                   atom_pair_size_of_meshcell_v2,
                                                                                   lgd);
            iter_num++;
        }
    }

    for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
    {
        for (int iat2 = 0; iat2 < GlobalC::ucell.nat; iat2++)
        {
            hamilt::AtomPair<double> *tmp_ap = hRGint->find_pair(iat1, iat2);
            if (tmp_ap == nullptr)
                continue;
            checkCuda(cudaMemcpy(tmp_ap->get_pointer(0), 
                                    GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2], 
                                    tmp_ap->get_row_size() * tmp_ap->get_col_size() * sizeof(double), 
                                    cudaMemcpyDeviceToHost));
        }
    }


    // free
    for (int i = 0; i < nStreams; ++i)
        checkCuda(cudaStreamDestroy(stream[i]));


    checkCuda(cudaFree(ucell_atom_nwl));
    checkCuda(cudaFree(psi_u));
    checkCuda(cudaFree(atom_iw2_new));
    checkCuda(cudaFree(atom_iw2_ylm));

    checkCuda(cudaFree(atom_pair_input_info_g_global));
    checkCuda(cudaFree(num_atom_pair_g_global));

    checkCuda(cudaFree(psi_input_double_g_global));
    checkCuda(cudaFree(psi_input_int_g_global));
    checkCuda(cudaFree(num_psir_g_global));

    checkCuda(cudaFreeHost(atom_pair_input_info_global));
    checkCuda(cudaFreeHost(num_atom_pair_global));
    checkCuda(cudaFreeHost(psi_input_double_global));
    checkCuda(cudaFreeHost(psi_input_int_global));
    checkCuda(cudaFreeHost(num_psir_global));

    checkCuda(cudaFreeHost(atom_pair_left_global));
    checkCuda(cudaFreeHost(atom_pair_right_global));
    checkCuda(cudaFreeHost(atom_pair_output_global));

    checkCuda(cudaFree(atom_pair_left_global_g));
    checkCuda(cudaFree(atom_pair_right_global_g));
    checkCuda(cudaFree(atom_pair_output_global_g));

    checkCuda(cudaFree(psir_ylm_left_global));
    checkCuda(cudaFree(psir_ylm_right_global));
    checkCuda(cudaFree(atom_nw_g));


    for (int i = 0; i < max_atom_pair_number; i++)
    {
        if (GridVlocal_v2_g[i] != nullptr)
        {
            checkCuda(cudaFree(GridVlocal_v2_g[i]));
        }
    }

    delete[] atom_nw_now;
    delete[] ucell_atom_nwl_now;
    delete[] psi_u_now;
    delete[] atom_iw2_new_now;
    delete[] atom_iw2_ylm_now;
}