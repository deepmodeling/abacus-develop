#include "omp.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "gint_vl.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include <fstream>
#include <sstream>
#define __DEBUG
__constant__ double ylmcoef[36];
__constant__ int bxyz_g[1];
__constant__ int max_size_g[1];
__constant__ int nwmax_g[1];
__constant__ double delta_r_g[1];

void dump_cuda_array_to_file(double *cuda_array, int width, int hight, const std::string &filename)
{
    double *h_data = new double[width * hight];
    cudaMemcpy(h_data, cuda_array, width * hight * sizeof(double), cudaMemcpyDeviceToHost);

    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
        std::cerr << "Failed to open file for writing." << std::endl;
    }
    for (int j = 0; j < hight; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            outFile << "hight" << j << "   width:" << i << "   " << h_data[j * width + i] << std::endl;
        }
    }
    outFile.close();
    delete[] h_data;
}

__global__ void get_psi_and_vldr3(double *input_double,
                        int *input_int,
                        int *num_psir,
                        int psi_size_max,
                        int *ucell_atom_nwl,
                        bool *atom_iw2_new,
                        int *atom_iw2_ylm,
                        int *atom_nw,
                        int nr_max,
                        double *psi_u,
                        double *psir_ylm_left,
                        double *psir_ylm_right)
{
    int size = num_psir[blockIdx.x];
    {
        int start_index = psi_size_max * blockIdx.x;
        int end_index = start_index + size;
        start_index += threadIdx.x;
        for (int index = start_index; index < end_index; index += blockDim.x)
        {

            double dr[3];
            int index_double = index * 5;
            dr[0] = input_double[index_double];
            dr[1] = input_double[index_double + 1];
            dr[2] = input_double[index_double + 2];
            double distance = input_double[index_double + 3];
            double vlbr3_value = input_double[index_double + 4];
            // begin calculation
            double ylma[49];    // Attention!!! At present, we only use L=5 at most. So (L+1) * (L+1)=36
            int index_int = index * 2;
            int it = input_int[index_int];
            int dist_tmp = input_int[index_int + 1];

            /***************************
            L = 0
            ***************************/
            ylma[0] = ylmcoef[0]; // l=0, m=0
            int nwl = ucell_atom_nwl[it];

            if (nwl == 0)
                goto YLM_END;

            /***************************
            L = 1
            ***************************/
            ylma[1] = ylmcoef[1] * dr[2];  // l=1, m=0
            ylma[2] = -ylmcoef[1] * dr[0]; // l=1, m=1
            ylma[3] = -ylmcoef[1] * dr[1]; // l=1, m=-1
            if (nwl == 1)
                goto YLM_END;

            /***************************
            L = 2
            ***************************/
            ylma[4] = ylmcoef[2] * dr[2] * ylma[1] - ylmcoef[3] * ylma[0]; // l=2, m=0
            {
                double tmp0 = ylmcoef[4] * dr[2];
                ylma[5] = tmp0 * ylma[2]; // l=2,m=1
                ylma[6] = tmp0 * ylma[3]; // l=2,m=-1

                tmp0 = ylmcoef[4] * dr[0];
                ylma[7] = ylmcoef[5] * ylma[4] - ylmcoef[6] * ylma[0] - tmp0 * ylma[2]; // l=2,m=2
                ylma[8] = -tmp0 * ylma[3];
                if (nwl == 2)
                    goto YLM_END;

                /***************************
                L = 3
                ***************************/
                ylma[9] = ylmcoef[7] * dr[2] * ylma[4] - ylmcoef[8] * ylma[1]; // l=3, m=0

                tmp0 = ylmcoef[9] * dr[2];
                ylma[10] = tmp0 * ylma[5] - ylmcoef[10] * ylma[2]; // l=3,m=1
                ylma[11] = tmp0 * ylma[6] - ylmcoef[10] * ylma[3]; // l=3,m=-1

                tmp0 = ylmcoef[11] * dr[2];
                ylma[12] = tmp0 * ylma[7]; // l=3,m=2
                ylma[13] = tmp0 * ylma[8]; // l=3,m=-2

                tmp0 = ylmcoef[14] * dr[0];
                ylma[14] = ylmcoef[12] * ylma[10] - ylmcoef[13] * ylma[2] - tmp0 * ylma[7]; // l=3,m=3
                ylma[15] = ylmcoef[12] * ylma[11] - ylmcoef[13] * ylma[3] - tmp0 * ylma[8]; // l=3,m=-3
                if (nwl == 3)
                    goto YLM_END;

                /***************************
                L = 4
                ***************************/
                ylma[16] = ylmcoef[15] * dr[2] * ylma[9] - ylmcoef[16] * ylma[4]; // l=4,m=0

                tmp0 = ylmcoef[17] * dr[2];
                ylma[17] = tmp0 * ylma[10] - ylmcoef[18] * ylma[5]; // l=4,m=1
                ylma[18] = tmp0 * ylma[11] - ylmcoef[18] * ylma[6]; // l=4,m=-1

                tmp0 = ylmcoef[19] * dr[2];
                ylma[19] = tmp0 * ylma[12] - ylmcoef[20] * ylma[7]; // l=4,m=2
                ylma[20] = tmp0 * ylma[13] - ylmcoef[20] * ylma[8]; // l=4,m=-2

                tmp0 = 3.0 * dr[2];
                ylma[21] = tmp0 * ylma[14]; // l=4,m=3
                ylma[22] = tmp0 * ylma[15]; // l=4,m=-3

                tmp0 = ylmcoef[23] * dr[0];
                ylma[23] = ylmcoef[21] * ylma[19] - ylmcoef[22] * ylma[7] - tmp0 * ylma[14]; // l=4,m=4
                ylma[24] = ylmcoef[21] * ylma[20] - ylmcoef[22] * ylma[8] - tmp0 * ylma[15]; // l=4,m=-4
                if (nwl == 4)
                    goto YLM_END;

                /***************************
                L = 5
                ***************************/
                ylma[25] = ylmcoef[24] * dr[2] * ylma[16] - ylmcoef[25] * ylma[9]; // l=5,m=0

                tmp0 = ylmcoef[26] * dr[2];
                ylma[26] = tmp0 * ylma[17] - ylmcoef[27] * ylma[10]; // l=5,m=1
                ylma[27] = tmp0 * ylma[18] - ylmcoef[27] * ylma[11]; // l=5,m=-1

                tmp0 = ylmcoef[28] * dr[2];
                ylma[28] = tmp0 * ylma[19] - ylmcoef[29] * ylma[12]; // l=5,m=2
                ylma[29] = tmp0 * ylma[20] - ylmcoef[29] * ylma[13]; // l=5,m=-2

                tmp0 = ylmcoef[30] * dr[2];
                ylma[30] = tmp0 * ylma[21] - ylmcoef[31] * ylma[14]; // l=5,m=3
                ylma[31] = tmp0 * ylma[22] - ylmcoef[31] * ylma[15]; // l=5,m=-3

                tmp0 = ylmcoef[32] * dr[2];
                ylma[32] = tmp0 * ylma[23]; // l=5,m=4
                ylma[33] = tmp0 * ylma[24]; // l=5,m=-4

                tmp0 = ylmcoef[35] * dr[0];
                ylma[34] = ylmcoef[33] * ylma[30] - ylmcoef[34] * ylma[14] - tmp0 * ylma[23]; // l=5,m=5
                ylma[35] = ylmcoef[33] * ylma[31] - ylmcoef[34] * ylma[15] - tmp0 * ylma[24]; // l=5,m=-5
                if (nwl == 5)
                    goto YLM_END;
                /*
                // if nwl > 5
                for (int il = 6; il <= nwl; il++)
                {
                    int istart = il * il;
                    int istart1 = (il - 1) * (il - 1);
                    int istart2 = (il - 2) * (il - 2);

                    double fac2 = sqrt(4.0 * istart - 1.0);
                    double fac4 = sqrt(4.0 * istart1 - 1.0);

                    for (int im = 0; im < 2 * il - 1; im++)
                    {
                        int imm = (im + 1) / 2;
                        ylma[istart + im] = fac2 / sqrt((double)istart - imm * imm) * (dr[2] * ylma[istart1 + im] - sqrt((double)istart1 - imm * imm) / fac4 * ylma[istart2 + im]);
                    }

                    double bl1 = sqrt(2.0 * il / (2.0 * il + 1.0));
                    double bl2 = sqrt((2.0 * il - 2.0) / (2.0 * il - 1.0));
                    double bl3 = sqrt(2.0) / fac2;

                    ylma[istart + 2 * il - 1] = (bl3 * ylma[istart + 2 * il - 5] - bl2 * ylma[istart2 + 2 * il - 5] - 2.0 * dr[0] * ylma[istart1 + 2 * il - 3]) / bl1;
                    ylma[istart + 2 * il] = (bl3 * ylma[istart + 2 * il - 4] - bl2 * ylma[istart2 + 2 * il - 4] - 2.0 * dr[0] * ylma[istart1 + 2 * il - 2]) / bl1;
                }*/
            }
        YLM_END:
            distance /= delta_r_g[0];

            int ip = (int)(distance);
            double dx = distance - ip;
            double dx2 = dx * dx;
            double dx3 = dx2 * dx;

            double c3 = 3.0 * dx2 - 2.0 * dx3;
            double c1 = 1.0 - c3;
            double c2 = (dx - 2.0 * dx2 + dx3) * delta_r_g[0];
            double c4 = (dx3 - dx2) * delta_r_g[0];

            double phi = 0.0;
            int it_nw = it * nwmax_g[0];
            int iw_nr = (it_nw * nr_max + ip) * 2;
            int it_nw_iw = it_nw;
            for (int iw = 0; iw < atom_nw[it]; ++iw)
            {
                if (atom_iw2_new[it_nw_iw])
                {
                    phi = c1 * psi_u[iw_nr] + c2 * psi_u[iw_nr + 1] + c3 * psi_u[iw_nr + 2] + c4 * psi_u[iw_nr + 3];
                }
                double temp = phi * ylma[atom_iw2_ylm[it_nw_iw]];
                psir_ylm_right[dist_tmp] = temp;
                psir_ylm_left[dist_tmp] = temp * vlbr3_value;
                dist_tmp++;
                iw_nr += nr_max;
                iw_nr += nr_max;
                it_nw_iw++;
            }
        }
    }
}

__global__ void psi_multiple(int *atom_pair_input_info_g,
                             int *num_atom_pair_g,
                             int grid_index,
                             double *psir_ylm_left,
                             double *psir_ylm_right,
                             int atom_pair_size_of_meshcell,
                             double *GridVlocal,
                             int lgd)
{
    //int k = blockIdx.x;
    grid_index += blockIdx.x;
    int atom_pair_num = num_atom_pair_g[blockIdx.x];
    int start_index = atom_pair_size_of_meshcell * blockIdx.x;
    int end_index = start_index + atom_pair_num;
    start_index += blockIdx.y * 6;
    int step = gridDim.y * 6;
    for (int atom_pair_index = start_index; atom_pair_index < end_index; atom_pair_index += step)
    {
        int atomnow1 = atom_pair_input_info_g[atom_pair_index];
        int atomnow2 = atom_pair_input_info_g[atom_pair_index + 1];
        int iw1 = threadIdx.x;
        int iw2 = threadIdx.y;
        if (iw1 >= atom_pair_input_info_g[atom_pair_index + 2] || iw2 >= atom_pair_input_info_g[atom_pair_index + 3])
        {
            return;
        }
        int lo1_iw1 = atom_pair_input_info_g[atom_pair_index + 4] + iw1;
        int lo2_iw2 = atom_pair_input_info_g[atom_pair_index + 5] + iw2;
        double v2 = 0.0;
        int vldr3_index = blockIdx.x * bxyz_g[0];

        for (int ib = 0; ib < bxyz_g[0]; ++ib)
        {
            int calc_index1 = vldr3_index * max_size_g[0];
            int calc_index2 = calc_index1 + atomnow2;
            calc_index1 += atomnow1;
            v2 += psir_ylm_left[calc_index1 * nwmax_g[0] + iw1] * psir_ylm_right[calc_index2 * nwmax_g[0] + iw2];
            vldr3_index++;
        }
        atomicAdd(&(GridVlocal[lo1_iw1 * lgd + lo2_iw2]), v2);
    }
}

void gint_gamma_vl_gpu(hamilt::HContainer<double>* hRGint,
                       const int lgd,
                       const int max_size,
                       const double vfactor,
                       const double *vlocal,
                       const double *ylmcoef_now,
                       const int bx,
                       const int by,
                       const int bz,
                       const int bxyz,
                       const int ncx,
                       const int ncy,
                       const int nczp,
                       const int NLOCAL_now,
                       const int nbxx,
                       int *start_ind,
                       const Grid_Technique &GridT)
{
    // printf("\n**************START GPU SEG***************\n");
#ifdef __DEBUG
    cudaEvent_t t1, t2, t3, t4;
    cudaEventCreate(&t1);
    cudaEventCreate(&t2);
    cudaEventCreate(&t3);
    cudaEventCreate(&t4);

    cudaEventRecord(t1);
#endif

    const Numerical_Orbital_Lm *pointer;
    // const double delta_r = GlobalC::ORB.dr_uniform;
    // const int total_atoms_on_grid = GridT.total_atoms_on_grid;
    const int nbx = GridT.nbx;
    const int nby = GridT.nby;
    const int nbz = GridT.nbzp;
    const int nwmax = GlobalC::ucell.nwmax;
    const int ntype = GlobalC::ucell.ntype;

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

    cudaMemcpyToSymbol(ylmcoef, ylmcoef_now, 36 * sizeof(double));
    cudaMemcpyToSymbol(bxyz_g, &bxyz, sizeof(int));
    cudaMemcpyToSymbol(max_size_g, &max_size, sizeof(int));
    cudaMemcpyToSymbol(nwmax_g, &nwmax, sizeof(int));
    cudaMemcpyToSymbol(delta_r_g, &GlobalC::ORB.dr_uniform, sizeof(double));

    int *atom_nw_g;
    cudaMalloc((void **)&atom_nw_g, ntype * sizeof(int));
    cudaMemcpy(atom_nw_g, atom_nw_now, ntype * sizeof(int), cudaMemcpyHostToDevice);

    int *ucell_atom_nwl;
    cudaMalloc((void **)&ucell_atom_nwl, ntype * sizeof(int));
    cudaMemcpy(ucell_atom_nwl, ucell_atom_nwl_now, ntype * sizeof(int), cudaMemcpyHostToDevice);

    double *psi_u;
    cudaMalloc((void **)&psi_u, ntype * nwmax * nr_max * sizeof(double) * 2);
    cudaMemcpy(psi_u, psi_u_now, ntype * nwmax * nr_max * sizeof(double) * 2, cudaMemcpyHostToDevice);

    bool *atom_iw2_new;
    int *atom_iw2_ylm;
    cudaMalloc((void **)&atom_iw2_new, ntype * nwmax * sizeof(bool));
    cudaMalloc((void **)&atom_iw2_ylm, ntype * nwmax * sizeof(int));
    cudaMemcpy(atom_iw2_new, atom_iw2_new_now, ntype * nwmax * sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(atom_iw2_ylm, atom_iw2_ylm_now, ntype * nwmax * sizeof(int), cudaMemcpyHostToDevice);

    double *psir_ylm_left;
    cudaMalloc((void **)&psir_ylm_left, nbz * max_size * bxyz * nwmax * sizeof(double));
    cudaMemset(psir_ylm_left, 0, nbz * max_size * bxyz * nwmax * sizeof(double));

    double *psir_ylm_right;
    cudaMalloc((void **)&psir_ylm_right, nbz * max_size * bxyz * nwmax * sizeof(double));
    cudaMemset(psir_ylm_right, 0, nbz * max_size * bxyz * nwmax * sizeof(double));

    double *GridVlocal_now = new double[lgd * lgd];

    double *GridVlocal;
    cudaMalloc((void **)&GridVlocal, lgd * lgd * sizeof(double));
    cudaMemset(GridVlocal, 0, lgd * lgd * sizeof(double));

    const int atom_pair_size_of_meshcell = max_size * max_size * 6;
    const int atom_pair_size_over_nbz = atom_pair_size_of_meshcell * nbz;

    int *atom_pair_input_info = new int[atom_pair_size_over_nbz];
    int *atom_pair_input_info_g;
    cudaMalloc((void **)&atom_pair_input_info_g, atom_pair_size_over_nbz * sizeof(int));

    int *num_atom_pair = new int[nbz];
    int *num_atom_pair_g;
    cudaMalloc((void **)&num_atom_pair_g, nbz * sizeof(int));

    const int psi_size_max = max_size * bxyz;
    double *psi_input_double = new double[psi_size_max * nbz * 5]; // [ x,y,z,distance, vlocal]
    double *psi_input_double_g;
    cudaMalloc((void **)&psi_input_double_g, psi_size_max * nbz * 5 * sizeof(double));

    int *psi_input_int = new int[psi_size_max * nbz * 2];
    int *psi_input_int_g;
    cudaMalloc((void **)&psi_input_int_g, psi_size_max * nbz * sizeof(int) * 2);

    int *num_psir = new int[nbz];
    int *num_psir_g;
    cudaMalloc((void **)&num_psir_g, nbz * sizeof(int));

#ifdef __DEBUG

    cudaEventRecord(t2);

    float copy_per_calc = 0;
    float calc_psi = 0;
    float calc_multiple = 0;
#endif
    for (int i = 0; i < nbx; i++)
    {
        for (int j = 0; j < nby; j++)
        {
#ifdef __DEBUG
            cudaEvent_t t1_5, t1_6, t1_7, t1_8;
            cudaEventCreate(&t1_5);
            cudaEventCreate(&t1_6);
            cudaEventCreate(&t1_7);
            cudaEventCreate(&t1_8);

            cudaEventRecord(t1_5);
#endif
            gpu_task_generate_vlocal(GridT, i, j, bx,
                                        by, bz, bxyz,
                                        atom_pair_size_of_meshcell, psi_size_max,
                                        max_size, ncx, ncy, nczp,
                                        vfactor,
                                        start_ind,
                                        vlocal, psi_input_double, psi_input_int,
                                        num_psir, atom_pair_input_info, num_atom_pair);

            cudaMemcpy(psi_input_double_g, psi_input_double, psi_size_max * nbz * 5 * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(psi_input_int_g, psi_input_int, psi_size_max * nbz * 2 * sizeof(int), cudaMemcpyHostToDevice);

            cudaMemcpy(num_psir_g, num_psir, nbz * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemset(psir_ylm_left, 0, nbz * max_size * bxyz * nwmax * sizeof(double));
            cudaMemset(psir_ylm_right, 0, nbz * max_size * bxyz * nwmax * sizeof(double));

            cudaMemcpy(atom_pair_input_info_g, atom_pair_input_info, atom_pair_size_over_nbz * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(num_atom_pair_g, num_atom_pair, nbz * sizeof(int), cudaMemcpyHostToDevice);
#ifdef __DEBUG
            cudaEventRecord(t1_6);
            cudaDeviceSynchronize();
#endif
            const int ALIGN_SIZE = 32;
            dim3 grid1(nbz);
            dim3 block1(ALIGN_SIZE);
            int shared_size = bxyz;
            get_psi_and_vldr3<<<grid1, block1, shared_size>>>(psi_input_double_g,
                                                                psi_input_int_g,
                                                                num_psir_g,
                                                                psi_size_max,
                                                                ucell_atom_nwl,
                                                                atom_iw2_new,
                                                                atom_iw2_ylm,
                                                                atom_nw_g,
                                                                nr_max,
                                                                psi_u,
                                                                psir_ylm_left,
                                                                psir_ylm_right);
#ifdef __DEBUG
            cudaDeviceSynchronize();
            cudaEventRecord(t1_7);
            cudaDeviceSynchronize();
#endif
            dim3 grid4(nbz, 256);
            dim3 block4(nwmax, nwmax);
            psi_multiple<<<grid4, block4>>>(atom_pair_input_info_g,
                                            num_atom_pair_g,
                                            i * nby * nbz + j * nbz,
                                            psir_ylm_left,
                                            psir_ylm_right,
                                            atom_pair_size_of_meshcell,
                                            GridVlocal,
                                            lgd);

#ifdef __DEBUG
            cudaDeviceSynchronize();
            cudaEventRecord(t1_8);
            float copy_per_calc_temp = 0;
            float calc_psi_temp = 0;
            float calc_multiple_temp = 0;
            cudaDeviceSynchronize();

            cudaEventElapsedTime(&copy_per_calc_temp, t1_5, t1_6);
            cudaEventElapsedTime(&calc_psi_temp, t1_6, t1_7);
            cudaEventElapsedTime(&calc_multiple_temp, t1_7, t1_8);
            copy_per_calc += copy_per_calc_temp;
            calc_psi += calc_psi_temp;
            calc_multiple += calc_multiple_temp;
#endif

        }
    }

    cudaMemcpy(GridVlocal_now, GridVlocal, lgd * lgd * sizeof(double), cudaMemcpyDeviceToHost);

    for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
    {
        for (int iat2 = 0; iat2 < GlobalC::ucell.nat; iat2++)
        {
            int it1 = GlobalC::ucell.iat2it[iat1];
            int it2 = GlobalC::ucell.iat2it[iat2];
            int lo1 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(it1, GlobalC::ucell.iat2ia[iat1], 0)];
            int lo2 = GridT.trace_lo[GlobalC::ucell.itiaiw2iwt(it2, GlobalC::ucell.iat2ia[iat2], 0)];
            if (lo1 <= lo2)
            {
				hamilt::AtomPair<double>* tmp_ap = hRGint->find_pair(iat1, iat2);
                int orb_index = 0;
                if (tmp_ap == NULL) continue;
                for(int orb_i = 0; orb_i < tmp_ap->get_row_size();orb_i++)
                {
                    for(int orb_j = 0;orb_j < tmp_ap->get_col_size();orb_j++)
                    {
                        tmp_ap->get_pointer(0)[orb_index] = GridVlocal_now[(lo1 + orb_i) * lgd + (lo2 + orb_j)];
                        orb_index++;
                    }
                }
            }
        }
    }
#ifdef __DEBUG

    cudaEventRecord(t3);
#endif

    // free
    cudaFree(psir_ylm_left);
    cudaFree(psir_ylm_right);
    cudaFree(atom_nw_g);

    cudaFree(ucell_atom_nwl);
    cudaFree(psi_u);
    cudaFree(atom_iw2_new);
    cudaFree(atom_iw2_ylm);
    cudaFree(GridVlocal);

    cudaFree(atom_pair_input_info_g);
    cudaFree(num_atom_pair_g);

    cudaFree(psi_input_double_g);
    cudaFree(psi_input_int_g);
    cudaFree(num_psir_g);

    delete[] atom_pair_input_info;
    delete[] num_atom_pair;

    delete[] psi_input_double;
    delete[] psi_input_int;
    delete[] num_psir;

    delete[] GridVlocal_now;
    delete[] atom_nw_now;
    delete[] ucell_atom_nwl_now;
    delete[] psi_u_now;
    delete[] atom_iw2_new_now;
    delete[] atom_iw2_ylm_now;

#ifdef __DEBUG
    cudaEventRecord(t4);
    float copy = 0;
    float calc = 0;
    float free = 0;
    cudaEventElapsedTime(&copy, t1, t2);
    cudaEventElapsedTime(&calc, t2, t3);
    cudaEventElapsedTime(&free, t3, t4);

    printf("copy time = %f\ncal time = %f\nfree time = %f\n", copy, calc, free);
    printf("copy_per calc time = %f\ncal psi time = %f\nmultiple time = %f\n", copy_per_calc, calc_psi, calc_multiple);
#endif

}