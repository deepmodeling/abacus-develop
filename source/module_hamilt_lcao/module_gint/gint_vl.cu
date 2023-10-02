#include "omp.h"
#include "gint_tools.h"
#include "gint_vl.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include <fstream>
#include <sstream>

__constant__ double ylmcoef[36];
__constant__ int bx_g[1];
__constant__ int by_g[1];
__constant__ int bz_g[1];
__constant__ int bxyz_g[1];
__constant__ int max_size_g[1];
__constant__ int nwmax_g[1];
__constant__ int namax_g[1];
__constant__ int nnnmax_g[1];
__constant__ int ntype_g[1];
__constant__ double delta_r_g[1];
__constant__ double vfactor_g[1];

void dump_cuda_array_to_file(double * cuda_array, int width, int hight, const std::string & filename)
{
    double *h_data = new double[width * hight];
    cudaMemcpy(h_data, cuda_array, width * hight * sizeof(double), cudaMemcpyDeviceToHost);

    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file for writing." << std::endl;
    }
    for (int j = 0; j < hight; ++j) {
        for (int i = 0; i < width; ++i) {
            outFile << "hight" << j << "   width:" << i << "   "<< h_data[j * width + i] << std::endl;
        }
    }
    outFile.close();
    delete[] h_data;
}

__global__ void get_psi(double *dr_all,
                        int *it_all,
                        int *psir_ylm_start,
                        int *num_psir,
                        int psi_size_up,
                        bool *cal_flag,
                        int *ucell_atom_nwl,
                        bool *atom_iw2_new,
                        int *atom_iw2_ylm,
                        int *atom_nw,
                        int nr_max,
                        double *psi_u,
                        double *dpsi_u,
                        double *psir_ylm)
{
    // int grid_index_now = blockIdx.x;
    int size = num_psir[blockIdx.x];
    int loop_size = size / blockDim.x;

    if (loop_size != 0)
    {
        int start_index = loop_size * threadIdx.x + psi_size_up * blockIdx.x;
        int end_index = start_index + loop_size;
        for (int index = start_index; index < end_index; index++)
        {
            int it = it_all[index];
            if (it < 0) continue;

            double dr[3];
            dr[0] = dr_all[index * 4];
            dr[1] = dr_all[index * 4 + 1];
            dr[2] = dr_all[index * 4 + 2];
            double distance = dr_all[index * 4 + 3];
            int dist_tmp = psir_ylm_start[index];

            cal_flag[dist_tmp] = true;

            // begin calculation
            double ylma[150];
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
                //	ylma[grid_index*nnnmax+8] = tmp1+tmp2*ylma[grid_index*nnnmax+3];//l=2,m=-2
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
                        //			if (im % 2 == 0) imm *= -1;

                        ylma[istart + im] = fac2 / sqrt((double)istart - imm * imm) * (dr[2] * ylma[istart1 + im] - sqrt((double)istart1 - imm * imm) / fac4 * ylma[istart2 + im]);
                    }

                    double bl1 = sqrt(2.0 * il / (2.0 * il + 1.0));
                    double bl2 = sqrt((2.0 * il - 2.0) / (2.0 * il - 1.0));
                    double bl3 = sqrt(2.0) / fac2;

                    ylma[istart + 2 * il - 1] = (bl3 * ylma[istart + 2 * il - 5] - bl2 * ylma[istart2 + 2 * il - 5] - 2.0 * dr[0] * ylma[istart1 + 2 * il - 3]) / bl1;
                    ylma[istart + 2 * il] = (bl3 * ylma[istart + 2 * il - 4] - bl2 * ylma[istart2 + 2 * il - 4] - 2.0 * dr[0] * ylma[istart1 + 2 * il - 2]) / bl1;
                }
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

            int iw;
            double phi = 0.0;
            for (iw = 0; iw < atom_nw[it]; ++iw)
            {
                if (atom_iw2_new[it * nwmax_g[0] + iw])
                {
                    phi = c1 * psi_u[it * nwmax_g[0] * nr_max + iw * nr_max + ip] + c2 * dpsi_u[it * nwmax_g[0] * nr_max + iw * nr_max + ip] + c3 * psi_u[it * nwmax_g[0] * nr_max + iw * nr_max + ip + 1] + c4 * dpsi_u[it * nwmax_g[0] * nr_max + iw * nr_max + ip + 1];
                }
                psir_ylm[dist_tmp * nwmax_g[0] + iw] = phi * ylma[atom_iw2_ylm[it * nwmax_g[0] + iw]];
            }
        }
    } // if size
}
__global__ void get_vldr3(int ij_index,
                          int nbx,
                          int nby,
                          int nbz,
                          int nbz_start,
                          int ncy,
                          int nczp,
                          double *vlocal,
                          int *start_ind_g,
                          double *vldr3)
{
    int k = blockIdx.x;
    int ii = threadIdx.x;
    int jj = threadIdx.y;
    int kk = threadIdx.z;
    int vindex = ii * ncy * nczp + jj * nczp + kk + start_ind_g[ij_index + k];
    vldr3[k * bx_g[0] * by_g[0] * bz_g[0] + ii * by_g[0] * bz_g[0] + jj * bz_g[0] + kk] = vlocal[vindex] * vfactor_g[0];
}

__global__ void psi_multiple(int grid_index,
                             int k,
                             int *how_many_atoms,
                             double *atom_pair_index1_g,
                             double *atom_pair_index2_g,
                             int *bcell_start,
                             int *which_atom,
                             int *iat2it,
                             int *iat2ia,
                             int *itiaiw2iwt,
                             bool *cal_flag,
                             double *psir_ylm,
                             int *trace_lo,
                             int *atom_nw,
                             double *vldr3,
                             double *GridVlocal,
                             int lgd)
{
    int atomnow1 = blockIdx.x;
    int atomnow2 = blockIdx.y;
    int iw1 = threadIdx.x;
    int iw2 = threadIdx.y;
    if (atomnow1 >= how_many_atoms[grid_index] || atomnow2 >= how_many_atoms[grid_index])
    {
        return;
    }
    int iat1 = which_atom[bcell_start[grid_index] + atomnow1];
    int iat2 = which_atom[bcell_start[grid_index] + atomnow2];
    int it1 = iat2it[iat1];
    int it2 = iat2it[iat2];
    if (iw1 >= atom_nw[it1] || iw2 >= atom_nw[it2])
    {
        return;
    }

    int lo1 = trace_lo[itiaiw2iwt[it1 * namax_g[0] + iat2ia[iat1]]];
    int lo2 = trace_lo[itiaiw2iwt[it2 * namax_g[0] + iat2ia[iat2]]];
    if (lo1 <= lo2)
    {
        int lo1_iw1 = lo1 + iw1;
        int lo2_iw2 = lo2 + iw2;
        double v2 = 0.0;
        for (int ib = 0; ib < bxyz_g[0]; ++ib)
        {
            int vldr3_index = k * bxyz_g[0] + ib;
            int calc_index1 = vldr3_index * max_size_g[0];
            int calc_index2 = calc_index1 + atomnow2;
            calc_index1 += atomnow1;
            if (cal_flag[calc_index1] &&
                cal_flag[calc_index2])
            {
                v2 += psir_ylm[calc_index1 * nwmax_g[0] + iw1] * vldr3[vldr3_index] * psir_ylm[calc_index2 * nwmax_g[0] + iw2];
            }
        }
        atomicAdd(&(GridVlocal[lo1_iw1 * lgd + lo2_iw2]), v2);
    }
}

void gint_gamma_vl_gpu(double *GridVlocal_now,
                       int lgd_now,
                       int nnnmax,
                       int max_size,
                       double vfactor,
                       const double *vlocal,
                       const double *ylmcoef_now,
                       int pwbx,
                       int pwby,
                       int pwbz,
                       int pwbxyz,
                       int pwncx,
                       int pwncy,
                       int pwnczp,
                       int NLOCAL_now,
                       int nbxx,
                       int *start_ind,
                       const Grid_Technique &GridT)
{
    // printf("\n**************START GPU SEG***************\n");

    cudaEvent_t t1, t2, t3, t4;
    cudaEventCreate(&t1);
    cudaEventCreate(&t2);
    cudaEventCreate(&t3);
    cudaEventCreate(&t4);

    cudaEventRecord(t1);

    const Numerical_Orbital_Lm *pointer;
    const double delta_r = GlobalC::ORB.dr_uniform;
    const int total_atoms_on_grid = GridT.total_atoms_on_grid;
    const int nbx = GridT.nbx;
    const int nby = GridT.nby;
    const int nbz_start = GridT.nbzp_start;
    const int nbz = GridT.nbzp;
    const int bx = pwbx;
    const int by = pwby;
    const int bz = pwbz;
    const int bxyz_now = pwbxyz;
    const int ncx = pwncx;
    const int ncy = pwncy;
    const int nczp = pwnczp;
    const int nwmax_now = GlobalC::ucell.nwmax;
    const int namax_now = GlobalC::ucell.namax;
    const int nype_now = GlobalC::ucell.ntype;

    size_t size_phi = GlobalC::ucell.ntype;
    double max_cut = 0;
    for (int i = 0; i < size_phi; i++)
    {
        if (GlobalC::ORB.Phi[i].getRcut() > max_cut)
        {
            max_cut = GlobalC::ORB.Phi[i].getRcut();
        }
    }

    size_t size_atom_nw = GlobalC::ucell.ntype;
    int *atom_nw_now;
    int *ucell_atom_nwl_now;
    atom_nw_now = new int[size_atom_nw];
    ucell_atom_nwl_now = new int[size_atom_nw];
    for (int i = 0; i < size_atom_nw; i++)
    {
        atom_nw_now[i] = GlobalC::ucell.atoms[i].nw;
        ucell_atom_nwl_now[i] = GlobalC::ucell.atoms[i].nwl;
    }

    int nr_max = static_cast<int>(1000 * max_cut) + 10;
    double *psi_u_now = new double[nype_now * nwmax_now * nr_max];
    double *dpsi_u_now = new double[nype_now * nwmax_now * nr_max];
    bool *atom_iw2_new_now = new bool[nype_now * nwmax_now];
    int *atom_iw2_ylm_now = new int[nype_now * nwmax_now];

    Atom *atomx;
    for (int i = 0; i < nype_now; i++)
    {
        atomx = &GlobalC::ucell.atoms[i];
        for (int j = 0; j < nwmax_now; j++)
        {
            if (j < atomx->nw)
            {

                atom_iw2_new_now[i * nwmax_now + j] = atomx->iw2_new[j];
                atom_iw2_ylm_now[i * nwmax_now + j] = atomx->iw2_ylm[j];
                pointer = &GlobalC::ORB.Phi[i].PhiLN(atomx->iw2l[j], atomx->iw2n[j]);
                for (int k = 0; k < nr_max; k++)
                {
                    if (k < pointer->nr_uniform)
                    {
                        psi_u_now[i * nwmax_now * nr_max + j * nr_max + k] = pointer->psi_uniform[k];
                        dpsi_u_now[i * nwmax_now * nr_max + j * nr_max + k] = pointer->dpsi_uniform[k];
                    }
                    else
                    {
                        psi_u_now[i * nwmax_now * nr_max + j * nr_max + k] = 0;
                        dpsi_u_now[i * nwmax_now * nr_max + j * nr_max + k] = 0;
                    }
                }
            }
            else
            {

                atom_iw2_new_now[i * nwmax_now + j] = false;
                atom_iw2_ylm_now[i * nwmax_now + j] = 0;
                for (int k = 0; k < nr_max; k++)
                {
                    psi_u_now[i * nwmax_now * nr_max + j * nr_max + k] = 0;
                    dpsi_u_now[i * nwmax_now * nr_max + j * nr_max + k] = 0;
                }
            }
        }
    }

    size_t size_itiaiw2iwt = nype_now * namax_now;
    int *itiaiw2iwt_now = new int[size_itiaiw2iwt];
    for (int i = 0; i < nype_now; i++)
    {
        for (int j = 0; j < namax_now; j++)
        {
            itiaiw2iwt_now[i * namax_now + j] = GlobalC::ucell.itiaiw2iwt(i, j, 0);
        }
    }

    cudaMemcpyToSymbol(ylmcoef, ylmcoef_now, 36 * sizeof(double));
    cudaMemcpyToSymbol(bx_g, &bx, sizeof(int));
    cudaMemcpyToSymbol(by_g, &by, sizeof(int));
    cudaMemcpyToSymbol(bz_g, &bz, sizeof(int));
    cudaMemcpyToSymbol(bxyz_g, &bxyz_now, sizeof(int));
    cudaMemcpyToSymbol(max_size_g, &max_size, sizeof(int));
    cudaMemcpyToSymbol(nwmax_g, &nwmax_now, sizeof(int));
    cudaMemcpyToSymbol(namax_g, &namax_now, sizeof(int));
    cudaMemcpyToSymbol(nnnmax_g, &nnnmax, sizeof(int));
    cudaMemcpyToSymbol(ntype_g, &nype_now, sizeof(int));
    cudaMemcpyToSymbol(delta_r_g, &delta_r, sizeof(double));
    cudaMemcpyToSymbol(vfactor_g, &vfactor, sizeof(double));

    // read only
    int *how_many_atoms;
    cudaError_t status = cudaMalloc((void **)&how_many_atoms, nbx * nby * nbz * sizeof(int));
    if (status != cudaSuccess)
    {
        fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(status));
        return;
    }

    cudaError_t status2 = cudaMemcpy(how_many_atoms, GridT.how_many_atoms, nbx * nby * nbz * sizeof(int), cudaMemcpyHostToDevice);
    if (status2 != cudaSuccess)
    {
        fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(status2));
        return;
    }

    int *bcell_start;
    cudaMalloc((void **)&bcell_start, nbx * nby * nbz * sizeof(int));
    cudaMemcpy(bcell_start, GridT.bcell_start, nbx * nby * nbz * sizeof(int), cudaMemcpyHostToDevice);

    int *which_bigcell;
    cudaMalloc((void **)&which_bigcell, total_atoms_on_grid * sizeof(int));
    cudaMemcpy(which_bigcell, GridT.which_bigcell, total_atoms_on_grid * sizeof(int), cudaMemcpyHostToDevice);

    int *which_atom;
    cudaMalloc((void **)&which_atom, total_atoms_on_grid * sizeof(int));
    cudaMemcpy(which_atom, GridT.which_atom, total_atoms_on_grid * sizeof(int), cudaMemcpyHostToDevice);

    int *iat2it;
    size_t size_iat2it = GlobalC::ucell.nat;
    cudaMalloc((void **)&iat2it, size_iat2it * sizeof(int));
    cudaMemcpy(iat2it, GlobalC::ucell.iat2it, size_iat2it * sizeof(int), cudaMemcpyHostToDevice);

    int *iat2ia;
    size_t size_iat2ia = GlobalC::ucell.nat;
    cudaMalloc((void **)&iat2ia, size_iat2ia * sizeof(int));
    cudaMemcpy(iat2ia, GlobalC::ucell.iat2ia, size_iat2ia * sizeof(int), cudaMemcpyHostToDevice);

    double *vlocal_cu;
    cudaMalloc((void **)&vlocal_cu, ncx * ncy * nczp * sizeof(double));
    cudaMemcpy(vlocal_cu, vlocal, ncx * ncy * nczp * sizeof(double), cudaMemcpyHostToDevice);

    int *atom_nw;
    cudaMalloc((void **)&atom_nw, size_atom_nw * sizeof(int));
    cudaMemcpy(atom_nw, atom_nw_now, size_atom_nw * sizeof(int), cudaMemcpyHostToDevice);

    int *ucell_atom_nwl;
    cudaMalloc((void **)&ucell_atom_nwl, size_atom_nw * sizeof(int));
    cudaMemcpy(ucell_atom_nwl, ucell_atom_nwl_now, size_atom_nw * sizeof(int), cudaMemcpyHostToDevice);

    double *psi_u;
    double *dpsi_u;
    cudaMalloc((void **)&psi_u, nype_now * nwmax_now * nr_max * sizeof(double));
    cudaMalloc((void **)&dpsi_u, nype_now * nwmax_now * nr_max * sizeof(double));
    cudaMemcpy(psi_u, psi_u_now, nype_now * nwmax_now * nr_max * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dpsi_u, dpsi_u_now, nype_now * nwmax_now * nr_max * sizeof(double), cudaMemcpyHostToDevice);

    bool *atom_iw2_new;
    int *atom_iw2_ylm;
    cudaMalloc((void **)&atom_iw2_new, nype_now * nwmax_now * sizeof(bool));
    cudaMalloc((void **)&atom_iw2_ylm, nype_now * nwmax_now * sizeof(int));
    cudaMemcpy(atom_iw2_new, atom_iw2_new_now, nype_now * nwmax_now * sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(atom_iw2_ylm, atom_iw2_ylm_now, nype_now * nwmax_now * sizeof(int), cudaMemcpyHostToDevice);

    int *trace_lo;
    size_t size_trace_lo = NLOCAL_now;
    cudaMalloc((void **)&trace_lo, size_trace_lo * sizeof(int));
    cudaMemcpy(trace_lo, GridT.trace_lo, size_trace_lo * sizeof(int), cudaMemcpyHostToDevice);

    int *itiaiw2iwt;
    cudaMalloc((void **)&itiaiw2iwt, size_itiaiw2iwt * sizeof(int));
    cudaMemcpy(itiaiw2iwt, itiaiw2iwt_now, size_itiaiw2iwt * sizeof(int), cudaMemcpyHostToDevice);

    int *start_ind_g;
    cudaMalloc((void **)&start_ind_g, nbxx * sizeof(int));
    cudaMemcpy(start_ind_g, start_ind, nbxx * sizeof(int), cudaMemcpyHostToDevice);

    double *vldr3;
    cudaMalloc((void **)&vldr3, nbz * bxyz_now * sizeof(double));
    cudaMemset(vldr3, 0, nbz * bxyz_now * sizeof(double));

    double *psir_ylm;
    cudaMalloc((void **)&psir_ylm, nbz * max_size * bxyz_now * nwmax_now * sizeof(double));
    cudaMemset(psir_ylm, 0, nbz * max_size * bxyz_now * nwmax_now * sizeof(double));

    bool *cal_flag;
    cudaMalloc((void **)&cal_flag, nbz * bxyz_now * max_size * sizeof(bool));
    cudaMemset(cal_flag, 0, nbz * bxyz_now * max_size * sizeof(bool));

    double *GridVlocal;
    cudaMalloc((void **)&GridVlocal, lgd_now * lgd_now * sizeof(double));
    cudaMemset(GridVlocal, 0, lgd_now * lgd_now * sizeof(double));

    const size_t atom_pair_num = (max_size * (max_size + 1)) / 2;

    double *atom_pair_index1 = new double[atom_pair_num];
    double *atom_pair_index2 = new double[atom_pair_num];
    for (int x = 0; x < max_size; x++)
    {
        for (int y = 0; y < max_size - x; y++)
        {
            int index = ((max_size * 2 - x + 1) * x) / 2 + y;
            atom_pair_index1[index] = x;
            atom_pair_index2[index] = y;
        }
    }
    double *atom_pair_index1_g;
    cudaMalloc((void **)&atom_pair_index1_g, atom_pair_num * sizeof(double));

    double *atom_pair_index2_g;
    cudaMalloc((void **)&atom_pair_index2_g, atom_pair_num * sizeof(double));

    cudaMemcpy(atom_pair_index1_g, atom_pair_index1, atom_pair_num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(atom_pair_index2_g, atom_pair_index2, atom_pair_num * sizeof(double), cudaMemcpyHostToDevice);

    int psi_size = max_size * bxyz_now;
    int ALIGN_SIZE = 32;
    int psi_size_up = ((psi_size + ALIGN_SIZE - 1) / ALIGN_SIZE) * ALIGN_SIZE;
    double *dr = new double[psi_size_up * nbz * 4]; // [ x,y,z,distance]
    int *it = new int[psi_size_up * nbz];
    int *psir_ylm_start = new int[psi_size_up * nbz];
    int *num_psir = new int[nbz];
    // begin kernel

    double *dr_g; // [ x,y,z,distance]
    cudaMalloc((void **)&dr_g, psi_size_up * nbz * 4 * sizeof(double));

    int *it_g;
    cudaMalloc((void **)&it_g, psi_size_up * nbz * sizeof(int));

    int *psir_ylm_start_g;
    cudaMalloc((void **)&psir_ylm_start_g, psi_size_up * nbz * sizeof(int));

    int *num_psir_g;
    cudaMalloc((void **)&num_psir_g, nbz * sizeof(int));

    cudaEventRecord(t2);

    // printf("maxsize=%d\n", max_size);

    for (int i = 0; i < nbx; i++)
    {
        for (int j = 0; j < nby; j++)
        {
            int num_psi_pos = 0;
            for (int z_index = 0; z_index < nbz; z_index++)
            {
                int num_get_psi = 0;
                int grid_index = i * nby * nbz + j * nbz + z_index;
                for (int id = 0; id < GridT.how_many_atoms[grid_index]; id++)
                {
                    for (int ib = 0; ib < bxyz_now; ib++)
                    {
                        int mcell_index = GridT.bcell_start[grid_index] + id;
                        int imcell = GridT.which_bigcell[mcell_index];
                        int iat = GridT.which_atom[mcell_index];
                        int it_temp = GlobalC::ucell.iat2it[iat];
                        double dr_temp[3];
                        dr_temp[0] = GridT.meshcell_pos[ib][0] + GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
                        dr_temp[1] = GridT.meshcell_pos[ib][1] + GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
                        dr_temp[2] = GridT.meshcell_pos[ib][2] + GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];

                        double distance = sqrt(dr_temp[0] * dr_temp[0] + dr_temp[1] * dr_temp[1] + dr_temp[2] * dr_temp[2]);
                        if (distance <= GlobalC::ORB.Phi[it_temp].getRcut())
                        {
                            int pos_temp = num_psi_pos + num_get_psi;
                            if (distance < 1.0E-9) distance += 1.0E-9;
                            dr[pos_temp * 4] = dr_temp[0] / distance;
                            dr[pos_temp * 4 + 1] = dr_temp[1] / distance;
                            dr[pos_temp * 4 + 2] = dr_temp[2] / distance;
                            dr[pos_temp * 4 + 3] = distance;
                            it[pos_temp] = it_temp;
                            int dist_tmp = z_index * bxyz_now * max_size + ib * max_size + id;
                            psir_ylm_start[pos_temp] = dist_tmp;
                            num_get_psi++;
                        }
                    }
                }
                int num_get_psi_up = ((num_get_psi + ALIGN_SIZE - 1) / ALIGN_SIZE) * ALIGN_SIZE;
                for (;num_get_psi < num_get_psi_up; num_get_psi++)
                {
                    int pos_temp = num_psi_pos + num_get_psi;
                    dr[pos_temp * 4] = 0.0;
                    dr[pos_temp * 4 + 1] = 0.0;
                    dr[pos_temp * 4 + 2] = 0.0;
                    dr[pos_temp * 4 + 3] = 0.0;
                    it[pos_temp] = -1;
                    psir_ylm_start[pos_temp] = psir_ylm_start[pos_temp - 1] + 1;
                }
                num_psir[z_index] = num_get_psi_up; // align to ALIGN_SIZE 32
                num_psi_pos += psi_size_up;
            }

            cudaMemcpy(dr_g, dr, psi_size_up * nbz * 4 * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(it_g, it, psi_size_up * nbz * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(psir_ylm_start_g, psir_ylm_start, psi_size_up * nbz * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(num_psir_g, num_psir, nbz * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemset(cal_flag, 0, nbz * bxyz_now * max_size * sizeof(bool));

            dim3 grid1(nbz);
            dim3 block1(32); // how_many_atoms,bxyz
            get_psi<<<grid1, block1>>>(dr_g,
                                       it_g,
                                       psir_ylm_start_g,
                                       num_psir_g,
                                       psi_size_up,
                                       cal_flag,
                                       ucell_atom_nwl,
                                       atom_iw2_new,
                                       atom_iw2_ylm,
                                       atom_nw,
                                       nr_max,
                                       psi_u,
                                       dpsi_u,
                                       psir_ylm);
            /* std::stringstream filename;
            filename << "psi" << i << "_" << j;
            dump_cuda_array_to_file(psir_ylm, max_size * bxyz_now * nwmax_now , nbz, filename.str());*/
            dim3 grid3(nbz);
            dim3 block3(bx, by, bz);
            get_vldr3<<<grid3, block3>>>(i * nby * nbz + j * nbz,
                                         nbx,
                                         nby,
                                         nbz,
                                         nbz_start,
                                         ncy,
                                         nczp,
                                         vlocal_cu,
                                         start_ind_g,
                                         vldr3);

            dim3 grid4(max_size, max_size);
            dim3 block4(nwmax_now, nwmax_now);
            for (int k = 0; k < nbz; k++)
            {
                psi_multiple<<<grid4, block4>>>(i * nby * nbz + j * nbz + k,
                                                k,
                                                how_many_atoms,
                                                atom_pair_index1_g,
                                                atom_pair_index2_g,
                                                bcell_start,
                                                which_atom,
                                                iat2it,
                                                iat2ia,
                                                itiaiw2iwt,
                                                cal_flag,
                                                psir_ylm,
                                                trace_lo,
                                                atom_nw,
                                                vldr3,
                                                GridVlocal,
                                                lgd_now);
            }
        } // j
    }     // i

    cudaDeviceSynchronize();
    cudaMemcpy(GridVlocal_now, GridVlocal, lgd_now * lgd_now * sizeof(double), cudaMemcpyDeviceToHost);
    // printf("GridVlocal_now[0]: %lf\n", GridVlocal_now[0]);
    cudaEventRecord(t3);
    cudaDeviceSynchronize();
    // free
    cudaFree(vldr3);
    cudaFree(psir_ylm);
    cudaFree(cal_flag);

    cudaFree(how_many_atoms);
    cudaFree(bcell_start);
    cudaFree(which_bigcell);
    cudaFree(which_atom);
    cudaFree(iat2it);
    cudaFree(iat2ia);
    cudaFree(vlocal_cu);
    cudaFree(ucell_atom_nwl);
    cudaFree(psi_u);
    cudaFree(dpsi_u);
    cudaFree(atom_iw2_new);
    cudaFree(atom_iw2_ylm);
    cudaFree(atom_nw);
    cudaFree(trace_lo);
    cudaFree(itiaiw2iwt);
    cudaFree(start_ind_g);
    cudaFree(GridVlocal);

    cudaFree(atom_pair_index1_g);
    cudaFree(atom_pair_index2_g);
    cudaFree(dr_g);
    cudaFree(it_g);
    cudaFree(psir_ylm_start_g);
    cudaFree(num_psir_g);

    delete[] atom_pair_index1;
    delete[] atom_pair_index2;
    delete[] dr;
    delete[] it;
    delete[] psir_ylm_start;
    delete[] num_psir;

    delete[] atom_nw_now;
    delete[] itiaiw2iwt_now;
    delete[] ucell_atom_nwl_now;
    delete[] psi_u_now;
    delete[] dpsi_u_now;
    delete[] atom_iw2_new_now;
    delete[] atom_iw2_ylm_now;

    cudaEventRecord(t4);
    float copy = 0;
    float calc = 0;
    float free = 0;
    cudaEventElapsedTime(&copy, t1, t2);
    cudaEventElapsedTime(&calc, t2, t3);
    cudaEventElapsedTime(&free, t3, t4);

    printf("copy time = %f\ncal time = %f\nfree time = %f\n", copy, calc, free);
}