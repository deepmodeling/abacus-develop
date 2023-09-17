#include "omp.h"
#include "gint_tools.h"
#include "gint_vl.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

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

__global__ void cu_gamma_vlocal_step1(int ij_index,
                                      int *how_many_atoms,
                                      int *bcell_start,
                                      int *which_bigcell,
                                      int *which_atom,
                                      int *iat2it,
                                      int *it,
                                      double *meshball_positions,
                                      double *tau_in_bigcell,
                                      double *meshcell_pos,
                                      double *dr,
                                      double *distance,
                                      double *ORB_Phi_getRcut,
                                      bool *cal_flag,
                                      int *ucell_atom_nwl,
                                      double *ylma,
                                      int *ip,
                                      double *dx,
                                      double *dx2,
                                      double *dx3,
                                      double *c1,
                                      double *c2,
                                      double *c3,
                                      double *c4,
                                      bool *atom_iw2_new,
                                      int *atom_iw2_ylm,
                                      int *atom_nw,
                                      int nr_max,
                                      double *psi_u,
                                      double *dpsi_u,
                                      double *psir_ylm)
{
    // int grid_index_now = blockIdx.x;
    int size = how_many_atoms[ij_index + blockIdx.x];
    if (size != 0)
    {
        double mt[] = {0, 0, 0};
        int id = threadIdx.x;
        int ib = threadIdx.y;
        if (id >= size)
        {
            return;
        }
        int dist_tmp = blockIdx.x * bxyz_g[0] * max_size_g[0] + ib * max_size_g[0] + id;
        int pos_tmp = dist_tmp * 3;
        int mcell_index = bcell_start[ij_index + blockIdx.x] + id;
        int imcell = which_bigcell[mcell_index];
        int iat = which_atom[mcell_index];
        it[dist_tmp] = iat2it[iat];
        mt[0] = meshball_positions[imcell * 3 + 0] - tau_in_bigcell[iat * 3 + 0];
        mt[1] = meshball_positions[imcell * 3 + 1] - tau_in_bigcell[iat * 3 + 1];
        mt[2] = meshball_positions[imcell * 3 + 2] - tau_in_bigcell[iat * 3 + 2];
        dr[pos_tmp + 0] = meshcell_pos[ib * 3 + 0] + mt[0];
        dr[pos_tmp + 1] = meshcell_pos[ib * 3 + 1] + mt[1];
        dr[pos_tmp + 2] = meshcell_pos[ib * 3 + 2] + mt[2];

        distance[dist_tmp] = sqrt(dr[pos_tmp + 0] * dr[pos_tmp + 0] + dr[pos_tmp + 1] * dr[pos_tmp + 1] + dr[pos_tmp + 2] * dr[pos_tmp + 2]);
        if (distance[dist_tmp] > ORB_Phi_getRcut[it[dist_tmp]])
        {
            cal_flag[dist_tmp] = false;
            return;
        }
        else
        {
            cal_flag[dist_tmp] = true;
        }
        if (distance[dist_tmp] < 1.0E-9)
            distance[dist_tmp] += 1.0E-9;

        //	Ylm::sph_harm (	ucell.atoms[it].nwl,dr[grid_index*bxyz_now*max_size_now*3+ib][id][0] /
        // distance[ib][id],dr[grid_index*bxyz_now*max_size_now*3+ib][id][1] /
        // distance[ib][id],dr[grid_index*bxyz_now*max_size_now*3+ib][id][2] / distance[ib][id],ylma);
        dr[pos_tmp + 0] /= distance[dist_tmp];
        dr[pos_tmp + 1] /= distance[dist_tmp];
        dr[pos_tmp + 2] /= distance[dist_tmp];

        // begin calculation

        int index_ylma = dist_tmp * nnnmax_g[0];

        /***************************
        L = 0
        ***************************/
        ylma[index_ylma + 0] = ylmcoef[0]; // l=0, m=0
        if (ucell_atom_nwl[it[dist_tmp]] == 0)
            goto YLM_END;

        /***************************
        L = 1
        ***************************/
        ylma[index_ylma + 1] = ylmcoef[1] * dr[pos_tmp + 2];  // l=1, m=0
        ylma[index_ylma + 2] = -ylmcoef[1] * dr[pos_tmp + 0]; // l=1, m=1
        ylma[index_ylma + 3] = -ylmcoef[1] * dr[pos_tmp + 1]; // l=1, m=-1
        if (ucell_atom_nwl[it[dist_tmp]] == 1)
            goto YLM_END;

        /***************************
        L = 2
        ***************************/
        ylma[index_ylma + 4] = ylmcoef[2] * dr[pos_tmp + 2] * ylma[index_ylma + 1] - ylmcoef[3] * ylma[index_ylma + 0]; // l=2, m=0
        {
            double tmp0 = ylmcoef[4] * dr[pos_tmp + 2];
            ylma[index_ylma + 5] = tmp0 * ylma[index_ylma + 2]; // l=2,m=1
            ylma[index_ylma + 6] = tmp0 * ylma[index_ylma + 3]; // l=2,m=-1

            tmp0 = ylmcoef[4] * dr[pos_tmp + 0];
            ylma[index_ylma + 7] = ylmcoef[5] * ylma[index_ylma + 4] - ylmcoef[6] * ylma[index_ylma + 0] - tmp0 * ylma[index_ylma + 2]; // l=2,m=2
            ylma[index_ylma + 8] = -tmp0 * ylma[index_ylma + 3];
            //	ylma[grid_index*nnnmax+8] = tmp1+tmp2*ylma[grid_index*nnnmax+3];//l=2,m=-2
            if (ucell_atom_nwl[it[dist_tmp]] == 2)
                goto YLM_END;

            /***************************
            L = 3
            ***************************/
            ylma[index_ylma + 9] = ylmcoef[7] * dr[pos_tmp + 2] * ylma[index_ylma + 4] - ylmcoef[8] * ylma[index_ylma + 1]; // l=3, m=0

            tmp0 = ylmcoef[9] * dr[pos_tmp + 2];
            ylma[index_ylma + 10] = tmp0 * ylma[index_ylma + 5] - ylmcoef[10] * ylma[index_ylma + 2]; // l=3,m=1
            ylma[index_ylma + 11] = tmp0 * ylma[index_ylma + 6] - ylmcoef[10] * ylma[index_ylma + 3]; // l=3,m=-1

            tmp0 = ylmcoef[11] * dr[pos_tmp + 2];
            ylma[index_ylma + 12] = tmp0 * ylma[index_ylma + 7]; // l=3,m=2
            ylma[index_ylma + 13] = tmp0 * ylma[index_ylma + 8]; // l=3,m=-2

            tmp0 = ylmcoef[14] * dr[pos_tmp + 0];
            ylma[index_ylma + 14] = ylmcoef[12] * ylma[index_ylma + 10] - ylmcoef[13] * ylma[index_ylma + 2] - tmp0 * ylma[index_ylma + 7]; // l=3,m=3
            ylma[index_ylma + 15] = ylmcoef[12] * ylma[index_ylma + 11] - ylmcoef[13] * ylma[index_ylma + 3] - tmp0 * ylma[index_ylma + 8]; // l=3,m=-3
            if (ucell_atom_nwl[it[dist_tmp]] == 3)
                goto YLM_END;

            /***************************
            L = 4
            ***************************/
            ylma[index_ylma + 16] = ylmcoef[15] * dr[pos_tmp + 2] * ylma[index_ylma + 9] - ylmcoef[16] * ylma[index_ylma + 4]; // l=4,m=0

            tmp0 = ylmcoef[17] * dr[pos_tmp + 2];
            ylma[index_ylma + 17] = tmp0 * ylma[index_ylma + 10] - ylmcoef[18] * ylma[index_ylma + 5]; // l=4,m=1
            ylma[index_ylma + 18] = tmp0 * ylma[index_ylma + 11] - ylmcoef[18] * ylma[index_ylma + 6]; // l=4,m=-1

            tmp0 = ylmcoef[19] * dr[pos_tmp + 2];
            ylma[index_ylma + 19] = tmp0 * ylma[index_ylma + 12] - ylmcoef[20] * ylma[index_ylma + 7]; // l=4,m=2
            ylma[index_ylma + 20] = tmp0 * ylma[index_ylma + 13] - ylmcoef[20] * ylma[index_ylma + 8]; // l=4,m=-2

            tmp0 = 3.0 * dr[pos_tmp + 2];
            ylma[index_ylma + 21] = tmp0 * ylma[index_ylma + 14]; // l=4,m=3
            ylma[index_ylma + 22] = tmp0 * ylma[index_ylma + 15]; // l=4,m=-3

            tmp0 = ylmcoef[23] * dr[pos_tmp + 0];
            ylma[index_ylma + 23] = ylmcoef[21] * ylma[index_ylma + 19] - ylmcoef[22] * ylma[index_ylma + 7] - tmp0 * ylma[index_ylma + 14]; // l=4,m=4
            ylma[index_ylma + 24] = ylmcoef[21] * ylma[index_ylma + 20] - ylmcoef[22] * ylma[index_ylma + 8] - tmp0 * ylma[index_ylma + 15]; // l=4,m=-4
            if (ucell_atom_nwl[it[dist_tmp]] == 4)
                goto YLM_END;

            /***************************
            L = 5
            ***************************/
            ylma[index_ylma + 25] = ylmcoef[24] * dr[pos_tmp + 2] * ylma[index_ylma + 16] - ylmcoef[25] * ylma[index_ylma + 9]; // l=5,m=0

            tmp0 = ylmcoef[26] * dr[pos_tmp + 2];
            ylma[index_ylma + 26] = tmp0 * ylma[index_ylma + 17] - ylmcoef[27] * ylma[index_ylma + 10]; // l=5,m=1
            ylma[index_ylma + 27] = tmp0 * ylma[index_ylma + 18] - ylmcoef[27] * ylma[index_ylma + 11]; // l=5,m=-1

            tmp0 = ylmcoef[28] * dr[pos_tmp + 2];
            ylma[index_ylma + 28] = tmp0 * ylma[index_ylma + 19] - ylmcoef[29] * ylma[index_ylma + 12]; // l=5,m=2
            ylma[index_ylma + 29] = tmp0 * ylma[index_ylma + 20] - ylmcoef[29] * ylma[index_ylma + 13]; // l=5,m=-2

            tmp0 = ylmcoef[30] * dr[pos_tmp + 2];
            ylma[index_ylma + 30] = tmp0 * ylma[index_ylma + 21] - ylmcoef[31] * ylma[index_ylma + 14]; // l=5,m=3
            ylma[index_ylma + 31] = tmp0 * ylma[index_ylma + 22] - ylmcoef[31] * ylma[index_ylma + 15]; // l=5,m=-3

            tmp0 = ylmcoef[32] * dr[pos_tmp + 2];
            ylma[index_ylma + 32] = tmp0 * ylma[index_ylma + 23]; // l=5,m=4
            ylma[index_ylma + 33] = tmp0 * ylma[index_ylma + 24]; // l=5,m=-4

            tmp0 = ylmcoef[35] * dr[pos_tmp + 0];
            ylma[index_ylma + 34] = ylmcoef[33] * ylma[index_ylma + 30] - ylmcoef[34] * ylma[index_ylma + 14] - tmp0 * ylma[index_ylma + 23]; // l=5,m=5
            ylma[index_ylma + 35] = ylmcoef[33] * ylma[index_ylma + 31] - ylmcoef[34] * ylma[index_ylma + 15] - tmp0 * ylma[index_ylma + 24]; // l=5,m=-5
            if (ucell_atom_nwl[it[dist_tmp]] == 5)
                goto YLM_END;

            // if ucell_atom_nwl[it] > 5
            for (int il = 6; il <= ucell_atom_nwl[it[dist_tmp]]; il++)
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

                    ylma[index_ylma + istart + im] = fac2 / sqrt((double)istart - imm * imm) * (dr[pos_tmp + 2] * ylma[index_ylma + istart1 + im] - sqrt((double)istart1 - imm * imm) / fac4 * ylma[index_ylma + istart2 + im]);
                }

                double bl1 = sqrt(2.0 * il / (2.0 * il + 1.0));
                double bl2 = sqrt((2.0 * il - 2.0) / (2.0 * il - 1.0));
                double bl3 = sqrt(2.0) / fac2;

                ylma[index_ylma + istart + 2 * il - 1] = (bl3 * ylma[index_ylma + istart + 2 * il - 5] - bl2 * ylma[index_ylma + istart2 + 2 * il - 5] - 2.0 * dr[pos_tmp + 0] * ylma[index_ylma + istart1 + 2 * il - 3]) / bl1;
                ylma[index_ylma + istart + 2 * il] = (bl3 * ylma[index_ylma + istart + 2 * il - 4] - bl2 * ylma[index_ylma + istart2 + 2 * il - 4] - 2.0 * dr[pos_tmp + 0] * ylma[index_ylma + istart1 + 2 * il - 2]) / bl1;
            }
        }
    YLM_END:
        distance[dist_tmp] /= delta_r_g[0];

        ip[dist_tmp] = (int)(distance[dist_tmp]);
        dx[dist_tmp] = distance[dist_tmp] - ip[dist_tmp];
        dx2[dist_tmp] = dx[dist_tmp] * dx[dist_tmp];
        dx3[dist_tmp] = dx2[dist_tmp] * dx[dist_tmp];

        c3[dist_tmp] = 3.0 * dx2[dist_tmp] - 2.0 * dx3[dist_tmp];
        c1[dist_tmp] = 1.0 - c3[dist_tmp];
        c2[dist_tmp] = (dx[dist_tmp] - 2.0 * dx2[dist_tmp] + dx3[dist_tmp]) * delta_r_g[0];
        c4[dist_tmp] = (dx3[dist_tmp] - dx2[dist_tmp]) * delta_r_g[0];

        int iw;
        double phi = 0.0;
        for (iw = 0; iw < atom_nw[it[dist_tmp]]; ++iw)
        {
            if (atom_iw2_new[it[dist_tmp] * nwmax_g[0] + iw])
            {
                phi = c1[dist_tmp] * psi_u[it[dist_tmp] * nwmax_g[0] * nr_max + iw * nr_max + ip[dist_tmp]] + c2[dist_tmp] * dpsi_u[it[dist_tmp] * nwmax_g[0] * nr_max + iw * nr_max + ip[dist_tmp]] + c3[dist_tmp] * psi_u[it[dist_tmp] * nwmax_g[0] * nr_max + iw * nr_max + ip[dist_tmp] + 1] + c4[dist_tmp] * dpsi_u[it[dist_tmp] * nwmax_g[0] * nr_max + iw * nr_max + ip[dist_tmp] + 1];
            }
            psir_ylm[dist_tmp * nwmax_g[0] + iw] = phi * ylma[dist_tmp * nnnmax_g[0] + atom_iw2_ylm[it[dist_tmp] * nwmax_g[0] + iw]];
        }
        //}//ib
        //}//id
    } // if size
}
__global__ void cu_gamma_vlocal_step3(int ij_index,
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

__global__ void cu_gamma_vlocal_step4(int i,
                                      int j,
                                      int k,
                                      int nby,
                                      int nbz,
                                      int *how_many_atoms,
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
    int atomnow1 = blockIdx.x * blockDim.x + threadIdx.x;
    int atomnow2 = blockIdx.y * blockDim.y + threadIdx.y;
    int grid_index = i * nby * nbz + j * nbz + k;
    if (atomnow1 >= how_many_atoms[grid_index] || atomnow2 >= how_many_atoms[grid_index])
    {
        return;
    }
    int iat1 = which_atom[bcell_start[grid_index] + atomnow1];
    int iat2 = which_atom[bcell_start[grid_index] + atomnow2];

    if (iat2it[iat1] > iat2it[iat2])
    {
        return;
    }

    int iw1, iw2, lo1, lo2, ib, v4;
    double v2;

    lo1 = trace_lo[itiaiw2iwt[iat2it[iat1] * namax_g[0] + iat2ia[iat1]]];
    lo2 = trace_lo[itiaiw2iwt[iat2it[iat2] * namax_g[0] + iat2ia[iat2]]];

    for (iw1 = 0; iw1 < atom_nw[iat2it[iat1]]; iw1++)
    {
        for (iw2 = 0; iw2 < atom_nw[iat2it[iat2]]; iw2++)
        {
            if ((lo1 + iw1) <= (lo2 + iw2))
            {
                // v4=(lo1+iw1)*lgd+lo2+iw2;
                v2 = GridVlocal[v4];
                for (ib = 0; ib < bxyz_g[0]; ib++)
                {
                    if (cal_flag[k * bxyz_g[0] * max_size_g[0] + ib * max_size_g[0] + atomnow1] && cal_flag[k * bxyz_g[0] * max_size_g[0] + ib * max_size_g[0] + atomnow2])
                    {
                        // v1=psir_ylm[k*bxyz_g[0]*max_size_g[0]*nwmax_g[0]+ib*max_size_g[0]*nwmax_g[0]+atomnow1*nwmax_g[0]+iw1]
                        // * vldr3[k*bxyz_g[0]+ib];
                        v2 += psir_ylm[k * bxyz_g[0] * max_size_g[0] * nwmax_g[0] + ib * max_size_g[0] * nwmax_g[0] + atomnow1 * nwmax_g[0] + iw1] * vldr3[k * bxyz_g[0] + ib] * psir_ylm[k * bxyz_g[0] * max_size_g[0] * nwmax_g[0] + ib * max_size_g[0] * nwmax_g[0] + atomnow2 * nwmax_g[0] + iw2];
                    }
                }
                GridVlocal[v4] = v2;
            }
        }
    }
}

__global__ void cu_gamma_vlocal_step4w(int grid_index,
                                       int *how_many_atoms,
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
    int k = blockIdx.x;
    int atom1 = threadIdx.x;
    int atom2 = threadIdx.y;
    if (atom1 >= how_many_atoms[grid_index] || atom2 >= how_many_atoms[grid_index])
    {
        return;
    }
    int iat1 = which_atom[bcell_start[grid_index] + atom1];
    int iat2 = which_atom[bcell_start[grid_index] + atom2];
    int it1 = iat2it[iat1];
    int it2 = iat2it[iat2];
    int lo1 = trace_lo[itiaiw2iwt[it1 * namax_g[0] + iat2ia[iat1]]];
    int lo2 = trace_lo[itiaiw2iwt[it2 * namax_g[0] + iat2ia[iat2]]];
    //if (lo1 <= lo2)
    {
        double v2 = 0.0;
        for (int ib = 0; ib < bxyz_g[0]; ++ib)
        {
            for (int iw1 = 0; iw1 < atom_nw[it1]; iw1++)
            {
                for (int iw2 = 0; iw2 < atom_nw[it2]; iw2++)
                {
                    int lo1_w = lo1 + iw1;
                    int lo2_w = lo2 + iw2;

                    int cal_flag_index = k * bxyz_g[0] * max_size_g[0] + ib * max_size_g[0];
                    int cal_flag_index1 = cal_flag_index + atom1;
                    int cal_flag_index2 = cal_flag_index + atom2;
                    if (cal_flag[cal_flag_index1] && cal_flag[cal_flag_index2])
                    {
                        // v1=psir_ylm[k*bxyz_g[0]*max_size_g[0]*nwmax_g[0]+ib*max_size_g[0]*nwmax_g[0]+atomnow1*nwmax_g[0]+iw1]
                        // * vldr3[k*bxyz_g[0]+ib];
                        v2 += psir_ylm[(cal_flag_index1)*nwmax_g[0] + iw1] * vldr3[k * bxyz_g[0] + ib] * psir_ylm[(cal_flag_index2)*nwmax_g[0] + iw2];
                    }
                    GridVlocal[lo2_w * lgd + lo1_w] += v2;
                }
            }
        }
    }
}

__global__ void cu_gamma_vlocal_step4ww(int gridxy,
                                        int nbz,
                                        int *how_many_atoms,
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

    int lo1 = blockIdx.x * blockDim.x + threadIdx.x;
    int lo2 = blockIdx.y * blockDim.y + threadIdx.y;
    if (lo1 > lo2)
    {
        return;
    }
    if (lo1 >= lgd || lo2 >= lgd)
    {
        return;
    }
    int k, grid_index;
    int size, atom1, atom2, iat1, iat2, it1, it2, iw1, iw2, ib, v2;
    for (k = 0; k < nbz; k++)
    {
        grid_index = gridxy + k;
        size = how_many_atoms[grid_index];
        for (atom1 = 0; atom1 < size; atom1++)
        {
            iat1 = which_atom[bcell_start[grid_index] + atom1];
            it1 = iat2it[iat1];

            for (atom2 = 0; atom2 < size; atom2++)
            {
                iat2 = which_atom[bcell_start[grid_index] + atom2];
                it2 = iat2it[iat2];
                if (it1 <= it2)
                {
                    for (iw1 = 0; iw1 < atom_nw[it1]; iw1++)
                    {
                        if (lo1 == trace_lo[itiaiw2iwt[it1 * namax_g[0] + iat2ia[iat1]]] + iw1)
                        {
                            for (iw2 = 0; iw2 < atom_nw[it2]; iw2++)
                            {
                                if (lo2 == trace_lo[itiaiw2iwt[it2 * namax_g[0] + iat2ia[iat2]]] + iw2)
                                {
                                    v2 = GridVlocal[lo1 * lgd + lo2];
                                    for (ib = 0; ib < bxyz_g[0]; ib++)
                                    {

                                        if (cal_flag[k * bxyz_g[0] * max_size_g[0] + ib * max_size_g[0] + atom1] && cal_flag[k * bxyz_g[0] * max_size_g[0] + ib * max_size_g[0] + atom2])
                                        {
                                            v2 += psir_ylm[k * bxyz_g[0] * max_size_g[0] * nwmax_g[0] + ib * max_size_g[0] * nwmax_g[0] + atom1 * nwmax_g[0] + iw1] * vldr3[k * bxyz_g[0] + ib] * psir_ylm[k * bxyz_g[0] * max_size_g[0] * nwmax_g[0] + ib * max_size_g[0] * nwmax_g[0] + atom2 * nwmax_g[0] + iw2];
                                        }
                                    }
                                    GridVlocal[lo1 * lgd + lo2] = v2;
                                    goto ENDK;
                                }
                            }
                        }
                    }
                }
            }
        }
    ENDK:
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
    const int max_size_now = max_size;
    const int nwmax_now = GlobalC::ucell.nwmax;
    const int namax_now = GlobalC::ucell.namax;
    const int nype_now = GlobalC::ucell.ntype;

    size_t size_meshball_positions = GridT.meshball_ncells * 3;
    size_t size_tau_in_bigcell = GlobalC::ucell.nat * 3;
    size_t size_meshcell_pos = bxyz_now * 3;

    double *meshball_positions_now = new double[size_meshball_positions];
    for (int i = 0; i < GridT.meshball_ncells; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            meshball_positions_now[i * 3 + j] = GridT.meshball_positions[i][j];
        }
    }

    double *tau_in_bigcell_now = new double[size_tau_in_bigcell];
    for (int i = 0; i < GlobalC::ucell.nat; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            tau_in_bigcell_now[i * 3 + j] = GridT.tau_in_bigcell[i][j];
        }
    }

    double *meshcell_pos_now = new double[size_meshcell_pos];
    for (int i = 0; i < bxyz_now; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            meshcell_pos_now[i * 3 + j] = GridT.meshcell_pos[i][j];
        }
    }

    size_t size_phi = GlobalC::ucell.ntype;
    double *phi_getRcut_now = new double[size_phi];
    double max_cut = 0;
    for (int i = 0; i < size_phi; i++)
    {
        phi_getRcut_now[i] = GlobalC::ORB.Phi[i].getRcut();
        if (phi_getRcut_now[i] > max_cut)
        {
            max_cut = phi_getRcut_now[i];
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
    cudaMemcpyToSymbol(max_size_g, &max_size_now, sizeof(int));
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

    double *meshball_positions;
    cudaMalloc((void **)&meshball_positions, size_meshball_positions * sizeof(double));
    cudaMemcpy(meshball_positions,
               meshball_positions_now,
               size_meshball_positions * sizeof(double),
               cudaMemcpyHostToDevice);

    double *tau_in_bigcell;
    cudaMalloc((void **)&tau_in_bigcell, size_tau_in_bigcell * sizeof(double));
    cudaMemcpy(tau_in_bigcell, tau_in_bigcell_now, size_tau_in_bigcell * sizeof(double), cudaMemcpyHostToDevice);

    double *meshcell_pos;
    cudaMalloc((void **)&meshcell_pos, size_meshcell_pos * sizeof(double));
    cudaMemcpy(meshcell_pos, meshcell_pos_now, size_meshcell_pos * sizeof(double), cudaMemcpyHostToDevice);

    double *ORB_Phi_getRcut;
    cudaMalloc((void **)&ORB_Phi_getRcut, size_phi * sizeof(double));
    cudaMemcpy(ORB_Phi_getRcut, phi_getRcut_now, size_phi * sizeof(double), cudaMemcpyHostToDevice);

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

    // para
    double *dr;
    cudaMalloc((void **)&dr, nbz * bxyz_now * max_size_now * 3 * sizeof(double));
    cudaMemset(dr, 0, nbz * bxyz_now * max_size_now * 3 * sizeof(double));

    double *distance;
    cudaMalloc((void **)&distance, nbz * bxyz_now * max_size_now * sizeof(double));
    cudaMemset(distance, 0, nbz * bxyz_now * max_size_now * sizeof(double));

    double *vldr3;
    cudaMalloc((void **)&vldr3, nbz * bxyz_now * sizeof(double));
    cudaMemset(vldr3, 0, nbz * bxyz_now * sizeof(double));

    int *it;
    cudaMalloc((void **)&it, nbz * bxyz_now * max_size_now * sizeof(int));
    cudaMemset(it, 0, nbz * bxyz_now * max_size_now * sizeof(int));

    double *ylma;
    cudaMalloc((void **)&ylma, nbz * max_size_now * bxyz_now * nnnmax * sizeof(double));
    cudaMemset(ylma, 0, nbz * max_size_now * bxyz_now * nnnmax * sizeof(double));

    double *psir_ylm;
    cudaMalloc((void **)&psir_ylm, nbz * max_size_now * bxyz_now * nwmax_now * sizeof(double));
    cudaMemset(psir_ylm, 0, nbz * max_size_now * bxyz_now * nwmax_now * sizeof(double));

    int *ip;
    cudaMalloc((void **)&ip, nbz * max_size_now * bxyz_now * sizeof(int));
    cudaMemset(ip, 0, nbz * max_size_now * bxyz_now * sizeof(int));

    double *dx, *dx2, *dx3, *c1, *c2, *c3, *c4;
    cudaMalloc((void **)&dx, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMemset(dx, 0, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMalloc((void **)&dx2, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMemset(dx2, 0, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMalloc((void **)&dx3, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMemset(dx3, 0, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMalloc((void **)&c1, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMemset(c1, 0, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMalloc((void **)&c2, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMemset(c2, 0, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMalloc((void **)&c3, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMemset(c3, 0, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMalloc((void **)&c4, nbz * max_size_now * bxyz_now * sizeof(double));
    cudaMemset(c4, 0, nbz * max_size_now * bxyz_now * sizeof(double));

    bool *cal_flag;
    cudaMalloc((void **)&cal_flag, nbz * bxyz_now * max_size_now * sizeof(bool));
    cudaMemset(cal_flag, 0, nbz * bxyz_now * max_size_now * sizeof(bool));

    double *GridVlocal;
    cudaMalloc((void **)&GridVlocal, lgd_now * lgd_now * sizeof(double));
    cudaMemset(GridVlocal, 0, lgd_now * lgd_now * sizeof(double));

    // begin kernel

    cudaEventRecord(t2);

    printf("maxsize=%d\n", max_size_now);

    for (int i = 0; i < nbx; i++)
    {
        for (int j = 0; j < nby; j++)
        {
            dim3 grid1(nbz);
            dim3 block1(max_size_now, bxyz_now); // how_many_atoms,bxyz
            cu_gamma_vlocal_step1<<<grid1, block1>>>(i * nby * nbz + j * nbz,
                                                     how_many_atoms,
                                                     bcell_start,
                                                     which_bigcell,
                                                     which_atom,
                                                     iat2it,
                                                     it,
                                                     meshball_positions,
                                                     tau_in_bigcell,
                                                     meshcell_pos,
                                                     dr,
                                                     distance,
                                                     ORB_Phi_getRcut,
                                                     cal_flag,
                                                     ucell_atom_nwl,
                                                     ylma,
                                                     ip,
                                                     dx,
                                                     dx2,
                                                     dx3,
                                                     c1,
                                                     c2,
                                                     c3,
                                                     c4,
                                                     atom_iw2_new,
                                                     atom_iw2_ylm,
                                                     atom_nw,
                                                     nr_max,
                                                     psi_u,
                                                     dpsi_u,
                                                     psir_ylm);

            dim3 grid3(nbz);
            dim3 block3(bx, by, bz);
            cu_gamma_vlocal_step3<<<grid3, block3>>>(i * nby * nbz + j * nbz,
                                                     nbx,
                                                     nby,
                                                     nbz,
                                                     nbz_start,
                                                     ncy,
                                                     nczp,
                                                     vlocal_cu,
                                                     start_ind_g,
                                                     vldr3);

            dim3 grid4(nbz);
            dim3 block4(max_size_now, max_size_now);
            cu_gamma_vlocal_step4w<<<grid4, block4>>>(i * nby * nbz + j * nbz,
                                                      how_many_atoms,
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
        } // j
    }     // i
    cudaDeviceSynchronize();
    cudaMemcpy(GridVlocal_now, GridVlocal, lgd_now * lgd_now * sizeof(double), cudaMemcpyDeviceToHost);
    printf("GridVlocal_now[0]: %lf\n", GridVlocal_now[0]);
    cudaEventRecord(t3);
    cudaDeviceSynchronize();
    // free
    cudaFree(dr);
    cudaFree(distance);
    cudaFree(ylma);
    cudaFree(vldr3);
    cudaFree(it);
    cudaFree(psir_ylm);
    cudaFree(ip);
    cudaFree(dx);
    cudaFree(dx2);
    cudaFree(dx3);
    cudaFree(c1);
    cudaFree(c2);
    cudaFree(c3);
    cudaFree(c4);
    cudaFree(cal_flag);

    cudaFree(how_many_atoms);
    cudaFree(bcell_start);
    cudaFree(which_bigcell);
    cudaFree(which_atom);
    cudaFree(iat2it);
    cudaFree(iat2ia);
    cudaFree(meshball_positions);
    cudaFree(tau_in_bigcell);
    cudaFree(meshcell_pos);
    cudaFree(ORB_Phi_getRcut);
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

    delete[] atom_nw_now;
    delete[] itiaiw2iwt_now;
    delete[] meshball_positions_now;
    delete[] tau_in_bigcell_now;
    delete[] meshcell_pos_now;
    delete[] phi_getRcut_now;
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