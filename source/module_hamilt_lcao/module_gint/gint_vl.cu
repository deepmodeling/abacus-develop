#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>

#include "gint_tools.h"

__device__ void sph_harm(const int& Lmax, // max momentum of l
                         const double& xdr,
                         const double& ydr,
                         const double& zdr,
                         double* rly,
                         double* ylmcoef)
{
    // begin calculation
    /***************************
             L = 0
    ***************************/
    rly[0] = ylmcoef[0]; // l=0, m=0
    if (Lmax == 0)
        return;

    /***************************
             L = 1
    ***************************/
    rly[1] = ylmcoef[1] * zdr;  // l=1, m=0
    rly[2] = -ylmcoef[1] * xdr; // l=1, m=1
    rly[3] = -ylmcoef[1] * ydr; // l=1, m=-1
    if (Lmax == 1)
        return;

    /***************************
             L = 2
    ***************************/
    rly[4] = ylmcoef[2] * zdr * rly[1] - ylmcoef[3] * rly[0]; // l=2, m=0

    double tmp0 = ylmcoef[4] * zdr;
    rly[5] = tmp0 * rly[2]; // l=2,m=1
    rly[6] = tmp0 * rly[3]; // l=2,m=-1

    double tmp2 = ylmcoef[4] * xdr;
    rly[7] = ylmcoef[5] * rly[4] - ylmcoef[6] * rly[0] - tmp2 * rly[2]; // l=2,m=2
    rly[8] = -tmp2 * rly[3];
    //	rly[8] = tmp1+tmp2*rly[3];//l=2,m=-2
    if (Lmax == 2)
        return;

    /***************************
             L = 3
    ***************************/
    rly[9] = ylmcoef[7] * zdr * rly[4] - ylmcoef[8] * rly[1]; // l=3, m=0

    double tmp3 = ylmcoef[9] * zdr;
    rly[10] = tmp3 * rly[5] - ylmcoef[10] * rly[2]; // l=3,m=1
    rly[11] = tmp3 * rly[6] - ylmcoef[10] * rly[3]; // l=3,m=-1

    double tmp4 = ylmcoef[11] * zdr;
    rly[12] = tmp4 * rly[7]; // l=3,m=2
    rly[13] = tmp4 * rly[8]; // l=3,m=-2

    double tmp5 = ylmcoef[14] * xdr;
    rly[14] = ylmcoef[12] * rly[10] - ylmcoef[13] * rly[2] - tmp5 * rly[7]; // l=3,m=3
    rly[15] = ylmcoef[12] * rly[11] - ylmcoef[13] * rly[3] - tmp5 * rly[8]; // l=3,m=-3
    if (Lmax == 3)
        return;

    /***************************
             L = 4
    ***************************/
    rly[16] = ylmcoef[15] * zdr * rly[9] - ylmcoef[16] * rly[4]; // l=4,m=0

    double tmp6 = ylmcoef[17] * zdr;
    rly[17] = tmp6 * rly[10] - ylmcoef[18] * rly[5]; // l=4,m=1
    rly[18] = tmp6 * rly[11] - ylmcoef[18] * rly[6]; // l=4,m=-1

    double tmp7 = ylmcoef[19] * zdr;
    rly[19] = tmp7 * rly[12] - ylmcoef[20] * rly[7]; // l=4,m=2
    rly[20] = tmp7 * rly[13] - ylmcoef[20] * rly[8]; // l=4,m=-2

    double tmp8 = 3.0 * zdr;
    rly[21] = tmp8 * rly[14]; // l=4,m=3
    rly[22] = tmp8 * rly[15]; // l=4,m=-3

    double tmp9 = ylmcoef[23] * xdr;
    rly[23] = ylmcoef[21] * rly[19] - ylmcoef[22] * rly[7] - tmp9 * rly[14]; // l=4,m=4
    rly[24] = ylmcoef[21] * rly[20] - ylmcoef[22] * rly[8] - tmp9 * rly[15]; // l=4,m=-4
    if (Lmax == 4)
        return;

    /***************************
             L = 5
    ***************************/
    rly[25] = ylmcoef[24] * zdr * rly[16] - ylmcoef[25] * rly[9]; // l=5,m=0

    double tmp10 = ylmcoef[26] * zdr;
    rly[26] = tmp10 * rly[17] - ylmcoef[27] * rly[10]; // l=5,m=1
    rly[27] = tmp10 * rly[18] - ylmcoef[27] * rly[11]; // l=5,m=-1

    double tmp11 = ylmcoef[28] * zdr;
    rly[28] = tmp11 * rly[19] - ylmcoef[29] * rly[12]; // l=5,m=2
    rly[29] = tmp11 * rly[20] - ylmcoef[29] * rly[13]; // l=5,m=-2

    double tmp12 = ylmcoef[30] * zdr;
    rly[30] = tmp12 * rly[21] - ylmcoef[31] * rly[14]; // l=5,m=3
    rly[31] = tmp12 * rly[22] - ylmcoef[31] * rly[15]; // l=5,m=-3

    double tmp13 = ylmcoef[32] * zdr;
    rly[32] = tmp13 * rly[23]; // l=5,m=4
    rly[33] = tmp13 * rly[24]; // l=5,m=-4

    double tmp14 = ylmcoef[35] * xdr;
    rly[34] = ylmcoef[33] * rly[30] - ylmcoef[34] * rly[14] - tmp14 * rly[23]; // l=5,m=5
    rly[35] = ylmcoef[33] * rly[31] - ylmcoef[34] * rly[15] - tmp14 * rly[24]; // l=5,m=-5
    if (Lmax == 5)
        return;

    // if Lmax > 5
    for (int il = 6; il <= Lmax; il++)
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

            rly[istart + im]
                = fac2 / sqrt((double)istart - imm * imm)
                  * (zdr * rly[istart1 + im] - sqrt((double)istart1 - imm * imm) / fac4 * rly[istart2 + im]);
        }

        double bl1 = sqrt(2.0 * il / (2.0 * il + 1.0));
        double bl2 = sqrt((2.0 * il - 2.0) / (2.0 * il - 1.0));
        double bl3 = sqrt(2.0) / fac2;

        rly[istart + 2 * il - 1]
            = (bl3 * rly[istart + 2 * il - 5] - bl2 * rly[istart2 + 2 * il - 5] - 2.0 * xdr * rly[istart1 + 2 * il - 3])
              / bl1;
        rly[istart + 2 * il]
            = (bl3 * rly[istart + 2 * il - 4] - bl2 * rly[istart2 + 2 * il - 4] - 2.0 * xdr * rly[istart1 + 2 * il - 2])
              / bl1;
    }

    return;
}
__global__ void cal_psir_ylm(
    const Grid_Technique& gt,
    const int bxyz,
    const int na_grid,            // number of atoms on this grid
    const int grid_index,         // 1d index of FFT index (i,j,k)
    const double delta_r,         // delta_r of the uniform FFT grid
    const int* const block_index, // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,  // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag,
    double* const* const psir_ylm,
    int* ucell_atoms_nwl,
    int* ucell_atoms_nw,
    bool** ucell_atoms_iw2_new,
    int** atom_iw2_ylm,
    int it,
    int id,
    int imcell,
    double** psi_uniform,
    double** dpsi_uniform,
    double** meshball_positions,
    double** tau_in_bigcell,
    double** meshcell_pos,
    int* which_bigcell,
    double* ylmcoef,
    int iat) // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
{
    double* ylma;
    double** it_psi_uniform;
    cudaMalloc((void**)&it_psi_uniform, ucell_atoms_nw[iat] * sizeof(double*));
    double** it_dpsi_uniform;
    cudaMalloc((void**)&it_dpsi_uniform, ucell_atoms_nw[iat] * sizeof(double*));

    // preprocess index
    for (int iw = 0; iw < ucell_atoms_nw[iat]; ++iw)
    {
        if (ucell_atoms_iw2_new[iat][iw])
        {
            it_psi_uniform[iw] = psi_uniform[iat];
            it_dpsi_uniform[iw] = dpsi_uniform[iat];
        }
    }

    // meshball_positions should be the bigcell position in meshball
    // to the center of meshball.
    // calculated in cartesian coordinates
    // the std::vector from the grid which is now being operated to the atom position.
    // in meshball language, is the std::vector from imcell to the center cel, plus
    // tau_in_bigcell.
    const double mt[3] = {meshball_positions[imcell][0] - tau_in_bigcell[iat][0],
                          meshball_positions[imcell][1] - tau_in_bigcell[iat][1],
                          meshball_positions[imcell][2] - tau_in_bigcell[iat][2]};

    // number of grids in each big cell (bxyz)
    for (int ib = 0; ib < bxyz; ib++)
    {
        double* p = &psir_ylm[ib][block_index[id]];
        if (!cal_flag[ib][id])
        {
            for (int i = 0; i < block_size[id]; i++)
            {
                p[i] = 0;
            }
        }
        else
        {
            // meshcell_pos: z is the fastest
            const double dr[3]
                = {meshcell_pos[ib][0] + mt[0], meshcell_pos[ib][1] + mt[1], meshcell_pos[ib][2] + mt[2]};
            double distance
                = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]); // distance between atom and grid

            if (distance < 1.0E-9)
                distance += 1.0E-9;

            //------------------------------------------------------
            // spherical harmonic functions Ylm
            //------------------------------------------------------
            //	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
            cudaMalloc((void**)&ylma, (ucell_atoms_nwl[it] + 1) * (ucell_atoms_nwl[it] + 1));

            sph_harm(ucell_atoms_nwl[it], dr[0] / distance, dr[1] / distance, dr[2] / distance, ylma, ylmcoef);
            // these parameters are related to interpolation
            // because once the distance from atom to grid point is known,
            // we can obtain the parameters for interpolation and
            // store them first! these operations can save lots of efforts.
            const double position = distance / delta_r;
            const int ip = static_cast<int>(position);
            const double dx = position - ip;
            const double dx2 = dx * dx;
            const double dx3 = dx2 * dx;

            const double c3 = 3.0 * dx2 - 2.0 * dx3;
            const double c1 = 1.0 - c3;
            const double c2 = (dx - 2.0 * dx2 + dx3) * delta_r;
            const double c4 = (dx3 - dx2) * delta_r;

            double phi = 0;
            for (int iw = 0; ucell_atoms_nw[iat]; ++iw)
            {
                if (ucell_atoms_iw2_new[iat][iw])
                {
                    auto psi_uniform = it_psi_uniform[iw];
                    auto dpsi_uniform = it_dpsi_uniform[iw];
                    phi = c1 * psi_uniform[ip] + c2 * dpsi_uniform[ip] // radial wave functions
                          + c3 * psi_uniform[ip + 1] + c4 * dpsi_uniform[ip + 1];
                }
                p[iw] = phi * ylma[(atom_iw2_ylm[iat])[iw]];
            } // end iw
        }     // end distance<=(GlobalC::ORB.Phi[it].getRcut()-1.0e-15)
    }         // end ib
    return;
}

__device__ double* get_vldr3(const double* const vlocal, // vlocal[ir]
                             const int bxyz,
                             const int bx,
                             const int by,
                             const int bz,
                             const int nplane,
                             const int start_ind,
                             const int ncyz,
                             const double dv)
{
    // set the index for obtaining local potentials
    double* vldr3 = new double[bxyz];

    int bindex = 0;

    for (int ii = 0; ii < bx; ii++)
    {
        const int ipart = ii * ncyz;
        for (int jj = 0; jj < by; jj++)
        {
            const int jpart = jj * nplane + ipart;
            for (int kk = 0; kk < bz; kk++)
            {
                vldr3[bindex] = vlocal[start_ind + kk + jpart] * dv;
                ++bindex;
            }
        }
    }
    return vldr3;
}
__device__ void get_meshcell_atom_pos(const int bxyz,
                                      double orb_phi,
                                      int iat,
                                      int id,
                                      int imcell,
                                      double** meshball_positions,
                                      double** tau_in_bigcell,
                                      double** meshcell_pos,
                                      bool**& cal_flag,
                                      double** dr)
{
    const double mt[3] = {meshball_positions[imcell][0] - tau_in_bigcell[iat][0],
                          meshball_positions[imcell][1] - tau_in_bigcell[iat][1],
                          meshball_positions[imcell][2] - tau_in_bigcell[iat][2]};

    for (int ib = 0; ib < bxyz; ib++)
    {
        dr[ib][0] = meshcell_pos[ib][0] + mt[0];
        dr[ib][1] = meshcell_pos[ib][1] + mt[1];
        dr[ib][2] = meshcell_pos[ib][2] + mt[2];
        const double distance = std::sqrt(dr[ib][0] * dr[ib][0] + dr[ib][1] * dr[ib][1]
                                          + dr[ib][2] * dr[ib][2]); // distance between atom and grid

        if (distance > orb_phi - 1.0e-10)
            cal_flag[ib][id] = false;
        else
            cal_flag[ib][id] = true;
    } // end ib
}
// CUDA内核函数，计算球谐函数并相乘
__global__ void gint_vlocal_cuda(int nbxx,
                                 int bxyz,
                                 int bx,
                                 int by,
                                 int bz,
                                 const int ncyz,
                                 int nplane,
                                 int* start_ind,
                                 int* how_many_atoms,
                                 int* which_atom,
                                 int* bcell_start,
                                 int* iat2it,
                                 int* which_bigcell,
                                 double* vlocal,
                                 double* orb_phi_rcut,
                                 double** meshball_positions,
                                 double** tau_in_bigcell,
                                 double** meshcell_pos,
                                 int* ucell_atoms_nwl,
                                 int* ucell_atoms_nw,
                                 bool** ucell_atoms_iw2_new,
                                 int** atom_iw2_ylm,
                                 int* trace_lo,
                                 int* iat2iatstart,
                                 const double dv)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = 0; i < nbxx; i++)
    {
        const int na_grid = how_many_atoms[i];

        double* vldr3 = get_vldr3(vlocal, // vlocal[ir]
                                  bxyz,
                                  bx,
                                  by,
                                  bz,
                                  nplane,
                                  start_ind[i],
                                  ncyz,
                                  dv);
        int* block_iw = new int[na_grid];
        int* block_index = new int[na_grid + 1];
        int* block_size = new int[na_grid];
        block_index[0] = 0;
        bool** cal_flag = new bool*[bxyz];
        ;
        for (int ib = 0; ib < bxyz; ib++)
        {
            cal_flag[ib] = new bool[na_grid];
        }
        double*** dr;

        for (int id = 0; id < na_grid; id++)
        {
            int iat = which_atom[bcell_start[i] + id]; // index of atom
            int it = iat2it[iat];                       // index of atom type
            const int start = iat2iatstart[iat];        // the index of the first wave function for atom (it,ia)
            block_iw[id] = trace_lo[start];
            block_index[id + 1] = block_index[id] + ucell_atoms_nw[it];
            block_size[id] = ucell_atoms_nw[it];

            int imcell = which_bigcell[bcell_start[i] + id];
            get_meshcell_atom_pos(bxyz,
                                  orb_phi_rcut[it],
                                  iat,
                                  id,
                                  imcell,
                                  meshball_positions,
                                  tau_in_bigcell,
                                  meshcell_pos,
                                  cal_flag,
                                  dr[id]);
        }
    }
}

// 主机函数，处理CUDA计算过程
void gint_vlocal_cuda(int nbxx,
                      int bxyz,
                      int bx,
                      int by,
                      int bz,
                      int nplane,
                      int* start_ind,
                      int* how_many_atoms,
                      double* vlocal,
                      const int ncyz,
                      const double dv,
                      const Grid_Technique& gt)
{

    // 分配和复制数据到设备内存
    int *start_ind_device, *how_many_atoms_device, bcell_start_device;
    int* which_atom_device;
    int* iat2it_device;
    int* iat2ia_device;
    int* trace_lo_device;
    int* ucell_atoms_nw_device;
    int* which_bigcell_device;

    double* tau_in_bigcell_device;
    double* meshball_positions_device;
    double* meshcell_pos_device;
    double* ucell_atoms_nwl;
    cudaMalloc((void**)&start_ind_device, nbxx * sizeof(int));
    cudaMalloc((void**)&how_many_atoms_device, nbxx * sizeof(int));

    cudaMemcpy(start_ind_device, start_ind, nbxx * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(how_many_atoms_device, how_many_atoms, nbxx * sizeof(int), cudaMemcpyHostToDevice);

    // 启动CUDA内核
    // gint_vlocal_cuda<<<100, 100>>>(d_results, numPoints, d_thetas, d_phis);

    // 将计算结果从设备内存复制回主机内存
    // cudaMemcpy(results, d_results, numPoints * sizeof(float), cudaMemcpyDeviceToHost);
}
