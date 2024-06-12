#include "gint_vl.cuh"
#include "interp.cuh"
#include "cuda_tools.cuh"
#include "sph.cuh"
namespace GintKernel
{

__global__ void get_psi_and_vldr3(double* ylmcoef,
                                  double delta_r_g,
                                  int bxyz_g,
                                  double nwmax_g,
                                  double* psi_input_double,
                                  int* psi_input_int,
                                  int* atom_num_per_bcell,
                                  int max_atom_per_bcell,
                                  int* ucell_atom_nwl,
                                  bool* atom_iw2_new,
                                  int* atom_iw2_ylm,
                                  int* atom_nw,
                                  int nr_max,
                                  double* psi_u,
                                  double* psi,
                                  double* psi_vldr3)
{
    const int size = atom_num_per_bcell[blockIdx.x];
    const int bcell_start = max_atom_per_bcell * blockIdx.x;
    const int end_index = bcell_start + size;
    const int start_index = bcell_start + threadIdx.x + blockDim.x * blockIdx.y;
    for (int index = start_index; index < end_index;
         index += blockDim.x * gridDim.y)
    {
        double dr[3];
        const int index_double = index * 5;
        dr[0] = psi_input_double[index_double];
        dr[1] = psi_input_double[index_double + 1];
        dr[2] = psi_input_double[index_double + 2];
        const double distance = psi_input_double[index_double + 3];
        const double vlbr3_value = psi_input_double[index_double + 4];
        double ylma[49];
        int index_int = index * 2;
        const int it = psi_input_int[index_int];
        int dist_tmp = psi_input_int[index_int + 1];
        const int nwl = ucell_atom_nwl[it];
        spherical_harmonics(dr, nwl, ylma, ylmcoef);

        interpolate(distance,
                    delta_r_g,
                    it,
                    nwmax_g,
                    nr_max,
                    atom_nw,
                    atom_iw2_new,
                    psi_u,
                    ylma,
                    atom_iw2_ylm,
                    psi,
                    dist_tmp,
                    bxyz_g);

        for (int iw = 0; iw < atom_nw[it]; ++iw)
        {
            psi_vldr3[dist_tmp] = psi[dist_tmp] * vlbr3_value;
            dist_tmp += bxyz_g;
        }
    }
}

} // namespace GintKernel