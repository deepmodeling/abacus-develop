#include "gint_vl.cuh"
#include "interp.cuh"
#include "module_hamilt_lcao/module_gint/kernels/cuda/cuda_tools.cuh"
#include "sph.cuh"
namespace GintKernel
{

__global__ void get_psi_and_vldr3(double* ylmcoef,
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
                                  int* atom_nw,
                                  int nr_max,
                                  double* psi_u,
                                  double* psir_ylm_left,
                                  double* psir_ylm_right)
{
    int size = num_psir[blockIdx.x];
    int start_index = psi_size_max * blockIdx.x;
    int end_index = start_index + size;
    start_index += threadIdx.x + blockDim.x * blockIdx.y;
    for (int index = start_index; index < end_index;
         index += blockDim.x * gridDim.y)
    {
        double dr[3];
        int index_double = index * 5;
        dr[0] = input_double[index_double];
        dr[1] = input_double[index_double + 1];
        dr[2] = input_double[index_double + 2];
        double distance = input_double[index_double + 3];
        double vlbr3_value = input_double[index_double + 4];
        double ylma[49];
        int index_int = index * 2;
        int it = input_int[index_int];
        int dist_tmp = input_int[index_int + 1];
        int nwl = ucell_atom_nwl[it];
        spherical_harmonics(dr, distance, nwl, ylma, ylmcoef);

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
                    psir_ylm_left,
                    dist_tmp,
                    bxyz_g);

        for (int iw = 0; iw < atom_nw[it]; ++iw)
        {
            psir_ylm_right[dist_tmp] = psir_ylm_left[dist_tmp] * vlbr3_value;
            dist_tmp += bxyz_g;
        }
    }
}

} // namespace GintKernel