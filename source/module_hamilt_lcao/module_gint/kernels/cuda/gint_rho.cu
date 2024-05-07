#include "interp.cuh"
#include "module_hamilt_lcao/module_gint/kernels/cuda/gint_rho.cuh"
#include "sph.cuh"

namespace GintKernel
{

__global__ void get_psi(const double* const ylmcoef,
                        double delta_r_g,
                        int bxyz_g,
                        double nwmax_g,
                        const double* const input_double,
                        const int* const input_int,
                        const int* const num_psir,
                        int psi_size_max,
                        const int* const ucell_atom_nwl,
                        const bool* const atom_iw2_new,
                        const int* const atom_iw2_ylm,
                        const int* const atom_nw,
                        int nr_max,
                        const double* const psi_u,
                        double* psir_ylm)
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
                    psir_ylm,
                    dist_tmp,
                    1);
    }
}

__global__ void psir_dot(const int* n,
                         double** vec_l_g,
                         int incl,
                         double** vec_r_g,
                         int incr,
                         double** results_g,
                         int batchcount)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = id; i < batchcount; i += stride)
    {
        double* sum = results_g[i];
        double* x = vec_l_g[i];
        double* y = vec_r_g[i];

        for (int j = 0; j < n[i]; j++)
        {
            sum[0] += x[j * incl] * y[j * incr];
        }
    }
}

} // namespace GintKernel