#include "interp.cuh"
#include "gint_rho.cuh"
#include "sph.cuh"

namespace GintKernel
{

__global__ void get_psi(const double* const ylmcoef,
                        double delta_r_g,
                        int bxyz_g,
                        double nwmax_g,
                        const double* const psi_input_double,
                        const int* const psi_input_int,
                        const int* const atom_num_per_bcell,
                        int max_atom_per_bcell,
                        const int* const ucell_atom_nwl,
                        const bool* const atom_iw2_new,
                        const int* const atom_iw2_ylm,
                        const int* const atom_nw,
                        int nr_max,
                        const double* const psi_u,
                        double* psi)
{
    const int size = atom_num_per_bcell[blockIdx.x];
    const int bcell_start = max_atom_per_bcell * blockIdx.x;
    const int end_index = bcell_start + size;
    const int start_index = bcell_start + threadIdx.x + blockDim.x * blockIdx.y;
    for (int index = start_index; index < end_index;
         index += blockDim.x * gridDim.y)
    {
        double dr[3];
        const int index_double = index * 4;
        dr[0] = psi_input_double[index_double];
        dr[1] = psi_input_double[index_double + 1];
        dr[2] = psi_input_double[index_double + 2];
        const double distance = psi_input_double[index_double + 3];
        double ylma[49];
        const int index_int = index * 2;
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
                    1);
    }
}

__global__ void psir_dot(const int bxyz,
                         const int vec_size,
                         const double* __restrict__ vec_a_g,
                         const double* __restrict__  vec_b_g,
                         double** results_g)
{
    extern __shared__ double s_data[];
    const int tid = threadIdx.x;
    const int offset = blockIdx.x * bxyz * vec_size + blockIdx.y * vec_size;
    const double* vec_a_mcell = vec_a_g + offset;
    const double* vec_b_mcell = vec_b_g + offset;

    s_data[tid] = 0.0;

    for(unsigned int k = tid; k < vec_size; k += blockDim.x)
    {
        s_data[tid] += vec_a_mcell[k] * vec_b_mcell[k];
    }

    __syncthreads();

    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s) {
            s_data[tid] += s_data[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        *results_g[blockIdx.x*bxyz + blockIdx.y] = s_data[0];
    }
}
} // namespace GintKernel