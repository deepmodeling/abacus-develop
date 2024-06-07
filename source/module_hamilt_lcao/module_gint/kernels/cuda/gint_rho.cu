#include "interp.cuh"
#include "gint_rho.cuh"
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

__global__ void psir_dot(const int nbzp,
                         const int bxyz,
                         const int vec_size,
                         double* vec_a_g,
                         double* vec_b_g,
                         double** results_g)
{
    extern __shared__ double s_data[];
    int tid = threadIdx.x;
    int offset = blockIdx.x * bxyz * vec_size + blockIdx.y * vec_size;
    double* vec_a_mcell = vec_a_g + offset;
    double* vec_b_mcell = vec_b_g + offset;

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