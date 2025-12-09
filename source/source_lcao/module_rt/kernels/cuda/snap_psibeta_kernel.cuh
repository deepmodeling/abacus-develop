#ifndef SNAP_PSIBETA_KERNEL_CUH
#define SNAP_PSIBETA_KERNEL_CUH

#include <cuComplex.h>
#include <cuda_runtime.h>

namespace module_rt
{
namespace gpu
{

// Grid constants
constexpr int RADIAL_GRID_NUM = 140;
constexpr int ANGULAR_GRID_NUM = 110;
constexpr int BLOCK_SIZE = 128;                         // >= ANGULAR_GRID_NUM for efficient reduction
constexpr int MAX_L = 4;                                // Maximum angular momentum supported
constexpr int MAX_YLM_SIZE = (MAX_L + 1) * (MAX_L + 1); // = 25

//=============================================================================
// Device Helper Functions
//=============================================================================

/**
 * @brief GPU version of cubic interpolation for radial functions
 */
__device__ __forceinline__ double interpolate_radial_gpu(const double* __restrict__ psi,
                                                         int mesh,
                                                         double inv_dk,
                                                         double distance)
{
    double position = distance * inv_dk;
    int iq = __double2int_rd(position); // floor

    if (iq > mesh - 4)
        return 0.0;
    if (iq < 0)
        return 0.0;

    double x0 = position - static_cast<double>(iq);
    double x1 = 1.0 - x0;
    double x2 = 2.0 - x0;
    double x3 = 3.0 - x0;

    return x1 * x2 * (psi[iq] * x3 + psi[iq + 3] * x0) / 6.0 + x0 * x3 * (psi[iq + 1] * x2 - psi[iq + 2] * x1) / 2.0;
}

/**
 * @brief Compute exp(i * theta) = cos(theta) + i * sin(theta)
 */
__device__ __forceinline__ cuDoubleComplex cu_exp_i(double theta)
{
    double s, c;
    sincos(theta, &s, &c);
    return make_cuDoubleComplex(c, s);
}

/**
 * @brief Complex multiplication
 */
__device__ __forceinline__ cuDoubleComplex cu_mul(cuDoubleComplex a, cuDoubleComplex b)
{
    return make_cuDoubleComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

/**
 * @brief Complex addition
 */
__device__ __forceinline__ cuDoubleComplex cu_add(cuDoubleComplex a, cuDoubleComplex b)
{
    return make_cuDoubleComplex(a.x + b.x, a.y + b.y);
}

/**
 * @brief Complex conjugate
 */
__device__ __forceinline__ cuDoubleComplex cu_conj(cuDoubleComplex a)
{
    return make_cuDoubleComplex(a.x, -a.y);
}

/**
 * @brief Complex * real
 */
__device__ __forceinline__ cuDoubleComplex cu_mul_real(cuDoubleComplex a, double r)
{
    return make_cuDoubleComplex(a.x * r, a.y * r);
}

/**
 * @brief Compute real spherical harmonics Y_lm at given direction (x, y, z)
 *        This is a simplified version for small L values
 *
 * @param L Angular momentum (0 <= L <= MAX_L)
 * @param x, y, z Direction (should be normalized)
 * @param ylm Output array of size (L+1)^2
 */
__device__ void compute_ylm_gpu(int L, double x, double y, double z, double* ylm);

//=============================================================================
// Kernel Input/Output Structures
//=============================================================================

/**
 * @brief Information about a single projector
 */
struct ProjectorData
{
    int L0;           // Angular momentum
    int beta_offset;  // Offset in flattened beta array
    int beta_mesh;    // Number of mesh points
    double beta_dk;   // Grid spacing
    double beta_rcut; // Cutoff radius
    double r_min;     // Minimum radial value
    double r_max;     // Maximum radial value
};

/**
 * @brief Information about a neighbor-orbital pair for atom-level batching
 *        Used to batch ALL neighbors for a center atom in a single kernel launch
 */
struct NeighborOrbitalData
{
    int neighbor_idx; // Which neighbor (ad index) this orbital belongs to
    double3 R1;       // Neighbor atom position (tau1 * lat0)

    // Orbital info
    int L1;          // Angular momentum
    int m1;          // Magnetic quantum number
    int N1;          // Radial quantum number
    int iw_index;    // Original orbital index for output mapping
    int psi_offset;  // Offset in flattened psi array
    int psi_mesh;    // Number of mesh points
    double psi_dk;   // Grid spacing
    double psi_rcut; // Cutoff radius
};

//=============================================================================
// Main Kernels
//=============================================================================

/**
 * @brief Atom-level batch kernel - processes ALL neighbors for a center atom
 *
 * Grid: (total_neighbor_orbitals, nproj, 1)
 * Block: (BLOCK_SIZE, 1, 1)
 *
 * This kernel processes all orbitals from all neighbors in a single launch,
 * reducing kernel launch overhead significantly.
 *
 * @param R0 Center atom position (projector location)
 * @param A Vector potential
 * @param neighbor_orbitals Array of neighbor-orbital pairs [total_neighbor_orbitals]
 * @param projectors Array of projector data [nproj]
 * @param psi_radial Flattened orbital radial functions
 * @param beta_radial Flattened projector radial functions
 * @param proj_m0_offset Offset for each projector's m=0 in output
 * @param total_neighbor_orbitals Total number of (neighbor, orbital) pairs
 * @param nproj Number of projectors
 * @param natomwfc Total projector components (sum of 2*L0+1)
 * @param nlm_dim 1 for no current, 4 for current
 * @param nlm_out Output [total_neighbor_orbitals * nlm_dim * natomwfc]
 */
__global__ void snap_psibeta_atom_batch_kernel(double3 R0,
                                               double3 A,
                                               const NeighborOrbitalData* __restrict__ neighbor_orbitals,
                                               const ProjectorData* __restrict__ projectors,
                                               const double* __restrict__ psi_radial,
                                               const double* __restrict__ beta_radial,
                                               const int* __restrict__ proj_m0_offset,
                                               int total_neighbor_orbitals,
                                               int nproj,
                                               int natomwfc,
                                               int nlm_dim,
                                               cuDoubleComplex* __restrict__ nlm_out);

//=============================================================================
// Host-side Helper Functions
//=============================================================================

/**
 * @brief Copy integration grids (Lebedev and Gauss-Legendre) to constant memory
 *        Should be called once before kernel execution in each calculate_HR call
 */
void copy_grids_to_device();

} // namespace gpu
} // namespace module_rt

#endif // SNAP_PSIBETA_KERNEL_CUH
