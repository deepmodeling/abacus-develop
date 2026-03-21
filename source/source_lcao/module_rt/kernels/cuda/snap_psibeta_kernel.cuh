/**
 * @file snap_psibeta_kernel.cuh
 * @brief CUDA kernel declarations for computing <psi|beta> overlap integrals
 *
 * This file provides GPU-accelerated computation of overlap integrals between
 * atomic orbitals (psi) and non-local projectors (beta) for real-time TDDFT
 * calculations. The implementation uses numerical integration on a combined
 * radial (Gauss-Legendre) and angular (Lebedev-Laikov) grid.
 *
 * Key Features:
 * - Atom-level batching: processes all neighbors for a center atom in single kernel
 * - Templated spherical harmonics for compile-time optimization
 * - Efficient memory access via constant memory for integration grids
 * - Warp-level reduction for high-performance summation
 */

#ifndef SNAP_PSIBETA_KERNEL_CUH
#define SNAP_PSIBETA_KERNEL_CUH

#include "source_base/tool_quit.h"
#include "source_base/kernels/cuda/sph_harm_gpu.cuh"
#include "source_base/module_device/device_check.h"

#include <cstdio>
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <string>

namespace module_rt
{
namespace gpu
{

//=============================================================================
// Configuration Constants
//=============================================================================

/// Number of points in radial Gauss-Legendre grid
constexpr int RADIAL_GRID_NUM = 140;

/// Number of points in angular Lebedev-Laikov grid (110-point rule)
constexpr int ANGULAR_GRID_NUM = 110;

/// Thread block size for kernel execution
constexpr int BLOCK_SIZE = 128;

/// Maximum supported angular momentum quantum number L
constexpr int MAX_L = 4;

/// Size of spherical harmonics array: (MAX_L + 1)^2 = 25
constexpr int MAX_YLM_SIZE = (MAX_L + 1) * (MAX_L + 1);

/// Maximum number of magnetic quantum numbers for a single L: 2*MAX_L + 1 = 9
constexpr int MAX_M0_SIZE = 2 * MAX_L + 1;

//=============================================================================
// Device Helper Functions - Complex Arithmetic
//=============================================================================

/**
 * @brief Compute exp(i * theta) = cos(theta) + i * sin(theta)
 * @param theta Phase angle in radians
 * @return Complex exponential as cuDoubleComplex
 */
__device__ __forceinline__ cuDoubleComplex cu_exp_i(double theta)
{
    double s, c;
    sincos(theta, &s, &c);
    return make_cuDoubleComplex(c, s);
}

/**
 * @brief Complex multiplication: a * b
 */
__device__ __forceinline__ cuDoubleComplex cu_mul(cuDoubleComplex a, cuDoubleComplex b)
{
    return make_cuDoubleComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

/**
 * @brief Complex addition: a + b
 */
__device__ __forceinline__ cuDoubleComplex cu_add(cuDoubleComplex a, cuDoubleComplex b)
{
    return make_cuDoubleComplex(a.x + b.x, a.y + b.y);
}

/**
 * @brief Complex conjugate: conj(a)
 */
__device__ __forceinline__ cuDoubleComplex cu_conj(cuDoubleComplex a)
{
    return make_cuDoubleComplex(a.x, -a.y);
}

/**
 * @brief Complex times real: a * r
 */
__device__ __forceinline__ cuDoubleComplex cu_mul_real(cuDoubleComplex a, double r)
{
    return make_cuDoubleComplex(a.x * r, a.y * r);
}

//=============================================================================
// Device Helper Functions - Radial Interpolation
//=============================================================================

/**
 * @brief Cubic spline interpolation for radial functions
 *
 * Implements cubic polynomial interpolation using 4 consecutive grid points.
 * This is the GPU equivalent of CPU-side PolyInt::Polynomial_Interpolation.
 *
 * @param psi     Radial function values on uniform grid
 * @param mesh    Number of grid points
 * @param inv_dk  Inverse of grid spacing (1/dk)
 * @param distance Radial distance r at which to interpolate
 * @return Interpolated function value
 */
__device__ __forceinline__ double interpolate_radial_gpu(const double* __restrict__ psi,
                                                         int mesh,
                                                         double inv_dk,
                                                         double distance)
{
    double position = distance * inv_dk;
    int iq = __double2int_rd(position); // floor(position)

    // Boundary checks
    if (iq > mesh - 4 || iq < 0)
    {
        return 0.0;
    }

    // Lagrange interpolation weights
    double x0 = position - static_cast<double>(iq);
    double x1 = 1.0 - x0;
    double x2 = 2.0 - x0;
    double x3 = 3.0 - x0;

    // 4-point Lagrange interpolation formula
    return x1 * x2 * (psi[iq] * x3 + psi[iq + 3] * x0) / 6.0 + x0 * x3 * (psi[iq + 1] * x2 - psi[iq + 2] * x1) / 2.0;
}

//=============================================================================
// Device Helper Functions - Spherical Harmonics
//=============================================================================


//=============================================================================
// Data Structures for Kernel Input
//=============================================================================

/**
 * @brief Non-local projector (beta function) information
 *
 * Contains all data needed to evaluate a single projector during integration.
 */
struct ProjectorData
{
    int L0;           ///< Angular momentum quantum number
    int beta_offset;  ///< Offset into flattened beta radial array
    int beta_mesh;    ///< Number of radial mesh points
    double beta_dk;   ///< Radial grid spacing
    double beta_rcut; ///< Cutoff radius for projector
    double r_min;     ///< Minimum radial grid value (integration start)
    double r_max;     ///< Maximum radial grid value (integration end)
};

/**
 * @brief Neighbor atom orbital information for atom-level batching
 *
 * Each structure represents one (neighbor_atom, orbital) pair that contributes
 * to the overlap integral. This enables processing ALL neighbors for a center
 * atom in a single kernel launch, minimizing launch overhead.
 */
struct NeighborOrbitalData
{
    int neighbor_idx; ///< Index of neighbor atom (ad index in adjacency list)
    double3 R1;       ///< Neighbor atom position in Cartesian coordinates (tau * lat0)

    // Orbital information
    int L1;          ///< Angular momentum of orbital
    int m1;          ///< Magnetic quantum number of orbital
    int N1;          ///< Radial quantum number of orbital
    int iw_index;    ///< Global orbital index for output mapping
    int psi_offset;  ///< Offset into flattened psi radial array
    int psi_mesh;    ///< Number of radial mesh points for orbital
    double psi_dk;   ///< Radial grid spacing for orbital
    double psi_rcut; ///< Cutoff radius for orbital
};

//=============================================================================
// Main CUDA Kernel Declaration
//=============================================================================

/**
 * @brief Atom-level batch kernel for <psi|beta> overlap computation
 *
 * This kernel processes ALL neighbor orbitals for a single center atom in one
 * launch, significantly reducing kernel launch overhead. Each thread block
 * handles the integration for one (neighbor_orbital, projector) pair.
 *
 * Grid Configuration:
 *   - gridDim.x = total_neighbor_orbitals (all orbitals from all neighbors)
 *   - gridDim.y = nproj (number of projectors on center atom)
 *
 * Block Configuration:
 *   - blockDim.x = BLOCK_SIZE threads for parallel integration
 *
 * Integration Strategy:
 *   - Angular loop (outer): each thread processes different angular points
 *   - Radial loop (inner): each thread accumulates over all radial points
 *   - Warp shuffle reduction for efficient summation
 *
 * @param R0                    Center atom position (projector location)
 * @param A                     Vector potential for phase factor
 * @param neighbor_orbitals     Array of neighbor-orbital data [total_neighbor_orbitals]
 * @param projectors            Array of projector data [nproj]
 * @param psi_radial            Flattened array of orbital radial functions
 * @param beta_radial           Flattened array of projector radial functions
 * @param proj_m0_offset        Starting index of each projector's m=0 component in output
 * @param total_neighbor_orbitals Total number of (neighbor, orbital) pairs
 * @param nproj                 Number of projectors on center atom
 * @param natomwfc              Total projector components: sum of (2*L0+1) for all projectors
 * @param nlm_dim               Output dimension: 1 for overlap only, 4 for overlap + current
 * @param nlm_out               Output array [total_neighbor_orbitals * nlm_dim * natomwfc]
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
// Host-side Initialization
//=============================================================================

/**
 * @brief Copy integration grids to GPU constant memory
 *
 * Copies the Lebedev-Laikov angular grid (110 points) and Gauss-Legendre
 * radial grid (140 points) to CUDA constant memory for fast access during
 * kernel execution.
 *
 * @note Must be called once before any kernel launches in a calculation session.
 */
void copy_grids_to_device();

} // namespace gpu
} // namespace module_rt

#endif // SNAP_PSIBETA_KERNEL_CUH
