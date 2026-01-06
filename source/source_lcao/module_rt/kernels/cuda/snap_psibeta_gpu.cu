/**
 * @file snap_psibeta_gpu.cu
 * @brief Host-side GPU interface for <psi|beta> overlap computation
 *
 * This file provides the high-level interface for GPU-accelerated computation
 * of overlap integrals between atomic orbitals (psi) and non-local projectors
 * (beta). It handles:
 * - GPU resource initialization and cleanup
 * - Data marshalling from ABACUS structures to GPU-friendly formats
 * - Kernel launch configuration
 * - Result unpacking back to ABACUS data structures
 */

#include "../snap_psibeta_gpu.h"
#include "snap_psibeta_kernel.cuh"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

#include <cstdio>
#include <cuda_runtime.h>
#include <string>
#include <vector>

namespace module_rt
{
namespace gpu
{

//=============================================================================
// GPU Resource Management
//=============================================================================

/**
 * @brief Initialize GPU resources for snap_psibeta computation
 *
 * Checks for available CUDA devices and copies integration grids
 * (Lebedev-Laikov angular and Gauss-Legendre radial) to constant memory.
 *
 * @note Call this once at the start of a calculation session before any
 *       snap_psibeta_atom_batch_gpu calls.
 */
void initialize_gpu_resources()
{
    // Verify CUDA device availability
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0)
    {
        ModuleBase::WARNING_QUIT("snap_psibeta_gpu", "No CUDA devices found or error getting device count!");
    }

    // Initialize integration grids in constant memory
    copy_grids_to_device();

    // Synchronize to ensure initialization is complete
    cudaDeviceSynchronize();
}

//=============================================================================
// Internal Helper Structures
//=============================================================================

/**
 * @brief Mapping structure for reconstructing output data
 *
 * Associates each orbital in the flattened GPU array with its original
 * neighbor and orbital indices for proper result placement.
 */
struct OrbitalMapping
{
    int neighbor_idx; ///< Index of neighbor atom in adjacency list
    int iw_index;     ///< Global orbital index for output mapping
};

//=============================================================================
// Main GPU Interface Function
//=============================================================================

/**
 * @brief Compute <psi|beta> overlap integrals on GPU
 *
 * This function processes ALL neighbor atoms for a single center atom (where
 * the projectors are located) in a single kernel launch, providing significant
 * performance improvement over per-neighbor processing.
 *
 * Workflow:
 * 1. Collect all (neighbor, orbital) pairs into flattened arrays
 * 2. Prepare projector data for the center atom
 * 3. Transfer data to GPU and launch kernel
 * 4. Retrieve results and reconstruct nlm_tot structure
 *
 * @param orb       LCAO orbital information
 * @param infoNL_   Non-local projector information
 * @param T0        Atom type of center atom (projector location)
 * @param R0        Position of center atom
 * @param A         Vector potential for phase factor
 * @param adjs      Adjacent atom information
 * @param ucell     Unit cell information
 * @param paraV     Parallel orbital distribution
 * @param npol      Number of spin polarizations
 * @param nlm_dim   Output dimension (1 for overlap only, 4 for overlap + current)
 * @param nlm_tot   Output: overlap integrals indexed as [neighbor][direction][orbital]
 */
void snap_psibeta_atom_batch_gpu(
    const LCAO_Orbitals& orb,
    const InfoNonlocal& infoNL_,
    const int T0,
    const ModuleBase::Vector3<double>& R0,
    const ModuleBase::Vector3<double>& A,
    const AdjacentAtomInfo& adjs,
    const UnitCell* ucell,
    const Parallel_Orbitals* paraV,
    const int npol,
    const int nlm_dim,
    std::vector<std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>>& nlm_tot)
{
    ModuleBase::timer::tick("module_rt", "snap_psibeta_gpu");

    //=========================================================================
    // Early exit if no projectors on center atom
    //=========================================================================

    const int nproj = infoNL_.nproj[T0];
    if (nproj == 0)
    {
        ModuleBase::timer::tick("module_rt", "snap_psibeta_gpu");
        return;
    }

    //=========================================================================
    // Compute projector output indices
    //=========================================================================

    int natomwfc = 0; // Total number of projector components
    std::vector<int> proj_m0_offset_h(nproj);

    for (int ip = 0; ip < nproj; ip++)
    {
        proj_m0_offset_h[ip] = natomwfc;
        int L0 = infoNL_.Beta[T0].Proj[ip].getL();

        // Validate angular momentum
        if (L0 > MAX_L)
        {
            ModuleBase::WARNING_QUIT("snap_psibeta_gpu",
                                     "L0=" + std::to_string(L0) + " exceeds MAX_L=" + std::to_string(MAX_L));
        }
        natomwfc += 2 * L0 + 1;
    }

    //=========================================================================
    // Collect all (neighbor, orbital) pairs
    //=========================================================================

    std::vector<NeighborOrbitalData> neighbor_orbitals_h;
    std::vector<double> psi_radial_h;
    std::vector<OrbitalMapping> orbital_mappings;

    for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
    {
        const int T1 = adjs.ntype[ad];
        const int I1 = adjs.natom[ad];
        const int iat1 = ucell->itia2iat(T1, I1);
        const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
        const Atom* atom1 = &ucell->atoms[T1];

        // Get unique orbital indices (union of row and column indices)
        auto all_indexes = paraV->get_indexes_row(iat1);
        auto col_indexes = paraV->get_indexes_col(iat1);
        all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
        std::sort(all_indexes.begin(), all_indexes.end());
        all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());

        // Process each orbital
        for (size_t iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
        {
            const int iw1 = all_indexes[iw1l] / npol;
            const int L1 = atom1->iw2l[iw1];
            const int m1 = atom1->iw2m[iw1];
            const int N1 = atom1->iw2n[iw1];

            // Skip orbitals with angular momentum beyond supported limit
            if (L1 > MAX_L)
            {
                continue;
            }

            // Get orbital radial function (use getPsi(), not getPsi_r())
            const double* phi_psi = orb.Phi[T1].PhiLN(L1, N1).getPsi();
            int mesh = orb.Phi[T1].PhiLN(L1, N1).getNr();
            double dk = orb.Phi[T1].PhiLN(L1, N1).getDk();
            double rcut = orb.Phi[T1].getRcut();

            // Append to flattened psi array
            size_t psi_offset = psi_radial_h.size();
            psi_radial_h.insert(psi_radial_h.end(), phi_psi, phi_psi + mesh);

            // Create neighbor-orbital data
            NeighborOrbitalData norb;
            norb.neighbor_idx = ad;
            norb.R1 = make_double3(tau1.x * ucell->lat0, tau1.y * ucell->lat0, tau1.z * ucell->lat0);
            norb.L1 = L1;
            norb.m1 = m1;
            norb.N1 = N1;
            norb.iw_index = all_indexes[iw1l];
            norb.psi_offset = static_cast<int>(psi_offset);
            norb.psi_mesh = mesh;
            norb.psi_dk = dk;
            norb.psi_rcut = rcut;

            neighbor_orbitals_h.push_back(norb);

            // Track mapping for result reconstruction
            OrbitalMapping mapping;
            mapping.neighbor_idx = ad;
            mapping.iw_index = all_indexes[iw1l];
            orbital_mappings.push_back(mapping);
        }
    }

    int total_neighbor_orbitals = static_cast<int>(neighbor_orbitals_h.size());
    if (total_neighbor_orbitals == 0)
    {
        ModuleBase::timer::tick("module_rt", "snap_psibeta_gpu");
        return;
    }

    //=========================================================================
    // Prepare projector data
    //=========================================================================

    std::vector<ProjectorData> projectors_h(nproj);
    std::vector<double> beta_radial_h;

    for (int ip = 0; ip < nproj; ip++)
    {
        const auto& proj = infoNL_.Beta[T0].Proj[ip];
        int L0 = proj.getL();
        int mesh = proj.getNr();
        double dk = proj.getDk();
        double rcut = proj.getRcut();
        const double* beta_r = proj.getBeta_r();
        const double* radial = proj.getRadial();

        projectors_h[ip].L0 = L0;
        projectors_h[ip].beta_offset = static_cast<int>(beta_radial_h.size());
        projectors_h[ip].beta_mesh = mesh;
        projectors_h[ip].beta_dk = dk;
        projectors_h[ip].beta_rcut = rcut;
        projectors_h[ip].r_min = radial[0];
        projectors_h[ip].r_max = radial[mesh - 1];

        beta_radial_h.insert(beta_radial_h.end(), beta_r, beta_r + mesh);
    }

    //=========================================================================
    // Allocate GPU memory
    //=========================================================================

    NeighborOrbitalData* neighbor_orbitals_d = nullptr;
    ProjectorData* projectors_d = nullptr;
    double* psi_radial_d = nullptr;
    double* beta_radial_d = nullptr;
    int* proj_m0_offset_d = nullptr;
    cuDoubleComplex* nlm_out_d = nullptr;

    size_t output_size = total_neighbor_orbitals * nlm_dim * natomwfc;

    CUDA_CHECK(cudaMalloc(&neighbor_orbitals_d, total_neighbor_orbitals * sizeof(NeighborOrbitalData)));
    CUDA_CHECK(cudaMalloc(&projectors_d, nproj * sizeof(ProjectorData)));
    CUDA_CHECK(cudaMalloc(&psi_radial_d, psi_radial_h.size() * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&beta_radial_d, beta_radial_h.size() * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&proj_m0_offset_d, nproj * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&nlm_out_d, output_size * sizeof(cuDoubleComplex)));

    //=========================================================================
    // Transfer data to GPU
    //=========================================================================

    CUDA_CHECK(cudaMemcpy(neighbor_orbitals_d,
                          neighbor_orbitals_h.data(),
                          total_neighbor_orbitals * sizeof(NeighborOrbitalData),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(projectors_d, projectors_h.data(), nproj * sizeof(ProjectorData), cudaMemcpyHostToDevice));
    CUDA_CHECK(
        cudaMemcpy(psi_radial_d, psi_radial_h.data(), psi_radial_h.size() * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(
        cudaMemcpy(beta_radial_d, beta_radial_h.data(), beta_radial_h.size() * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(proj_m0_offset_d, proj_m0_offset_h.data(), nproj * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(nlm_out_d, 0, output_size * sizeof(cuDoubleComplex)));

    //=========================================================================
    // Launch kernel
    //=========================================================================

    double3 R0_d3 = make_double3(R0.x, R0.y, R0.z);
    double3 A_d3 = make_double3(A.x, A.y, A.z);

    dim3 grid(total_neighbor_orbitals, nproj, 1);
    dim3 block(BLOCK_SIZE, 1, 1);

    snap_psibeta_atom_batch_kernel<<<grid, block>>>(R0_d3,
                                                    A_d3,
                                                    neighbor_orbitals_d,
                                                    projectors_d,
                                                    psi_radial_d,
                                                    beta_radial_d,
                                                    proj_m0_offset_d,
                                                    total_neighbor_orbitals,
                                                    nproj,
                                                    natomwfc,
                                                    nlm_dim,
                                                    nlm_out_d);

    // Check for launch errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        cudaFree(neighbor_orbitals_d);
        cudaFree(projectors_d);
        cudaFree(psi_radial_d);
        cudaFree(beta_radial_d);
        cudaFree(proj_m0_offset_d);
        cudaFree(nlm_out_d);
        ModuleBase::WARNING_QUIT("snap_psibeta_gpu",
                                 std::string("Atom batch kernel launch error: ") + cudaGetErrorString(err));
    }

    CUDA_CHECK(cudaDeviceSynchronize());

    //=========================================================================
    // Retrieve results
    //=========================================================================

    std::vector<cuDoubleComplex> nlm_out_h(output_size);
    CUDA_CHECK(cudaMemcpy(nlm_out_h.data(), nlm_out_d, output_size * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));

    //=========================================================================
    // Reconstruct output structure
    //=========================================================================

    for (int i = 0; i < total_neighbor_orbitals; i++)
    {
        int ad = orbital_mappings[i].neighbor_idx;
        int iw_index = orbital_mappings[i].iw_index;

        std::vector<std::vector<std::complex<double>>> nlm(nlm_dim);
        for (int d = 0; d < nlm_dim; d++)
        {
            nlm[d].resize(natomwfc);
            for (int k = 0; k < natomwfc; k++)
            {
                size_t idx = i * nlm_dim * natomwfc + d * natomwfc + k;
                nlm[d][k] = std::complex<double>(nlm_out_h[idx].x, nlm_out_h[idx].y);
            }
        }

        // Insert into nlm_tot[neighbor][direction][orbital]
        for (int dir = 0; dir < nlm_dim; dir++)
        {
            nlm_tot[ad][dir].insert({iw_index, nlm[dir]});
        }
    }

    //=========================================================================
    // Cleanup GPU memory
    //=========================================================================

    cudaFree(neighbor_orbitals_d);
    cudaFree(projectors_d);
    cudaFree(psi_radial_d);
    cudaFree(beta_radial_d);
    cudaFree(proj_m0_offset_d);
    cudaFree(nlm_out_d);

    ModuleBase::timer::tick("module_rt", "snap_psibeta_gpu");
}

} // namespace gpu
} // namespace module_rt
