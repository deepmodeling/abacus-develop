#include "../snap_psibeta_gpu.h"
#include "snap_psibeta_kernel.cuh"

#include <cstdio>
#include <cuda_runtime.h>
#include <vector>

namespace module_rt
{
namespace gpu
{

//=============================================================================
// CUDA Error Checking Macro
//=============================================================================

#define CUDA_CHECK(call)                                                                                               \
    do                                                                                                                 \
    {                                                                                                                  \
        cudaError_t err = (call);                                                                                      \
        if (err != cudaSuccess)                                                                                        \
        {                                                                                                              \
            fprintf(stderr, "[CUDA] Error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err));              \
            return false;                                                                                              \
        }                                                                                                              \
    } while (0)

//=============================================================================
// Resource Management - Simplified (using constant memory for grids)
//=============================================================================

void initialize_gpu_resources()
{
    // Check for CUDA device
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0)
    {
        fprintf(stderr, "[CUDA] No CUDA devices found or error getting device count!\n");
        return;
    }

    // Copy integration grids (Lebedev and Gauss-Legendre) to constant memory
    copy_grids_to_device();

    // Sync to ensure all initialization is complete
    cudaDeviceSynchronize();
}

void finalize_gpu_resources()
{
    // No explicit cleanup needed - constant memory is managed by CUDA runtime
    // Just reset any CUDA errors that may have accumulated
    cudaGetLastError();
}

//=============================================================================
// Main GPU Interface Function
//=============================================================================

bool snap_psibeta_neighbor_batch_gpu(const LCAO_Orbitals& orb,
                                     const InfoNonlocal& infoNL_,
                                     const int T1,
                                     const ModuleBase::Vector3<double>& R1,
                                     const Atom* atom1,
                                     const std::vector<int>& all_indexes,
                                     const int npol,
                                     const int T0,
                                     const ModuleBase::Vector3<double>& R0,
                                     const ModuleBase::Vector3<double>& A,
                                     std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>& nlm_all,
                                     const bool calc_r)
{

    // Get projector info
    const int nproj = infoNL_.nproj[T0];
    if (nproj == 0)
    {
        return true;
    }

    // Calculate number of orbitals
    int num_orbitals = all_indexes.size() / npol;
    if (num_orbitals == 0)
    {
        return true;
    }

    // Calculate natomwfc (total projector components)
    int natomwfc = 0;
    std::vector<int> proj_m0_offset_h(nproj);
    for (int ip = 0; ip < nproj; ip++)
    {
        proj_m0_offset_h[ip] = natomwfc;
        int L0 = infoNL_.Beta[T0].Proj[ip].getL();

        // Check L0 against MAX_L
        if (L0 > MAX_L)
        {
            fprintf(stderr, "[CUDA] Warning: L0=%d exceeds MAX_L=%d, falling back to CPU\n", L0, MAX_L);
            return false;
        }

        natomwfc += 2 * L0 + 1;
    }

    int nlm_dim = calc_r ? 4 : 1;

    // Resize output
    if (nlm_all.size() != nlm_dim)
    {
        nlm_all.resize(nlm_dim);
    }

    //=========================================================================
    // Prepare Host Data
    //=========================================================================

    // Orbital data
    std::vector<OrbitalData> orbitals_h(num_orbitals);

    // Calculate total psi size and offsets
    std::vector<int> psi_offsets_h(num_orbitals);
    int total_psi_size = 0;

    for (int i = 0; i < num_orbitals; i++)
    {
        int iw1l = i * npol;
        int iw1 = all_indexes[iw1l] / npol;

        int L1 = atom1->iw2l[iw1];
        int m1 = atom1->iw2m[iw1];
        int N1 = atom1->iw2n[iw1];

        // Check L1 against MAX_L
        if (L1 > MAX_L)
        {
            fprintf(stderr, "[CUDA] Warning: L1=%d exceeds MAX_L=%d, falling back to CPU\n", L1, MAX_L);
            return false;
        }

        const auto& phi_ln = orb.Phi[T1].PhiLN(L1, N1);

        orbitals_h[i].L1 = L1;
        orbitals_h[i].m1 = m1;
        orbitals_h[i].N1 = N1;
        orbitals_h[i].iw_index = all_indexes[iw1l];
        orbitals_h[i].psi_offset = total_psi_size;
        orbitals_h[i].psi_mesh = phi_ln.getNr();
        orbitals_h[i].psi_dk = phi_ln.getDk();
        orbitals_h[i].psi_rcut = phi_ln.getRcut();

        psi_offsets_h[i] = total_psi_size;
        total_psi_size += phi_ln.getNr();
    }

    // Copy psi radial data
    std::vector<double> psi_radial_h(total_psi_size);
    for (int i = 0; i < num_orbitals; i++)
    {
        int iw1l = i * npol;
        int iw1 = all_indexes[iw1l] / npol;

        int L1 = atom1->iw2l[iw1];
        int N1 = atom1->iw2n[iw1];

        const auto& phi_ln = orb.Phi[T1].PhiLN(L1, N1);
        int mesh = phi_ln.getNr();
        const double* psi = phi_ln.getPsi();

        for (int j = 0; j < mesh; j++)
        {
            psi_radial_h[orbitals_h[i].psi_offset + j] = psi[j];
        }
    }

    // Projector data
    std::vector<ProjectorData> projectors_h(nproj);
    std::vector<int> beta_offsets_h(nproj);
    int total_beta_size = 0;

    for (int ip = 0; ip < nproj; ip++)
    {
        const auto& proj = infoNL_.Beta[T0].Proj[ip];

        projectors_h[ip].L0 = proj.getL();
        projectors_h[ip].beta_offset = total_beta_size;
        projectors_h[ip].beta_mesh = proj.getNr();
        projectors_h[ip].beta_dk = proj.getDk();
        projectors_h[ip].beta_rcut = proj.getRcut();

        // Get r_min and r_max from radial grid
        const double* radial = proj.getRadial();
        int mesh = proj.getNr();
        projectors_h[ip].r_min = radial[0];
        projectors_h[ip].r_max = radial[mesh - 1];

        beta_offsets_h[ip] = total_beta_size;
        total_beta_size += mesh;
    }

    // Copy beta radial data
    std::vector<double> beta_radial_h(total_beta_size);
    for (int ip = 0; ip < nproj; ip++)
    {
        const auto& proj = infoNL_.Beta[T0].Proj[ip];
        int mesh = proj.getNr();
        const double* beta_r = proj.getBeta_r();

        for (int j = 0; j < mesh; j++)
        {
            beta_radial_h[projectors_h[ip].beta_offset + j] = beta_r[j];
        }
    }

    //=========================================================================
    // Allocate Device Memory
    //=========================================================================

    OrbitalData* orbitals_d = nullptr;
    ProjectorData* projectors_d = nullptr;
    double* psi_radial_d = nullptr;
    double* beta_radial_d = nullptr;
    int* proj_m0_offset_d = nullptr;
    cuDoubleComplex* nlm_out_d = nullptr;

    size_t output_size = num_orbitals * nlm_dim * natomwfc;

    CUDA_CHECK(cudaMalloc(&orbitals_d, num_orbitals * sizeof(OrbitalData)));
    CUDA_CHECK(cudaMalloc(&projectors_d, nproj * sizeof(ProjectorData)));
    CUDA_CHECK(cudaMalloc(&psi_radial_d, total_psi_size * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&beta_radial_d, total_beta_size * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&proj_m0_offset_d, nproj * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&nlm_out_d, output_size * sizeof(cuDoubleComplex)));

    //=========================================================================
    // Copy Data to Device
    //=========================================================================

    CUDA_CHECK(cudaMemcpy(orbitals_d, orbitals_h.data(), num_orbitals * sizeof(OrbitalData), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(projectors_d, projectors_h.data(), nproj * sizeof(ProjectorData), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(psi_radial_d, psi_radial_h.data(), total_psi_size * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(
        cudaMemcpy(beta_radial_d, beta_radial_h.data(), total_beta_size * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(proj_m0_offset_d, proj_m0_offset_h.data(), nproj * sizeof(int), cudaMemcpyHostToDevice));

    // Initialize output to zero
    CUDA_CHECK(cudaMemset(nlm_out_d, 0, output_size * sizeof(cuDoubleComplex)));

    //=========================================================================
    // Launch Kernel
    //=========================================================================

    double3 R1_d3 = make_double3(R1.x, R1.y, R1.z);
    double3 R0_d3 = make_double3(R0.x, R0.y, R0.z);
    double3 A_d3 = make_double3(A.x, A.y, A.z);

    dim3 grid(num_orbitals, nproj, 1);
    dim3 block(BLOCK_SIZE, 1, 1);

    snap_psibeta_neighbor_batch_kernel<<<grid, block>>>(R1_d3,
                                                        R0_d3,
                                                        A_d3,
                                                        orbitals_d,
                                                        projectors_d,
                                                        psi_radial_d,
                                                        beta_radial_d,
                                                        proj_m0_offset_d,
                                                        num_orbitals,
                                                        nproj,
                                                        natomwfc,
                                                        nlm_dim,
                                                        nlm_out_d);

    // Check for kernel errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "[CUDA] Kernel launch error: %s\n", cudaGetErrorString(err));
        // Cleanup and return
        cudaFree(orbitals_d);
        cudaFree(projectors_d);
        cudaFree(psi_radial_d);
        cudaFree(beta_radial_d);
        cudaFree(proj_m0_offset_d);
        cudaFree(nlm_out_d);
        return false;
    }

    CUDA_CHECK(cudaDeviceSynchronize());

    // Check for errors after synchronization
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "[CUDA] Kernel execution error: %s\n", cudaGetErrorString(err));
        cudaFree(orbitals_d);
        cudaFree(projectors_d);
        cudaFree(psi_radial_d);
        cudaFree(beta_radial_d);
        cudaFree(proj_m0_offset_d);
        cudaFree(nlm_out_d);
        return false;
    }

    //=========================================================================
    // Copy Results Back
    //=========================================================================

    std::vector<cuDoubleComplex> nlm_out_h(output_size);
    CUDA_CHECK(cudaMemcpy(nlm_out_h.data(), nlm_out_d, output_size * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));

    // Reconstruct output into nlm_all format
    for (int i = 0; i < num_orbitals; i++)
    {
        int iw_index = orbitals_h[i].iw_index;

        // Create nlm vectors for this orbital
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

        // Insert into output maps
        for (int dir = 0; dir < nlm_dim; dir++)
        {
            nlm_all[dir].insert({iw_index, nlm[dir]});
        }
    }

    //=========================================================================
    // Cleanup
    //=========================================================================

    cudaFree(orbitals_d);
    cudaFree(projectors_d);
    cudaFree(psi_radial_d);
    cudaFree(beta_radial_d);
    cudaFree(proj_m0_offset_d);
    cudaFree(nlm_out_d);
    return true;
}

} // namespace gpu
} // namespace module_rt
