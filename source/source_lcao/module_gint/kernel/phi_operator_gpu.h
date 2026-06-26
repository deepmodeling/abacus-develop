#pragma once
#include <memory>
#include <vector>
#include <cstdint>
#include <cuda_runtime.h>

#include "source_lcao/module_gint/gint_type.h"   // Vec3i, used by PairInfo
#include "source_lcao/module_gint/batch_biggrid.h"
#include "gint_gpu_vars.h"
#include "cuda_mem_wrapper.h"

namespace ModuleGint
{

// Per-atom-pair metadata cached in PhiOperatorGpu::set_bgrid_batch() so that
// phi_mul_phi / phi_mul_dm skip the O(bgrid * atoms^2) enumeration on every
// call. Holds only the fields both callers need; HContainer lookups still
// happen lazily in the hot path because they depend on hRGint / dm.
struct PairInfo
{
    int phi_1_offset;
    int phi_2_offset;
    int phi_len_mgrid;
    int iat_1;
    int iat_2;
    Vec3i r_diff;
    uint16_t nw1;
    uint16_t nw2;
    uint8_t ia_le;   // (ia_1 <= ia_2) within the bgrid, for is_symm filter
    uint8_t is_diag; // (ia_1 == ia_2), for phi_mul_dm is_symm alpha
};

// One non-empty (nw1, nw2) bucket produced by the counting sort in
// phi_mul_phi / phi_mul_dm: a flattened shape key, the bucket's start offset in
// the flat gemm_* arrays, and how many atom pairs landed in it.
struct GemmShapeBucket
{
    int key;  // nw1 * nw_stride_ + nw2
    int off;  // start offset in the gemm_* host/device arrays
    int cnt;  // number of atom pairs in this bucket
};

template<typename Real = double>
class PhiOperatorGpu
{

public:
    PhiOperatorGpu(std::shared_ptr<const GintGpuVars> gint_gpu_vars, cudaStream_t stream = 0);
    ~PhiOperatorGpu();

    void set_bgrid_batch(std::shared_ptr<BatchBigGrid> bgrid_batch);

    void set_phi(Real* phi_d) const;

    // These remain double-only (for force/stress paths)
    void set_phi_dphi(double* phi_d, double* dphi_x_d, double* dphi_y_d, double* dphi_z_d) const;

    void set_ddphi(double* ddphi_xx_d, double* ddphi_xy_d, double* ddphi_xz_d,
                   double* ddphi_yy_d, double* ddphi_yz_d, double* ddphi_zz_d) const;

    void phi_mul_vldr3(
        const Real* vl_d,
        const Real dr3,
        const Real* phi_d,
        Real* result_d) const;
    
    // The GEMM output buffers (hr in phi_mul_phi, phi_dm in phi_mul_dm) are
    // always double, independent of Real. When Real=float the per-pair inner
    // products are reduced in fp32 (cheap); the cross-pair accumulation into a
    // shared hr/phi_dm element is what runs in fp64, via an atomicAdd into
    // these double buffers, so summing many atom-pair contributions doesn't
    // drift.
    void phi_mul_phi(
        const Real* phi_d,
        const Real* phi_vldr3_d,
        HContainer<double>& hRGint,
        double* hr_d) const;

    void phi_mul_dm(
        const Real* phi_d,
        const Real* dm_d,
        const HContainer<Real>& dm,
        const bool is_symm,
        double* phi_dm_d);

    // phi_j_d is the output of phi_mul_dm and therefore always double.
    void phi_dot_phi(
        const Real* phi_i_d,
        const double* phi_j_d,
        double* rho_d) const;
    
    // These remain double-only (for force/stress paths)
    void phi_dot_dphi(
        const double* phi_d,
        const double* dphi_x_d,
        const double* dphi_y_d,
        const double* dphi_z_d,
        double* fvl_d) const;
    
    void phi_dot_dphi_r(
        const double* phi_d,
        const double* dphi_x_d,
        const double* dphi_y_d,
        const double* dphi_z_d,
        double* svl_d) const;

private:
    std::shared_ptr<BatchBigGrid> bgrid_batch_;
    std::shared_ptr<const GintGpuVars> gint_gpu_vars_;

    // the number of meshgrids on a biggrid
    int mgrids_num_;
    
    int phi_len_;

    // Stride for flattening a (nw1, nw2) pair into a single dense bucket key
    // (`nw1 * nw_stride_ + nw2`), so shape-exact bucketing of phi_mul_phi /
    // phi_mul_dm can index a flat table instead of hashing. Set once in the
    // ctor to ucell.nwmax + 1 -- nwmax is the largest per-atom orbital count in
    // the cell, so there is no artificial ceiling: the table is sized to the
    // actual basis (typical nwmax ~25).
    int nw_stride_ = 0;

    cudaStream_t stream_ = 0;
    cudaEvent_t event_;

    // The first number in every group of two represents the number of atoms on that bigcell.
    // The second number represents the cumulative number of atoms up to that bigcell.
    CudaMemWrapper<int2> atoms_num_info_;

    // the iat of each atom
    CudaMemWrapper<int> atoms_iat_;

    // atoms_bgrids_rcoords_ here represents the relative coordinates from the big grid to the atoms
    CudaMemWrapper<double3> atoms_bgrids_rcoords_;

    // the start index of the phi array for each atom
    CudaMemWrapper<int> atom_phi_start_;
    // The length of phi for a single meshgrid on each big grid.
    CudaMemWrapper<int> bgrid_phi_len_;
    // The start index of the phi array for each big grid.
    CudaMemWrapper<int> bgrid_phi_start_;
    // Mapping of the index of meshgrid in the batch of biggrids to the index of meshgrid in the local cell
    CudaMemWrapper<int> batch_mgrid_lidx_;

    mutable CudaMemWrapper<int> gemm_lda_;
    mutable CudaMemWrapper<int> gemm_ldb_;
    mutable CudaMemWrapper<int> gemm_ldc_;
    mutable CudaMemWrapper<const Real*> gemm_A_;
    mutable CudaMemWrapper<const Real*> gemm_B_;
    // C accumulator pointers are always double*: both phi_mul_phi (hr) and
    // phi_mul_dm (phi_dm) write into fp64 buffers via the GEMM's fp64 atomicAdd.
    mutable CudaMemWrapper<double*> gemm_C_;
    mutable CudaMemWrapper<Real> gemm_alpha_;

    // Full (ia_1, ia_2) pair enumeration, rebuilt in set_bgrid_batch().
    // Consumed by phi_mul_phi (TN, iat_1 <= iat_2 filter) and phi_mul_dm
    // (NN, optional is_symm upper-triangle filter).
    std::vector<PairInfo> pair_cache_;

    // Scratch buffer reused across phi_mul_phi / phi_mul_dm calls to cache
    // per-pair HContainer offsets from Pass 1 and replay them in Pass 2
    // without a second find_matrix_offset() call.
    mutable std::vector<int> pair_scratch_offset_;

    // Dense (nw_stride_ * nw_stride_) counting-sort scratch shared by
    // phi_mul_phi / phi_mul_dm. Sized once in the ctor and just re-zeroed per
    // call, so the hot path never reallocates.
    mutable std::vector<int> bucket_counts_;
    mutable std::vector<int> bucket_base_;
    mutable std::vector<int> bucket_cursor_;

    // Compact list of non-empty buckets for the current call. Reused (cleared,
    // not reallocated) by both phi_mul_phi and phi_mul_dm.
    mutable std::vector<GemmShapeBucket> buckets_;
};

}