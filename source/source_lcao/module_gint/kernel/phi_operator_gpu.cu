#include "phi_operator_gpu.h"
#include "phi_operator_kernel.cuh"
#include "dgemm_vbatch.h"
#include <cuda_runtime.h>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <cassert>
#include "source_base/module_device/device_check.h"

namespace ModuleGint
{

template<typename Real>
PhiOperatorGpu<Real>::PhiOperatorGpu(std::shared_ptr<const GintGpuVars> gint_gpu_vars, cudaStream_t stream)
:gint_gpu_vars_(gint_gpu_vars), stream_(stream),
mgrids_num_(BatchBigGrid::get_bgrid_info()->get_mgrids_num()),
atoms_num_info_(BatchBigGrid::get_max_batch_size(), stream_, true),
bgrid_phi_len_(BatchBigGrid::get_max_batch_size(), stream_, true),
bgrid_phi_start_(BatchBigGrid::get_max_batch_size(), stream_, true),
atoms_iat_(BatchBigGrid::get_max_atoms_num(), stream_, true),
atoms_bgrids_rcoords_(BatchBigGrid::get_max_atoms_num(), stream_, true),
atom_phi_start_(BatchBigGrid::get_max_atoms_num(), stream_, true),
batch_mgrid_lidx_(BatchBigGrid::get_max_batch_size()
    * BatchBigGrid::get_bgrid_info()->get_mgrids_num(), stream_, true),
gemm_lda_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_ldb_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_ldc_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_A_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_B_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_C_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_alpha_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true)
{
    CHECK_CUDA(cudaEventCreateWithFlags(&event_, cudaEventDisableTiming));
    // nwmax is the largest per-atom orbital count in the cell, so a stride of
    // nwmax + 1 lets every valid (nw1, nw2) flatten to a distinct bucket key
    // with no artificial cap. Allocate the bucketing scratch once here; the hot
    // path only re-zeroes it.
    nw_stride_ = gint_gpu_vars_->nwmax + 1;
    bucket_counts_.assign(nw_stride_ * nw_stride_, 0);
    bucket_base_.assign(nw_stride_ * nw_stride_, 0);
    bucket_cursor_.assign(nw_stride_ * nw_stride_, 0);
    buckets_.reserve(32);
}

template<typename Real>
PhiOperatorGpu<Real>::~PhiOperatorGpu()
{
    CHECK_CUDA(cudaEventDestroy(event_));
}

template<typename Real>
void PhiOperatorGpu<Real>::set_bgrid_batch(std::shared_ptr<BatchBigGrid> bgrid_batch)
{
    bgrid_batch_ = bgrid_batch;
    auto atoms_num_info_h = atoms_num_info_.get_host_ptr();
    auto bgrid_phi_len_h = bgrid_phi_len_.get_host_ptr();
    auto bgrid_phi_start_h = bgrid_phi_start_.get_host_ptr();
    auto atoms_iat_h = atoms_iat_.get_host_ptr();
    auto atom_rcoords_h = atoms_bgrids_rcoords_.get_host_ptr();
    auto atom_phi_start_h = atom_phi_start_.get_host_ptr();
    auto batch_mgrid_lidx_h = batch_mgrid_lidx_.get_host_ptr();
    int i = 0;
    int j = 0;
    int atoms_accum = 0;
    phi_len_ = 0;
    int phi_start = 0;
    std::vector<int> mgrid_lidx;
    CHECK_CUDA(cudaEventSynchronize(event_));
    for (const auto& bgrid : bgrid_batch->get_bgrids())
    {
        atoms_num_info_h[i] = make_int2(bgrid->get_atoms_num(), atoms_accum);
        atoms_accum += bgrid->get_atoms_num();
        bgrid_phi_start_h[i] = phi_start;
        bgrid->set_mgrids_local_idx(mgrid_lidx);
        std::copy(mgrid_lidx.begin(), mgrid_lidx.end(),
            batch_mgrid_lidx_h + i * mgrids_num_);
        int phi_len_bgrid = 0;
        for (const auto& atom : bgrid->get_atoms())
        {
            atoms_iat_h[j] = atom->get_iat();
            Vec3d rcoord = bgrid->get_bgrid_atom_rcoord(atom);
            atom_rcoords_h[j] = make_double3(rcoord.x, rcoord.y, rcoord.z);
            atom_phi_start_h[j] = phi_len_ + phi_len_bgrid;
            phi_len_bgrid += atom->get_nw();
            j++;
        }
        bgrid_phi_len_h[i] = phi_len_bgrid;
        phi_len_ += phi_len_bgrid * bgrid->get_mgrids_num();
        phi_start += phi_len_bgrid * bgrid->get_mgrids_num();
        i++;
    }

    atoms_num_info_.copy_host_to_device_async(bgrid_batch->get_batch_size());
    bgrid_phi_len_.copy_host_to_device_async(bgrid_batch->get_batch_size());
    bgrid_phi_start_.copy_host_to_device_async(bgrid_batch->get_batch_size());
    atoms_iat_.copy_host_to_device_async(bgrid_batch->get_atoms_num());
    atoms_bgrids_rcoords_.copy_host_to_device_async(bgrid_batch->get_atoms_num());
    atom_phi_start_.copy_host_to_device_async(bgrid_batch->get_atoms_num());
    batch_mgrid_lidx_.copy_host_to_device_async(bgrid_batch->get_batch_size() * mgrids_num_);
    CHECK_CUDA(cudaEventRecord(event_, stream_));

    // Pre-enumerate every (ia_1, ia_2) pair so phi_mul_phi / phi_mul_dm can
    // skip this O(sum atoms_i^2) walk on every call. HContainer lookups and
    // per-bgrid filters stay in the hot path since they depend on the
    // specific caller (hRGint vs dm, symmetric vs not).
    pair_cache_.clear();
    for (int i = 0; i < bgrid_batch->get_batch_size(); ++i)
    {
        const auto& bgrid = bgrid_batch->get_bgrids()[i];
        const int pre_atoms = atoms_num_info_h[i].y;
        const int atoms_num = bgrid->get_atoms_num();
        const int phi_len_mgrid = bgrid->get_phi_len();
        const auto& atoms = bgrid->get_atoms();
        for (int ia_1 = 0; ia_1 < atoms_num; ++ia_1)
        {
            const auto& atom_1 = atoms[ia_1];
            for (int ia_2 = 0; ia_2 < atoms_num; ++ia_2)
            {
                const auto& atom_2 = atoms[ia_2];
                PairInfo p;
                p.phi_1_offset  = atom_phi_start_h[pre_atoms + ia_1];
                p.phi_2_offset  = atom_phi_start_h[pre_atoms + ia_2];
                p.phi_len_mgrid = phi_len_mgrid;
                p.iat_1         = atom_1->get_iat();
                p.iat_2         = atom_2->get_iat();
                p.r_diff        = atom_1->get_R() - atom_2->get_R();
                p.nw1           = static_cast<uint16_t>(atom_1->get_nw());
                p.nw2           = static_cast<uint16_t>(atom_2->get_nw());
                // The shape key (nw1 * nw_stride_ + nw2) indexes a dense
                // nw_stride_^2 table in phi_mul_phi / phi_mul_dm. nw_stride_ =
                // nwmax + 1 and nwmax is by construction the largest per-atom nw
                // in the cell, so this can only trip on an upstream
                // inconsistency rather than an undersized cap.
                assert(p.nw1 < nw_stride_ && p.nw2 < nw_stride_);
                p.ia_le         = static_cast<uint8_t>(ia_1 <= ia_2);
                p.is_diag       = static_cast<uint8_t>(ia_1 == ia_2);
                pair_cache_.push_back(p);
            }
        }
    }
    // Sized to match pair_cache_; Pass 1 of phi_mul_phi / phi_mul_dm overwrites
    // every entry before Pass 2 reads it, so no initial value is needed.
    pair_scratch_offset_.resize(pair_cache_.size());
}

template<typename Real>
void PhiOperatorGpu<Real>::set_phi(Real* phi_d) const
{
    dim3 grid_dim(mgrids_num_, bgrid_batch_->get_batch_size());
    dim3 threads_per_block(64);
    set_phi_kernel<Real><<<grid_dim, threads_per_block, 0, stream_>>>(
        gint_gpu_vars_->nwmax,
        mgrids_num_,
        gint_gpu_vars_->nr_max,
        gint_gpu_vars_->dr_uniform,
        gint_gpu_vars_->ucell_atom_nwl_d,
        gint_gpu_vars_->atom_iw2_new_d,
        gint_gpu_vars_->atom_iw2_ylm_d,
        gint_gpu_vars_->atom_nw_d,
        gint_gpu_vars_->iat2it_d,
        gint_gpu_vars_->rcut_d,
        gint_gpu_vars_->psi_u_d,
        gint_gpu_vars_->dpsi_u_d,
        gint_gpu_vars_->mgrids_pos_d,
        atoms_iat_.get_device_ptr(),
        atoms_bgrids_rcoords_.get_device_ptr(),
        atoms_num_info_.get_device_ptr(),
        atom_phi_start_.get_device_ptr(),
        bgrid_phi_len_.get_device_ptr(),
        phi_d);
    CHECK_LAST_CUDA_ERROR("kernel launch");
}

template<typename Real>
void PhiOperatorGpu<Real>::set_phi_dphi(double* phi_d, double* dphi_x_d, double* dphi_y_d, double* dphi_z_d) const
{
    dim3 grid_dim(mgrids_num_, bgrid_batch_->get_batch_size());
    dim3 threads_per_block(64);
    // Dispatch the WantPhi template based on whether phi is requested.
    // Lets the compiler drop the phi[] stores entirely in the dphi-only case
    // (gint_tau) without paying a per-iw `phi != nullptr` branch in the loop.
    auto launch = [&](auto want_phi) {
        constexpr bool WantPhi = decltype(want_phi)::value;
        set_phi_dphi_kernel<WantPhi><<<grid_dim, threads_per_block, 0, stream_>>>(
            gint_gpu_vars_->nwmax,
            mgrids_num_,
            gint_gpu_vars_->nr_max,
            gint_gpu_vars_->dr_uniform,
            gint_gpu_vars_->ucell_atom_nwl_d,
            gint_gpu_vars_->atom_iw2_new_d,
            gint_gpu_vars_->atom_iw2_ylm_d,
            gint_gpu_vars_->atom_iw2_l_d,
            gint_gpu_vars_->atom_nw_d,
            gint_gpu_vars_->iat2it_d,
            gint_gpu_vars_->rcut_d,
            gint_gpu_vars_->psi_u_d,
            gint_gpu_vars_->dpsi_u_d,
            gint_gpu_vars_->mgrids_pos_d,
            atoms_iat_.get_device_ptr(),
            atoms_bgrids_rcoords_.get_device_ptr(),
            atoms_num_info_.get_device_ptr(),
            atom_phi_start_.get_device_ptr(),
            bgrid_phi_len_.get_device_ptr(),
            phi_d,
            dphi_x_d,
            dphi_y_d,
            dphi_z_d);
    };
    if (phi_d != nullptr) {
        launch(std::true_type{});
    } else {
        launch(std::false_type{});
    }
    CHECK_LAST_CUDA_ERROR("kernel launch");
}

template<typename Real>
void PhiOperatorGpu<Real>::set_ddphi(double* ddphi_xx_d, double* ddphi_xy_d, double* ddphi_xz_d,
                               double* ddphi_yy_d, double* ddphi_yz_d, double* ddphi_zz_d) const
{
    // Since the underlying implementation of `set_ddphi` uses `ddphi +=` instead of `ddphi =`,
    // the ddphi array needs to be zeroed out at the beginning of the function.
    CHECK_CUDA(cudaMemsetAsync(ddphi_xx_d, 0, phi_len_ * sizeof(double), stream_));
    CHECK_CUDA(cudaMemsetAsync(ddphi_xy_d, 0, phi_len_ * sizeof(double), stream_));
    CHECK_CUDA(cudaMemsetAsync(ddphi_xz_d, 0, phi_len_ * sizeof(double), stream_));
    CHECK_CUDA(cudaMemsetAsync(ddphi_yy_d, 0, phi_len_ * sizeof(double), stream_));
    CHECK_CUDA(cudaMemsetAsync(ddphi_yz_d, 0, phi_len_ * sizeof(double), stream_));
    CHECK_CUDA(cudaMemsetAsync(ddphi_zz_d, 0, phi_len_ * sizeof(double), stream_));
    dim3 grid_dim(mgrids_num_, bgrid_batch_->get_batch_size());
    dim3 threads_per_block(64);
    set_ddphi_kernel<<<grid_dim, threads_per_block, 0, stream_>>>(
        gint_gpu_vars_->nwmax,
        mgrids_num_,
        gint_gpu_vars_->nr_max,
        gint_gpu_vars_->dr_uniform,
        gint_gpu_vars_->ucell_atom_nwl_d,
        gint_gpu_vars_->atom_iw2_new_d,
        gint_gpu_vars_->atom_iw2_ylm_d,
        gint_gpu_vars_->atom_iw2_l_d,
        gint_gpu_vars_->atom_nw_d,
        gint_gpu_vars_->iat2it_d,
        gint_gpu_vars_->rcut_d,
        gint_gpu_vars_->psi_u_d,
        gint_gpu_vars_->dpsi_u_d,
        gint_gpu_vars_->mgrids_pos_d,
        atoms_iat_.get_device_ptr(),
        atoms_bgrids_rcoords_.get_device_ptr(),
        atoms_num_info_.get_device_ptr(),
        atom_phi_start_.get_device_ptr(),
        bgrid_phi_len_.get_device_ptr(),
        ddphi_xx_d,
        ddphi_xy_d,
        ddphi_xz_d,
        ddphi_yy_d,
        ddphi_yz_d,
        ddphi_zz_d);
    CHECK_LAST_CUDA_ERROR("kernel launch");
}

template<typename Real>
void PhiOperatorGpu<Real>::phi_mul_vldr3(
    const Real* vl_d,
    const Real dr3,
    const Real* phi_d,
    Real* result_d) const
{
    dim3 grid_dim(mgrids_num_, bgrid_batch_->get_batch_size());
    dim3 threads_per_block(64);
    phi_mul_vldr3_kernel<Real><<<grid_dim, threads_per_block, 0, stream_>>>(
        vl_d,
        dr3,
        phi_d,
        mgrids_num_,
        batch_mgrid_lidx_.get_device_ptr(),
        bgrid_phi_len_.get_device_ptr(),
        bgrid_phi_start_.get_device_ptr(),
        result_d);
    CHECK_LAST_CUDA_ERROR("kernel launch");
}

template<typename Real>
void PhiOperatorGpu<Real>::phi_mul_phi(
    const Real* phi_d,
    const Real* phi_vldr3_d,
    HContainer<double>& hRGint,
    double* hr_d) const
{
    // Shape-exact bucketing: group atom pairs by (nw1, nw2). K = mgrids_num_
    // is already batch-wide constant, so (nw1, nw2) fully determines the GEMM
    // shape. Each bucket hands gemm_tn_vbatch scalar (nw1, nw2, mgrids_num_),
    // so the tile ladder picks the tightest tile for every shape and the
    // wrapper sizes the grid exactly -- no cross-shape tile waste, no
    // over-launched blocks.
    //
    // Algorithm: counting-sort-style two-pass over the pre-enumerated
    // pair_cache_ populated in set_bgrid_batch().
    //   Pass 1: HContainer lookup -> stash hr_offset, count items per shape.
    //   Prefix sum: build the list of non-empty buckets + their flat offsets.
    //   Pass 2: scatter A/B/C pointers + lda/ldb/ldc into the flat host arrays
    //           at each bucket's slot, then one H2D copy per array and one
    //           vbatch launch per bucket. (m/n/k arrays are no longer
    //           scattered -- the wrapper fills them on-device from the
    //           scalar bucket shape.)

    auto& counts = bucket_counts_;
    std::fill(counts.begin(), counts.end(), 0);

    // Pass 1: filter + HContainer lookup + per-shape count.
    for (size_t i = 0; i < pair_cache_.size(); ++i)
    {
        const auto& p = pair_cache_[i];
        if (p.iat_1 > p.iat_2) { pair_scratch_offset_[i] = -1; continue; }
        const int hr = hRGint.find_matrix_offset(p.iat_1, p.iat_2, p.r_diff);
        pair_scratch_offset_[i] = hr;
        if (hr == -1) { continue; }
        counts[p.nw1 * nw_stride_ + p.nw2]++;
    }

    // Prefix sum over dense keys -> compact bucket list.
    auto& buckets = buckets_;
    buckets.clear();
    auto& key_to_base = bucket_base_;
    std::fill(key_to_base.begin(), key_to_base.end(), 0);
    int ap_num = 0;
    for (int k = 0; k < nw_stride_ * nw_stride_; ++k)
    {
        if (counts[k] == 0) { continue; }
        buckets.push_back({k, ap_num, counts[k]});
        key_to_base[k] = ap_num;
        ap_num += counts[k];
    }

    auto* h_A   = gemm_A_.get_host_ptr();
    auto* h_B   = gemm_B_.get_host_ptr();
    auto* h_C   = gemm_C_.get_host_ptr();
    auto* h_lda = gemm_lda_.get_host_ptr();
    auto* h_ldb = gemm_ldb_.get_host_ptr();
    auto* h_ldc = gemm_ldc_.get_host_ptr();

    CHECK_CUDA(cudaEventSynchronize(event_));

    // Pass 2: scatter into the flat host arrays at per-bucket cursors.
    auto& cursor = bucket_cursor_;
    std::fill(cursor.begin(), cursor.end(), 0);
    for (size_t i = 0; i < pair_cache_.size(); ++i)
    {
        const int hr = pair_scratch_offset_[i];
        if (hr == -1) { continue; }
        const auto& p = pair_cache_[i];
        const int key = p.nw1 * nw_stride_ + p.nw2;
        const int pos = key_to_base[key] + cursor[key]++;
        h_A[pos]   = phi_d + p.phi_1_offset;
        h_B[pos]   = phi_vldr3_d + p.phi_2_offset;
        h_C[pos]   = hr_d + hr;
        h_lda[pos] = p.phi_len_mgrid;
        h_ldb[pos] = p.phi_len_mgrid;
        h_ldc[pos] = p.nw2;
    }

    gemm_A_.copy_host_to_device_async(ap_num);
    gemm_B_.copy_host_to_device_async(ap_num);
    gemm_C_.copy_host_to_device_async(ap_num);
    gemm_lda_.copy_host_to_device_async(ap_num);
    gemm_ldb_.copy_host_to_device_async(ap_num);
    gemm_ldc_.copy_host_to_device_async(ap_num);
    CHECK_CUDA(cudaEventRecord(event_, stream_));

    for (const auto& b : buckets)
    {
        const int nw1 = b.key / nw_stride_;
        const int nw2 = b.key % nw_stride_;
        gemm_tn_vbatch<Real>(nw1,
                        nw2,
                        mgrids_num_,
                        gemm_A_.get_device_ptr() + b.off,
                        gemm_lda_.get_device_ptr() + b.off,
                        gemm_B_.get_device_ptr() + b.off,
                        gemm_ldb_.get_device_ptr() + b.off,
                        gemm_C_.get_device_ptr() + b.off,
                        gemm_ldc_.get_device_ptr() + b.off,
                        b.cnt,
                        stream_,
                        nullptr);
    }
}

template<typename Real>
void PhiOperatorGpu<Real>::phi_mul_dm(
    const Real* phi_d,
    const Real* dm_d,
    const HContainer<Real>& dm,
    const bool is_symm,
    double* phi_dm_d)
{
    CHECK_CUDA(cudaMemsetAsync(phi_dm_d, 0, phi_len_ * sizeof(double), stream_));

    // Shape-exact bucketing: same structure as phi_mul_phi, but NN-flavored.
    //   M = mgrids_num_ (batch-wide constant), N = nw2, K = nw1.
    // Shape key is still (nw1, nw2); M is absent from the key since it's
    // identical across every pair in the batch. is_symm selects the
    // upper-triangle (ia_1 <= ia_2) subset and fills per-pair alpha.

    auto& counts = bucket_counts_;
    std::fill(counts.begin(), counts.end(), 0);

    // Pass 1: filter + HContainer lookup + per-shape count.
    for (size_t i = 0; i < pair_cache_.size(); ++i)
    {
        const auto& p = pair_cache_[i];
        if (is_symm && !p.ia_le) { pair_scratch_offset_[i] = -1; continue; }
        const int dm_offset = dm.find_matrix_offset(p.iat_1, p.iat_2, p.r_diff);
        pair_scratch_offset_[i] = dm_offset;
        if (dm_offset == -1) { continue; }
        counts[p.nw1 * nw_stride_ + p.nw2]++;
    }

    // Prefix sum over dense keys -> compact bucket list.
    auto& buckets = buckets_;
    buckets.clear();
    auto& key_to_base = bucket_base_;
    std::fill(key_to_base.begin(), key_to_base.end(), 0);
    int ap_num = 0;
    for (int k = 0; k < nw_stride_ * nw_stride_; ++k)
    {
        if (counts[k] == 0) { continue; }
        buckets.push_back({k, ap_num, counts[k]});
        key_to_base[k] = ap_num;
        ap_num += counts[k];
    }

    auto* h_A     = gemm_A_.get_host_ptr();
    auto* h_B     = gemm_B_.get_host_ptr();
    auto* h_C     = gemm_C_.get_host_ptr();
    auto* h_lda   = gemm_lda_.get_host_ptr();
    auto* h_ldb   = gemm_ldb_.get_host_ptr();
    auto* h_ldc   = gemm_ldc_.get_host_ptr();
    auto* h_alpha = gemm_alpha_.get_host_ptr();

    CHECK_CUDA(cudaEventSynchronize(event_));

    // Pass 2: scatter.
    auto& cursor = bucket_cursor_;
    std::fill(cursor.begin(), cursor.end(), 0);
    for (size_t i = 0; i < pair_cache_.size(); ++i)
    {
        const int dm_offset = pair_scratch_offset_[i];
        if (dm_offset == -1) { continue; }
        const auto& p = pair_cache_[i];
        const int key = p.nw1 * nw_stride_ + p.nw2;
        const int pos = key_to_base[key] + cursor[key]++;
        h_A[pos]   = phi_d + p.phi_1_offset;
        h_B[pos]   = dm_d + dm_offset;
        h_C[pos]   = phi_dm_d + p.phi_2_offset;
        h_lda[pos] = p.phi_len_mgrid;
        h_ldb[pos] = p.nw2;
        h_ldc[pos] = p.phi_len_mgrid;
        if (is_symm)
        {
            h_alpha[pos] = p.is_diag ? Real(1.0) : Real(2.0);
        }
    }

    gemm_A_.copy_host_to_device_async(ap_num);
    gemm_B_.copy_host_to_device_async(ap_num);
    gemm_C_.copy_host_to_device_async(ap_num);
    gemm_lda_.copy_host_to_device_async(ap_num);
    gemm_ldb_.copy_host_to_device_async(ap_num);
    gemm_ldc_.copy_host_to_device_async(ap_num);
    if (is_symm)
    {
        // if is_symm == false, gemm_alpha_ is always 1.0 and is skipped on device
        gemm_alpha_.copy_host_to_device_async(ap_num);
    }
    CHECK_CUDA(cudaEventRecord(event_, stream_));

    for (const auto& b : buckets)
    {
        const int nw1 = b.key / nw_stride_;
        const int nw2 = b.key % nw_stride_;
        auto alpha_ptr = is_symm ? (gemm_alpha_.get_device_ptr() + b.off) : nullptr;
        gemm_nn_vbatch<Real>(mgrids_num_,
                        nw2,
                        nw1,
                        gemm_A_.get_device_ptr() + b.off,
                        gemm_lda_.get_device_ptr() + b.off,
                        gemm_B_.get_device_ptr() + b.off,
                        gemm_ldb_.get_device_ptr() + b.off,
                        gemm_C_.get_device_ptr() + b.off,
                        gemm_ldc_.get_device_ptr() + b.off,
                        b.cnt,
                        stream_,
                        alpha_ptr);
    }
}

template<typename Real>
void PhiOperatorGpu<Real>::phi_dot_phi(
    const Real* phi_i_d,
    const double* phi_j_d,
    double* rho_d) const
{
    dim3 grid_dim(mgrids_num_, bgrid_batch_->get_batch_size());
    dim3 threads_per_block(64);
    phi_dot_phi_kernel<Real, double><<<grid_dim, threads_per_block, sizeof(double) * 32, stream_>>>(
        phi_i_d,
        phi_j_d,
        mgrids_num_,
        batch_mgrid_lidx_.get_device_ptr(),
        bgrid_phi_len_.get_device_ptr(),
        bgrid_phi_start_.get_device_ptr(),
        rho_d);
    CHECK_LAST_CUDA_ERROR("kernel launch");
}

template<typename Real>
void PhiOperatorGpu<Real>::phi_dot_dphi(
    const double* phi_d,
    const double* dphi_x_d,
    const double* dphi_y_d,
    const double* dphi_z_d,
    double* fvl_d) const
{
    dim3 grid_dim(bgrid_batch_->get_max_atoms_per_bgrid(),
                  bgrid_batch_->get_batch_size());
    // Kernel reduce is single-warp; blockDim.x MUST stay 32.
    dim3 threads_per_block(32);
    phi_dot_dphi_kernel<<<grid_dim, threads_per_block, 0, stream_>>>(
        phi_d,
        dphi_x_d,
        dphi_y_d,
        dphi_z_d,
        mgrids_num_,
        bgrid_phi_len_.get_device_ptr(),
        atoms_num_info_.get_device_ptr(),
        atom_phi_start_.get_device_ptr(),
        atoms_iat_.get_device_ptr(),
        gint_gpu_vars_->iat2it_d,
        gint_gpu_vars_->atom_nw_d,
        fvl_d);
    CHECK_LAST_CUDA_ERROR("kernel launch");
}

template<typename Real>
void PhiOperatorGpu<Real>::phi_dot_dphi_r(
    const double* phi_d,
    const double* dphi_x_d,
    const double* dphi_y_d,
    const double* dphi_z_d,
    double* svl_d) const
{
    dim3 grid_dim(mgrids_num_,
                  bgrid_batch_->get_batch_size());
    // Kernel reduce is single-warp; blockDim.x MUST stay 32.
    dim3 threads_per_block(32);
    phi_dot_dphi_r_kernel<<<grid_dim, threads_per_block, 0, stream_>>>(
        phi_d,
        dphi_x_d,
        dphi_y_d,
        dphi_z_d,
        mgrids_num_,
        bgrid_phi_len_.get_device_ptr(),
        atoms_num_info_.get_device_ptr(),
        atom_phi_start_.get_device_ptr(),
        atoms_iat_.get_device_ptr(),
        atoms_bgrids_rcoords_.get_device_ptr(),
        gint_gpu_vars_->mgrids_pos_d,
        gint_gpu_vars_->iat2it_d,
        gint_gpu_vars_->atom_nw_d,
        svl_d);
    CHECK_LAST_CUDA_ERROR("kernel launch");
}

// Explicit instantiations
template class PhiOperatorGpu<double>;
template class PhiOperatorGpu<float>;

}