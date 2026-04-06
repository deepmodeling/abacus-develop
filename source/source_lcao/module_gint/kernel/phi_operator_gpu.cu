#include "phi_operator_gpu.h"
#include "phi_operator_kernel.cuh"
#include "dgemm_vbatch.h"
#include <cuda_runtime.h>
#include "source_base/module_device/device_check.h"

namespace ModuleGint
{

template<typename Real>
PhiOperatorGpu<Real>::PhiOperatorGpu(std::shared_ptr<const GintGpuVars> gint_gpu_vars, cudaStream_t stream)
:gint_gpu_vars_(gint_gpu_vars), stream_(stream),
mgrids_num_(BatchBigGrid::get_bgrid_info()->get_mgrids_num()),
atoms_num_info_(BatchBigGrid::get_max_batch_size(), stream_, true),
bgrids_phi_len_(BatchBigGrid::get_max_batch_size(), stream_, true),
bgrids_phi_start_(BatchBigGrid::get_max_batch_size(), stream_, true),
atoms_iat_(BatchBigGrid::get_max_atoms_num(), stream_, true),
atoms_bgrids_rcoords_(BatchBigGrid::get_max_atoms_num(), stream_, true),
atoms_phi_start_(BatchBigGrid::get_max_atoms_num(), stream_, true),
mgrids_local_idx_batch_(BatchBigGrid::get_max_batch_size() 
    * BatchBigGrid::get_bgrid_info()->get_mgrids_num(), stream_, true),
gemm_m_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_n_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_k_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_lda_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_ldb_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_ldc_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_A_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_B_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_C_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true),
gemm_alpha_(BatchBigGrid::get_max_atom_pairs_num(), stream_, true)
{
    CHECK_CUDA(cudaEventCreateWithFlags(&event_, cudaEventDisableTiming));
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
    auto bgrids_phi_len_h = bgrids_phi_len_.get_host_ptr();
    auto bgrids_phi_start_h = bgrids_phi_start_.get_host_ptr();
    auto atoms_iat_h = atoms_iat_.get_host_ptr();
    auto atoms_bgrids_rcoords_h = atoms_bgrids_rcoords_.get_host_ptr();
    auto atoms_phi_start_h = atoms_phi_start_.get_host_ptr();
    auto mgrids_local_idx_batch_h = mgrids_local_idx_batch_.get_host_ptr();
    int i = 0;
    int j = 0;
    int atoms_accum = 0;
    phi_len_ = 0;
    int phi_start = 0;
    std::vector<int> mgrids_local_idx;
    CHECK_CUDA(cudaEventSynchronize(event_));
    for (const auto& bgrid : bgrid_batch->get_bgrids())
    {
        atoms_num_info_h[i] = make_int2(bgrid->get_atoms_num(), atoms_accum);
        atoms_accum += bgrid->get_atoms_num();
        bgrids_phi_start_h[i] = phi_start;
        bgrid->set_mgrids_local_idx(mgrids_local_idx);
        std::copy(mgrids_local_idx.begin(), mgrids_local_idx.end(),
            mgrids_local_idx_batch_h + i * mgrids_num_);
        int phi_len_bgrid = 0;
        for (const auto& atom : bgrid->get_atoms())
        {
            atoms_iat_h[j] = atom->get_iat();
            Vec3d rcoord = bgrid->get_bgrid_atom_rcoord(atom);
            atoms_bgrids_rcoords_h[j] = make_double3(rcoord.x, rcoord.y, rcoord.z);
            atoms_phi_start_h[j] = phi_len_ + phi_len_bgrid;
            phi_len_bgrid += atom->get_nw();
            j++;
        }
        bgrids_phi_len_h[i] = phi_len_bgrid;
        phi_len_ += phi_len_bgrid * bgrid->get_mgrids_num();
        phi_start += phi_len_bgrid * bgrid->get_mgrids_num();
        i++;
    }

    atoms_num_info_.copy_host_to_device_async(bgrid_batch->get_batch_size());
    bgrids_phi_len_.copy_host_to_device_async(bgrid_batch->get_batch_size());
    bgrids_phi_start_.copy_host_to_device_async(bgrid_batch->get_batch_size());
    atoms_iat_.copy_host_to_device_async(bgrid_batch->get_atoms_num());
    atoms_bgrids_rcoords_.copy_host_to_device_async(bgrid_batch->get_atoms_num());
    atoms_phi_start_.copy_host_to_device_async(bgrid_batch->get_atoms_num());
    mgrids_local_idx_batch_.copy_host_to_device_async(bgrid_batch->get_batch_size() * mgrids_num_);
    CHECK_CUDA(cudaEventRecord(event_, stream_));
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
        atoms_phi_start_.get_device_ptr(),
        bgrids_phi_len_.get_device_ptr(),
        phi_d);
    CHECK_LAST_CUDA_ERROR("kernel launch");
}

template<typename Real>
void PhiOperatorGpu<Real>::set_phi_dphi(double* phi_d, double* dphi_x_d, double* dphi_y_d, double* dphi_z_d) const
{
    dim3 grid_dim(mgrids_num_, bgrid_batch_->get_batch_size());
    dim3 threads_per_block(64);
    set_phi_dphi_kernel<<<grid_dim, threads_per_block, 0, stream_>>>(
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
        atoms_phi_start_.get_device_ptr(),
        bgrids_phi_len_.get_device_ptr(),
        phi_d,
        dphi_x_d,
        dphi_y_d,
        dphi_z_d);
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
        atoms_phi_start_.get_device_ptr(),
        bgrids_phi_len_.get_device_ptr(),
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
        mgrids_local_idx_batch_.get_device_ptr(),
        bgrids_phi_len_.get_device_ptr(),
        bgrids_phi_start_.get_device_ptr(),
        result_d);
    CHECK_LAST_CUDA_ERROR("kernel launch");
}

template<typename Real>
void PhiOperatorGpu<Real>::phi_mul_phi(
    const Real* phi_d,
    const Real* phi_vldr3_d,
    HContainer<Real>& hRGint,
    Real* hr_d) const
{
    // ap_num means number of atom pairs
    int ap_num = 0;
    int max_m = 0;
    int max_n = 0;
    int max_k = mgrids_num_;
    CHECK_CUDA(cudaEventSynchronize(event_));
    for (int i = 0; i < bgrid_batch_->get_batch_size(); i++)
    {
        auto bgrid = bgrid_batch_->get_bgrids()[i];
        // the length of phi on a mesh grid
        const int phi_len_mgrid = bgrid->get_phi_len();
        const int pre_atoms = atoms_num_info_.get_host_ptr()[i].y;
        for (int ia_1 = 0; ia_1 < bgrid->get_atoms_num(); ia_1++)
        {
            auto atom_1 = bgrid->get_atoms()[ia_1];
            const int iat_1 = atom_1->get_iat();
            const auto& r_1 = atom_1->get_R();
            const int nw1 = atom_1->get_nw();
            const int phi_1_offset = atoms_phi_start_.get_host_ptr()[pre_atoms + ia_1];

            for (int ia_2 = 0; ia_2 < bgrid->get_atoms_num(); ia_2++)
            {
                auto atom_2 = bgrid->get_atoms()[ia_2];
                const int iat_2 = atom_2->get_iat();
                const auto& r_2 = atom_2->get_R();
                const int nw2 = atom_2->get_nw();

                if(iat_1 > iat_2)
                { continue; }
                
                int hr_offset = hRGint.find_matrix_offset(iat_1, iat_2, r_1 - r_2);
                if (hr_offset == -1)
                { continue; }

                const int phi_2_offset = atoms_phi_start_.get_host_ptr()[pre_atoms + ia_2];

                gemm_A_.get_host_ptr()[ap_num] = phi_d + phi_1_offset;
                gemm_B_.get_host_ptr()[ap_num] = phi_vldr3_d + phi_2_offset;
                gemm_C_.get_host_ptr()[ap_num] = hr_d + hr_offset;
                gemm_lda_.get_host_ptr()[ap_num] = phi_len_mgrid;
                gemm_ldb_.get_host_ptr()[ap_num] = phi_len_mgrid;
                gemm_ldc_.get_host_ptr()[ap_num] = nw2;
                gemm_m_.get_host_ptr()[ap_num] = nw1;
                gemm_n_.get_host_ptr()[ap_num] = nw2;
                gemm_k_.get_host_ptr()[ap_num] = bgrid->get_mgrids_num();
                ap_num++;

                max_m = std::max(max_m, nw1);
                max_n = std::max(max_n, nw2);
            }
        }
    }

    gemm_A_.copy_host_to_device_async(ap_num);
    gemm_B_.copy_host_to_device_async(ap_num);
    gemm_C_.copy_host_to_device_async(ap_num);
    gemm_lda_.copy_host_to_device_async(ap_num);
    gemm_ldb_.copy_host_to_device_async(ap_num);
    gemm_ldc_.copy_host_to_device_async(ap_num);
    gemm_m_.copy_host_to_device_async(ap_num);
    gemm_n_.copy_host_to_device_async(ap_num);
    gemm_k_.copy_host_to_device_async(ap_num);
    CHECK_CUDA(cudaEventRecord(event_, stream_));
    
    gemm_tn_vbatch<Real>(max_m,
                    max_n,
                    max_k,
                    gemm_m_.get_device_ptr(),
                    gemm_n_.get_device_ptr(),
                    gemm_k_.get_device_ptr(),
                    gemm_A_.get_device_ptr(),
                    gemm_lda_.get_device_ptr(),
                    gemm_B_.get_device_ptr(),
                    gemm_ldb_.get_device_ptr(),
                    gemm_C_.get_device_ptr(),
                    gemm_ldc_.get_device_ptr(),
                    ap_num,
                    stream_,
                    nullptr);
}

template<typename Real>
void PhiOperatorGpu<Real>::phi_mul_dm(
    const Real* phi_d,
    const Real* dm_d,
    const HContainer<Real>& dm,
    const bool is_symm,
    Real* phi_dm_d)
{
    CHECK_CUDA(cudaMemsetAsync(phi_dm_d, 0, phi_len_ * sizeof(Real), stream_));
    // ap_num means number of atom pairs
    int ap_num = 0;
    int max_m = mgrids_num_;
    int max_n = 0;
    int max_k = 0;
    CHECK_CUDA(cudaEventSynchronize(event_));
    for (int i = 0; i < bgrid_batch_->get_batch_size(); i++)
    {
        auto bgrid = bgrid_batch_->get_bgrids()[i];
        // the length of phi on a mesh grid
        const int phi_len_mgrid = bgrid->get_phi_len();
        const int pre_atoms = atoms_num_info_.get_host_ptr()[i].y;
        for (int ia_1 = 0; ia_1 < bgrid->get_atoms_num(); ia_1++)
        {
            auto atom_1 = bgrid->get_atoms()[ia_1];
            const int iat_1 = atom_1->get_iat();
            const auto& r_1 = atom_1->get_R();
            const int nw1 = atom_1->get_nw();
            const int phi_1_offset = atoms_phi_start_.get_host_ptr()[pre_atoms + ia_1];
            int ia_2 = is_symm ? ia_1 : 0;
            for (; ia_2 < bgrid->get_atoms_num(); ia_2++)
            {
                auto atom_2 = bgrid->get_atoms()[ia_2];
                const int iat_2 = atom_2->get_iat();
                const auto& r_2 = atom_2->get_R();
                const int nw2 = atom_2->get_nw();

                int dm_offset = dm.find_matrix_offset(iat_1, iat_2, r_1-r_2);
                if (dm_offset == -1)
                { continue; }

                const int phi_dm_offset = atoms_phi_start_.get_host_ptr()[pre_atoms + ia_2];

                gemm_A_.get_host_ptr()[ap_num] = phi_d + phi_1_offset;
                gemm_B_.get_host_ptr()[ap_num] = dm_d + dm_offset;
                gemm_C_.get_host_ptr()[ap_num] = phi_dm_d + phi_dm_offset;
                gemm_lda_.get_host_ptr()[ap_num] = phi_len_mgrid;
                gemm_ldb_.get_host_ptr()[ap_num] = nw2;
                gemm_ldc_.get_host_ptr()[ap_num] = phi_len_mgrid;
                gemm_m_.get_host_ptr()[ap_num] = mgrids_num_;
                gemm_n_.get_host_ptr()[ap_num] = nw2;
                gemm_k_.get_host_ptr()[ap_num] = nw1;
                gemm_alpha_.get_host_ptr()[ap_num] = ia_1 == ia_2 ? Real(1.0) : Real(2.0);
                ap_num++;

                max_n = std::max(max_n, nw2);
                max_k = std::max(max_k, nw1);
            }
        }
    }

    gemm_A_.copy_host_to_device_async(ap_num);
    gemm_B_.copy_host_to_device_async(ap_num);
    gemm_C_.copy_host_to_device_async(ap_num);
    gemm_lda_.copy_host_to_device_async(ap_num);
    gemm_ldb_.copy_host_to_device_async(ap_num);
    gemm_ldc_.copy_host_to_device_async(ap_num);
    gemm_m_.copy_host_to_device_async(ap_num);
    gemm_n_.copy_host_to_device_async(ap_num);
    gemm_k_.copy_host_to_device_async(ap_num);
    if(is_symm)
    {
        // if is_symm == false, gemm_alpha_ always equals 1.0,
        // so we don't need to copy it to device
        gemm_alpha_.copy_host_to_device_async(ap_num);
    }
    CHECK_CUDA(cudaEventRecord(event_, stream_));

    auto alpha_ptr = is_symm ? gemm_alpha_.get_device_ptr() : nullptr;
    gemm_nn_vbatch<Real>(max_m,
                    max_n,
                    max_k,
                    gemm_m_.get_device_ptr(),
                    gemm_n_.get_device_ptr(),
                    gemm_k_.get_device_ptr(),
                    gemm_A_.get_device_ptr(),
                    gemm_lda_.get_device_ptr(),
                    gemm_B_.get_device_ptr(),
                    gemm_ldb_.get_device_ptr(),
                    gemm_C_.get_device_ptr(),
                    gemm_ldc_.get_device_ptr(),
                    ap_num,
                    stream_,
                    alpha_ptr);
}

template<typename Real>
void PhiOperatorGpu<Real>::phi_dot_phi(
    const Real* phi_i_d,
    const Real* phi_j_d,
    Real* rho_d) const
{
    dim3 grid_dim(mgrids_num_, bgrid_batch_->get_batch_size());
    dim3 threads_per_block(64);
    phi_dot_phi_kernel<Real><<<grid_dim, threads_per_block, sizeof(Real) * 32, stream_>>>(
        phi_i_d,
        phi_j_d,
        mgrids_num_,
        mgrids_local_idx_batch_.get_device_ptr(),
        bgrids_phi_len_.get_device_ptr(),
        bgrids_phi_start_.get_device_ptr(),
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
    dim3 grid_dim(bgrid_batch_->get_max_atoms_num_per_bgrid(),
                  bgrid_batch_->get_batch_size());
    dim3 threads_per_block(32);
    phi_dot_dphi_kernel<<<grid_dim, threads_per_block, sizeof(double) * 32 * 3, stream_>>>(
        phi_d,
        dphi_x_d,
        dphi_y_d,
        dphi_z_d,
        mgrids_num_,
        bgrids_phi_len_.get_device_ptr(),
        atoms_num_info_.get_device_ptr(),
        atoms_phi_start_.get_device_ptr(),
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
    dim3 threads_per_block(32);
    phi_dot_dphi_r_kernel<<<grid_dim, threads_per_block, sizeof(double) * 32 * 6, stream_>>>(
        phi_d,
        dphi_x_d,
        dphi_y_d,
        dphi_z_d,
        mgrids_num_,
        bgrids_phi_len_.get_device_ptr(),
        atoms_num_info_.get_device_ptr(),
        atoms_phi_start_.get_device_ptr(),
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