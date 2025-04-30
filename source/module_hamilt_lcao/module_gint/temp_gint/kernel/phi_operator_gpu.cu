#include "phi_operator_gpu.h"
#include "phi_operator_kernel.cuh"
#include <cuda_runtime.h>

namespace ModuleGint
{
PhiOperatorGpu::PhiOperatorGpu(std::shared_ptr<const GintGpuVars> gint_gpu_vars, cudaStream_t stream)
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
    checkCuda(cudaEventCreateWithFlags(&event_, cudaEventDisableTiming));
}

PhiOperatorGpu::~PhiOperatorGpu()
{
    checkCuda(cudaEventDestroy(event_));
}

void PhiOperatorGpu::set_bgrid_batch(std::shared_ptr<BatchBigGrid> bgrid_batch)
{
    bgrid_batch_ = bgrid_batch;
    // The set_size here is to determine how many bytes to transfer in the subsequent copy_host_to_device
    atoms_num_info_.set_size(bgrid_batch->get_batch_size());
    bgrids_phi_len_.set_size(bgrid_batch->get_batch_size());
    bgrids_phi_start_.set_size(bgrid_batch->get_batch_size());
    atoms_iat_.set_size(bgrid_batch->get_atoms_num());
    atoms_bgrids_rcoords_.set_size(bgrid_batch->get_atoms_num());
    atoms_phi_start_.set_size(bgrid_batch->get_atoms_num());
    mgrids_local_idx_batch_.set_size(bgrid_batch->get_batch_size() * mgrids_num_);
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
    checkCuda(cudaEventSynchronize(event_));
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

    atoms_num_info_.copy_host_to_device_async();
    bgrids_phi_len_.copy_host_to_device_async();
    bgrids_phi_start_.copy_host_to_device_async();
    atoms_iat_.copy_host_to_device_async();
    atoms_bgrids_rcoords_.copy_host_to_device_async();
    atoms_phi_start_.copy_host_to_device_async();
    mgrids_local_idx_batch_.copy_host_to_device_async();
    checkCuda(cudaEventRecord(event_, stream_));
}

void PhiOperatorGpu::set_phi(double* phi_d)
{
    checkCuda(cudaMemsetAsync(phi_d, 0, phi_len_ * sizeof(double), stream_));
    dim3 grid_dim(mgrids_num_, bgrid_batch_->get_batch_size());
    dim3 threads_per_block(64);
    set_phi_kernel<<<grid_dim, threads_per_block, 0, stream_>>>(
        gint_gpu_vars_->nwmax,
        mgrids_num_,
        gint_gpu_vars_->nr_max,
        gint_gpu_vars_->dr_uniform,
        gint_gpu_vars_->ylmcoef_d,
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
}

void PhiOperatorGpu::phi_mul_vldr3(
    const double* vl_d,
    const double dr3,
    const double* phi_d,
    double* result_d)
{
    checkCuda(cudaMemsetAsync(result_d, 0, phi_len_ * sizeof(double), stream_));
    dim3 grid_dim(mgrids_num_, bgrid_batch_->get_batch_size());
    dim3 threads_per_block(64);
    phi_mul_vldr3_kernel<<<grid_dim, threads_per_block, 0, stream_>>>(
        vl_d,
        dr3,
        phi_d,
        mgrids_num_,
        mgrids_local_idx_batch_.get_device_ptr(),
        bgrids_phi_len_.get_device_ptr(),
        bgrids_phi_start_.get_device_ptr(),
        result_d);
}

void PhiOperatorGpu::phi_mul_phi_vldr3(
    const double* phi_d,
    const double* phi_vldr3_d,
    std::shared_ptr<HContainer<double>> hRGint,
    double* hr_d)
{
    // ap_num means number of atom pairs
    int ap_num = 0;
    int max_m = 0;
    int max_n = 0;
    checkCuda(cudaEventSynchronize(event_));
    for (int i = 0; i < bgrid_batch_->get_batch_size(); i++)
    {
        auto bgrid = bgrid_batch_->get_bgrids()[i];
        // the length of phi on a mesh grid
        const int phi_len_mgrid = bgrid->get_phi_len();
        const int pre_atoms = atoms_num_info_.get_host_ptr()[i].y;
        for (int j = 0; j < bgrid->get_atoms_num(); j++)
        {
            auto atom_1 = bgrid->get_atoms()[j];
            const int iat_1 = atom_1->get_iat();
            const auto& r_1 = atom_1->get_R();
            const int m = atom_1->get_nw();
            const int phi_1_offset = atoms_phi_start_.get_host_ptr()[pre_atoms + j];

            for (int k = 0; k < bgrid->get_atoms_num(); k++)
            {
                auto atom_2 = bgrid->get_atoms()[k];
                const int iat_2 = atom_2->get_iat();
                const auto& r_2 = atom_2->get_R();
                const int n = atom_2->get_nw();

                if(iat_1 > iat_2)
                { continue; }
                
                int hr_offset = hRGint->find_matrix_offset(iat_1, iat_2, r_1 - r_2);
                if (hr_offset == -1)
                { continue; }

                const int phi_2_offset = atoms_phi_start_.get_host_ptr()[pre_atoms + k];


                gemm_A_.get_host_ptr()[ap_num] = phi_d + phi_1_offset;
                gemm_B_.get_host_ptr()[ap_num] = phi_vldr3_d + phi_2_offset;
                gemm_C_.get_host_ptr()[ap_num] = hr_d + hr_offset;
                gemm_lda_.get_host_ptr()[ap_num] = phi_len_mgrid;
                gemm_ldb_.get_host_ptr()[ap_num] = phi_len_mgrid;
                gemm_ldc_.get_host_ptr()[ap_num] = n;
                gemm_m_.get_host_ptr()[ap_num] = m;
                gemm_n_.get_host_ptr()[ap_num] = n;
                gemm_k_.get_host_ptr()[ap_num] = bgrid->get_mgrids_num();
                ap_num++;

                max_m = std::max(max_m, m);
                max_n = std::max(max_n, n);
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
    checkCuda(cudaEventRecord(event_, stream_));
    
    gint_gpu_vars_->fastest_matrix_mul(max_m,
                                       max_n,
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
    checkCudaLastError();
}

}