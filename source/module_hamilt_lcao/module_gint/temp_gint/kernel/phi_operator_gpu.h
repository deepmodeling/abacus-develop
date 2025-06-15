#pragma once
#include <memory>
#include <cuda_runtime.h>

#include "module_hamilt_lcao/module_gint/temp_gint/batch_biggrid.h"
#include "gint_helper.cuh"
#include "gint_gpu_vars.h"
#include "cuda_mem_wrapper.h"

namespace ModuleGint
{

class PhiOperatorGpu
{

public:
    PhiOperatorGpu(std::shared_ptr<const GintGpuVars> gint_gpu_vars, cudaStream_t stream = 0);
    ~PhiOperatorGpu();

    void set_bgrid_batch(std::shared_ptr<BatchBigGrid> bgrid_batch);

    void set_phi(double* phi_d) const;

    void set_phi_dphi(double* phi_d, double* dphi_x_d, double* dphi_y_d, double* dphi_z_d) const;

    void set_ddphi(double* ddphi_xx_d, double* ddphi_xy_d, double* ddphi_xz_d,
                   double* ddphi_yy_d, double* ddphi_yz_d, double* ddphi_zz_d) const;

    void phi_mul_vldr3(
        const double* vl_d,
        const double dr3,
        const double* phi_d,
        double* result_d) const;
    
    void phi_mul_phi(
        const double* phi_d,
        const double* phi_vldr3_d,
        HContainer<double>& hRGint,
        double* hr_d) const;
    
    void phi_mul_dm(
        const double* phi_d,
        const double* dm_d,
        const HContainer<double>& dm,
        const bool is_symm,
        double* phi_dm_d);

    void phi_dot_phi(
        const double* phi_i_d,
        const double* phi_j_d,
        double* rho_d) const;
    
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
    CudaMemWrapper<int> atoms_phi_start_;
    // The length of phi for a single meshgrid on each big grid.
    CudaMemWrapper<int> bgrids_phi_len_;
    // The start index of the phi array for each big grid.
    CudaMemWrapper<int> bgrids_phi_start_;
    // Mapping of the index of meshgrid in the batch of biggrids to the index of meshgrid in the local cell
    CudaMemWrapper<int> mgrids_local_idx_batch_;

    mutable CudaMemWrapper<int> gemm_m_;
    mutable CudaMemWrapper<int> gemm_n_;
    mutable CudaMemWrapper<int> gemm_k_;
    mutable CudaMemWrapper<int> gemm_lda_;
    mutable CudaMemWrapper<int> gemm_ldb_;
    mutable CudaMemWrapper<int> gemm_ldc_;
    mutable CudaMemWrapper<const double*> gemm_A_;
    mutable CudaMemWrapper<const double*> gemm_B_;
    mutable CudaMemWrapper<double*> gemm_C_; 
    mutable CudaMemWrapper<double> gemm_alpha_;
};

}