#pragma once
#include <memory>
#include <cuda_runtime.h>

#include "source_lcao/module_gint/batch_biggrid.h"
#include "gint_gpu_vars.h"
#include "cuda_mem_wrapper.h"

namespace ModuleGint
{

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
    
    // All GEMM accumulators (hr in phi_mul_phi, phi_dm in phi_mul_dm) are
    // double-typed regardless of Real: when Real=float the multiplies stay in
    // fp32 (cheap) but per-block reductions and device-side atomicAdd run in
    // fp64 so the global reductions don't drift.
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

    mutable CudaMemWrapper<int> gemm_m_;
    mutable CudaMemWrapper<int> gemm_n_;
    mutable CudaMemWrapper<int> gemm_k_;
    mutable CudaMemWrapper<int> gemm_lda_;
    mutable CudaMemWrapper<int> gemm_ldb_;
    mutable CudaMemWrapper<int> gemm_ldc_;
    mutable CudaMemWrapper<const Real*> gemm_A_;
    mutable CudaMemWrapper<const Real*> gemm_B_;
    // Single C-pointer buffer: both phi_mul_phi (output hr) and phi_mul_dm
    // (output phi_dm) write into double* accumulators, so a single shared
    // gemm_C_ device buffer can serve both call sites.
    mutable CudaMemWrapper<double*> gemm_C_;
    mutable CudaMemWrapper<Real> gemm_alpha_;
};

}