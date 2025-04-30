#pragma once
#include <memory>
#include <cuda_runtime.h>

#include "module_hamilt_lcao/module_gint/temp_gint/batch_biggrid.h"
#include "module_hamilt_lcao/module_gint/kernels/cuda/cuda_tools.cuh"
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

    void set_phi(double* phi_d);

    void phi_mul_vldr3(
        const double* vl_d,
        const double dr3,
        const double* phi_d,
        double* result_d);
    
    void phi_mul_phi_vldr3(
        const double* phi_d,
        const double* phi_vldr3_d,
        std::shared_ptr<HContainer<double>> hRGint,
        double* hr_d);
    

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

    CudaMemWrapper<int> gemm_m_;
    CudaMemWrapper<int> gemm_n_;
    CudaMemWrapper<int> gemm_k_;
    CudaMemWrapper<int> gemm_lda_;
    CudaMemWrapper<int> gemm_ldb_;
    CudaMemWrapper<int> gemm_ldc_;
    CudaMemWrapper<const double*> gemm_A_;
    CudaMemWrapper<const double*> gemm_B_;
    CudaMemWrapper<double*> gemm_C_; 
    CudaMemWrapper<int> gemm_alpha_;
};

}