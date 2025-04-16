#pragma once
#include <memory>
#include <cuda_runtime.h>

#include "temp_gint/batch_biggrid.h"
#include "module_gint/kernels/cuda/cuda_tools.cuh"
#include "gint_gpu_vars.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace ModuleGint
{

class PhiOperatorGpu
{

public:
    PhiOperatorGpu(cudaStream_t stream);

    void set_bgrid_batch(std::shared_ptr<BatchBigGrid> bgrid_batch);
    static void set_gint_gpu_vars(std::shared_ptr<GintGpuVars> gint_gpu_vars)
    {
        gint_gpu_vars_ = gint_gpu_vars;
    };

    void 

private:
    std::shared_ptr<BatchBigGrid> bgrid_batch_;
    static std::shared_ptr<GintGpuVars> gint_gpu_vars_;
    cudaStream_t stream_ = nullptr;

    // The first number in every group of two represents the number of atoms on that bigcell.
    // The second number represents the cumulative number of atoms up to that bigcell.
    thrust::host_vector<int2> atoms_num_info_h_;
    thrust::device_vector<int2> atoms_num_info_d_;

    // the iat of each atom
    thrust::host_vector<int> atoms_iat_h_;
    thrust::device_vector<int> atoms_iat_d_;

    // atoms_bgrids_rcoords_ here represents the relative coordinates from the big grid to the atoms
    thrust::host_vector<double3> atoms_bgrids_rcoords_h_;
    thrust::device_vector<double3> atoms_bgrids_rcoords_d_;
    
};

}