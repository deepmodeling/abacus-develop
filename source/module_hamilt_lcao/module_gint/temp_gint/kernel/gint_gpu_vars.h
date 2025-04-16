#pragma once

#include <cuda_runtime.h>
#include "set_const_mem.cuh"
#include "module_base/ylm.h"
#include "module_cell/unitcell.h"
#include "module_cell/atom_spec.h"
#include "module_gint/kernels/cuda/gemm_selector.cuh"
#include "module_gint/temp_gint/gint_info.h"
#include "module_gint/kernels/cuda/cuda_tools.cuh"

namespace ModuleGint
{

class GintGpuVars
{
    public:
    GintGpuVars(const std::shared_ptr<GintInfo> gint_info,
                const UnitCell& ucell,
                Numerical_Orbital* Phi);
    ~GintGpuVars();

    // ylmcoef_d is __constant__ memory
    double* ylmcoef_d = nullptr;
    int* atom_nw_d = nullptr;
    int* ucell_atom_nwl_d = nullptr;
    bool* atom_iw2_new_d = nullptr;
    int* atom_iw2_ylm_d = nullptr;
    int* atom_iw2_l_d = nullptr;
    double* psi_u_d = nullptr;
    double* dpsi_u_d = nullptr;
    double* d2psi_u_d = nullptr;
    double3* mgrid_pos_d = nullptr;
    int* iat2it_d = nullptr;
    int dev_id = 0;
    matrix_multiple_func_type fastest_matrix_mul;

}

}