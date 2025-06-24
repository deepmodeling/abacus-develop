#pragma once

#include <cuda_runtime.h>
#include "set_const_mem.cuh"
#include "source_base/ylm.h"
#include "source_cell/unitcell.h"
#include "source_cell/atom_spec.h"
#include "module_hamilt_lcao/module_gint/temp_gint/biggrid_info.h"
#include "gint_helper.cuh"
#include "module_hamilt_lcao/module_gint/kernels/cuda/gemm_selector.cuh"

namespace ModuleGint
{

class GintGpuVars
{
    public:
    GintGpuVars(std::shared_ptr<const BigGridInfo> bgrid_info,
                const UnitCell& ucell,
                const Numerical_Orbital* Phi);
    ~GintGpuVars();
    
    int nwmax;
    double dr_uniform;
    double nr_max;
    // ylmcoef_d is __constant__ memory, no need to cudaFree
    double* ylmcoef_d = nullptr;
    double* rcut_d = nullptr;
    int* atom_nw_d = nullptr;
    int* ucell_atom_nwl_d = nullptr;
    bool* atom_iw2_new_d = nullptr;
    int* atom_iw2_ylm_d = nullptr;
    int* atom_iw2_l_d = nullptr;
    double* psi_u_d = nullptr;
    double* dpsi_u_d = nullptr;
    double* d2psi_u_d = nullptr;
    double3* mgrids_pos_d = nullptr;
    int* iat2it_d = nullptr;

    // the index of gpu device
    int dev_id_ = 0;
    matrix_multiple_func_type fastest_matrix_mul;

};

}