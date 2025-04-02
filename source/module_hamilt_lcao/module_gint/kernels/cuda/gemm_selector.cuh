#ifndef GEMM_SELECTOR_H
#define GEMM_SELECTOR_H

#include "module_cell/unitcell.h"
#include "cuda_runtime.h"
typedef std::function<void(int,
                           int,
                           int*,
                           int*,
                           int*,
                           double**,
                           int*,
                           double**,
                           int*,
                           double**,
                           int*,
                           int,
                           cudaStream_t,
                           double* alpha)>
matrix_multiple_func_type;

void gemm_algo_selector(int k, matrix_multiple_func_type& func,const UnitCell& ucell);

#endif