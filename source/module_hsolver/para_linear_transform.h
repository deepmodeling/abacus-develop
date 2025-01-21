#ifndef __PARA_LINEAR_TRANSFORM_H__
#define __PARA_LINEAR_TRANSFORM_H__
#include "module_base/kernels/math_kernel_op.h"
#include "module_base/module_device/device.h"
#include "module_base/module_device/memory_op.h"
#include "module_base/parallel_device.h"
#ifdef __MPI
#include "mpi.h"
#endif
namespace hsolver
{

template <typename T, typename Device>
struct para_linear_transform_op
{
    using syncmem_dev_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using resmem_dev_op = base_device::memory::resize_memory_op<T, Device>;
    using setmem_dev_op = base_device::memory::set_memory_op<T, Device>;
    using delmem_dev_op = base_device::memory::delete_memory_op<T, Device>;
    /**
     * @brief A_global =  alpha * A_global * U_global + beta * A_global
     *        A is a local matrix with nrow rows and ncol_loc columns
     *        U_global is a matrix with ncol_glo rows and ncol_glo columns
     * @example rotate wave functions: A = A * U
     *          orthogonalize wave functions: A = A - A * U
     *
     * @param A : input/output matrix
     * @param alpha : alpha
     * @param beta : beta
     * @param U_global : input matrix
     * @param nrow : number of rows of A
     * @param LDA : leading dimension of A
     * @param ncol_loc : number of columns of A
     * @param ncol_glo : number of columns and rows of U_global
     * @param col_world : column communicator world
     * @param rank_col : rank of col_world
     * @param nproc_col : number of processes in col_world
     *
     */
    void operator()(T* A,
                    const T alpha,
                    const T beta,
                    const T* U_global,
                    const int& nrow,
                    const int& LDA,
                    const int& ncol_loc,
                    const int& ncol_glo,
#ifdef __MPI
                    MPI_Comm col_world,
#endif
                    const int rank_col,
                    const int nproc_col);
};
} // namespace hsolver
#endif