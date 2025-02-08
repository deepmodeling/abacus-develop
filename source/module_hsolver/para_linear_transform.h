#ifndef __PARA_LINEAR_TRANSFORM_H__
#define __PARA_LINEAR_TRANSFORM_H__
#include "module_base/kernels/math_kernel_op.h"
#include "module_base/module_device/device.h"
#include "module_base/module_device/memory_op.h"
#include "module_base/parallel_device.h"
#include <vector>
#ifdef __MPI
#include "mpi.h"
#endif
namespace hsolver
{

/**
 * @brief B =  alpha * A * U + beta * B
 *        A and B are local matrice
 *        U can be a local matrix or a global matrix
 */
template <typename T, typename Device>
class PLinearTransform
{
  public:
    using syncmem_dev_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using resmem_dev_op = base_device::memory::resize_memory_op<T, Device>;
    using setmem_dev_op = base_device::memory::set_memory_op<T, Device>;
    using delmem_dev_op = base_device::memory::delete_memory_op<T, Device>;
    int nproc_col = 1;
    int rank_col = 0;
    int nrowA = 0;
    int ncolA = 0;
    int ncolB = 0;
    int LDA = 0;
    bool localU = false;
#ifdef __MPI
    MPI_Comm col_world = MPI_COMM_NULL;
    std::vector<int> colA_loc;
    std::vector<int> start_colA;
    std::vector<int> start_colB;
    int max_colA = 0;
    int ncolA_glo = 0;
    int max_colB = 0;
#endif

    /**
     * @brief set the dimension of A, B, and U
     *        A: LDA * nrow, U_global: ncolA_global * ncolB_global, U_local: ncolA_global * ncolB
     *        B: LDA * ncolB
     */
    void set_dimension(const int nrowA,
                       const int ncolA,
                       const int ncolB,
                       const int LDA,
#ifdef __MPI
                       MPI_Comm col_world,
#endif
                       const bool localU);

    /**
     * @brief B =  alpha * A * U + beta * B
     *        A is a local matrix with nrow rows and ncolA_loc columns
     *        B is a local matrix with nrow rows and ncolB_loc columns
     *        U can be a local matrix or a global matrix
     * @example rotate wave functions: B = A * U
     *          orthogonalize wave functions: B = - A * U + B
     *
     * @param alpha : alpha
     * @param A : input matrix
     * @param U_global : input matrix
     * @param beta : beta
     * @param B : input/output matrix
     *
     */
    void act(const T alpha, const T* A, const T* U_global, const T beta, T* B);
};
} // namespace hsolver
#endif