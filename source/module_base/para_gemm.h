#ifndef PARA_GEMM_H
#define PARA_GEMM_H
#include "module_base/module_device/device.h"
#include "module_base/module_device/memory_op.h"

#include <vector>
#ifdef __MPI
#include "mpi.h"
#endif

namespace ModuleBase
{
/**
 * @brief this class is used to perform parallel matrix multiplication
 *        C = alpha * A^H * B + beta * C
 *        Here, A and B are local matrices in each proc,
 *        C can be C_local or C_global, depending on the value of gatherC
 *        C_local is a local matrix in each proc
 *        C_global is a global matrix gathered from all procs and all procs have their own C_global matrix with the same
 *        C_global and C_local have the same LDC, but different column numbers
 * values.
 */
template <typename T, typename Device = base_device::DEVICE_CPU>
class PGemmCN
{
  public:
    PGemmCN();
    ~PGemmCN();

    /**
     * @brief set the dimension of A, B, and C
     *
     * @param ncolA number of columns of A, which is a local matrix in each proc
     * @param LDA leading dimension of A in each proc
     * @param ncolB number of columns of B, which is a local matrix in each proc
     * @param LDB leading dimension of B in each proc
     * @param nrow number of rows of A or B
     * @param LDC leading dimension of C. C can be C_local or C_global
     * @param mode 1: gather C_local to C_global, 2:C_local(nrow * ncol_loc), 3:C_global(nrow_loc * ncol)
     */
    void set_dimension(
#ifdef __MPI
        MPI_Comm comm_col,
        MPI_Comm comm_row,
#endif
        const int ncolA,
        const int LDA,
        const int ncolB,
        const int LDB,
        const int nrow,
        const int LDC,
        const int mode = 1);

    /**
     * @brief calculate C = alpha * A^H * B + beta * C
     *
     */
    void multiply(const T alpha, const T* A, const T* B, const T beta, T* C);
#ifdef __MPI
    MPI_Comm col_world = MPI_COMM_NULL; ///< column communicator world
    MPI_Comm row_world = MPI_COMM_NULL; ///< row communicator world

    int col_rank = 0;  ///< rank in col_world
    int col_nproc = 1; ///< number of procs in col_world
    int row_rank = 0;  ///< rank in row_world
    int row_nproc = 1; ///< number of procs in row_world

    std::vector<int> colA_loc; ///< [col_nproc] number of columns of A matrix in each proc
    int max_colA = 0;          ///< maximum number of columns of A matrix in all procs
    std::vector<int> colB_loc; ///< [col_nproc] number of columns of B matrix in each proc
    int max_colB = 0;          ///< maximum number of columns of B matrix in all procs

    std::vector<MPI_Request> requests; ///< MPI request
    std::vector<int> recv_counts;      ///< receive counts for gathering C_local to C_global
    std::vector<int> displs;           ///< displacements for gathering C_local to C_global
    int size_C_local = 0;              ///< size of C_local, which is a local matrix in each proc
    int size_C_global = 0;             ///< size of C_global, which is the global C matrix gathered from all procs
    bool gatherC = true;               ///< whether gather C_local to C_global
    bool divideCrow = false;           ///< whether divide C_global to C_local
#endif
    int ncolA = 0; ///< number of columns of A, which is a local matrix in each proc
    int ncolB = 0; ///< number of columns of B, which is a local matrix in each proc
    int nrow = 0;  ///< number of rows of A or B
    int LDA = 0;   ///< leading dimension of A in each proc
    int LDB = 0;   ///< leading dimension of B in each proc
    int LDC = 0;   ///< leading dimension of C, which can be C_local or C_global
  private:
    /// @brief for col_nproc == 1
    void multiply_single(const T alpha, const T* A, const T* B, const T beta, T* C);
#ifdef __MPI
    /// @brief for mode = 1 or 2
    void multiply_col(const T alpha, const T* A, const T* B, const T beta, T* C);
    /// @brief for mode = 3
    void multiply_row(const T alpha, const T* A, const T* B, const T beta, T* C);
#endif
    using resmem_dev_op = base_device::memory::resize_memory_op<T, Device>;
    using delmem_dev_op = base_device::memory::delete_memory_op<T, Device>;
    using syncmem_dev_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_d2h_op = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;
    using syncmem_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
};
} // namespace ModuleBase
#endif