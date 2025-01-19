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
 *        C_global = alpha * A^+ * B + beta * C_global
 *        Here, A and B are local matrices in each proc, and C_global is a global matrix gathered from all procs
 *        All procs have their own C_global matrix with the same values.
 */
template <typename T, typename Device = base_device::DEVICE_CPU>
class PGemmCN
{
  public:
    PGemmCN();
    ~PGemmCN();

    /**
     * @brief set the dimension of A, B, and C_global
     *
     * @param ncolA number of columns of A, which is a local matrix in each proc
     * @param LDA leading dimension of A in each proc
     * @param ncolB number of columns of B, which is a local matrix in each proc
     * @param LDB leading dimension of B in each proc
     * @param nrow number of rows of A or B
     * @param LDC_global leading dimension of C_global, which is the global C matrix gathered from all procs
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
        const int LDC_global);
    /**
     * @brief calculate C_global = alpha * A^+ * B + beta * C_global
     *
     * @param alpha
     * @param A
     * @param B
     * @param beta
     * @param C_global
     */
    void multiply(const T alpha, const T* A, const T* B, const T beta, T* C_global);
#ifdef __MPI
    MPI_Comm col_world = MPI_COMM_NULL; ///< column communicator world
    MPI_Comm row_world = MPI_COMM_NULL; ///< row communicator world

    int col_rank = 0;  ///< rank in col_world
    int col_nproc = 1; ///< number of procs in col_world
    int row_rank = 0;  ///< rank in row_world

    std::vector<int> colA_loc; ///< [col_nproc] number of columns of A matrix in each proc
    int max_colA = 0;          ///< maximum number of columns of A matrix in all procs
    std::vector<int> colB_loc; ///<[col_nproc] number of columns of B matrix in each proc
    std::vector<int> row_loc;  ///<[col_nproc] number of rows of C matrix in each proc

    std::vector<MPI_Request> requests; ///< MPI request
    std::vector<int> recv_counts;      ///< receive counts for gathering C_local to C_global
    std::vector<int> displs;           ///< displacements for gathering C_local to C_global
    int send_counts = 0;               ///< send counts for gathering C_local to C_global
    int size_C_global = 0;             ///< size of C_global, which is the global C matrix gathered from all procs
#endif
    int ncolA = 0;      ///< number of columns of A, which is a local matrix in each proc
    int ncolB = 0;      ///< number of columns of B, which is a local matrix in each proc
    int nrow = 0;       ///< number of rows of A or B
    int LDA = 0;        ///< leading dimension of A in each proc
    int LDB = 0;        ///< leading dimension of B in each proc
    int LDC_global = 0; ///< leading dimension of C_global, which is the global C matrix gathered from all procs
  private:
    using resmem_dev_op = base_device::memory::resize_memory_op<T, Device>; 
    using delmem_dev_op = base_device::memory::delete_memory_op<T, Device>;
    using syncmem_dev_op = base_device::memory::synchronize_memory_op<T, Device, Device>;

};
} // namespace ModuleBase
#endif