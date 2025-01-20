#include "para_gemm.h"

#include "kernels/math_kernel_op.h"
#include "parallel_device.h"
namespace ModuleBase
{
template <typename T, typename Device>
PGemmCN<T, Device>::PGemmCN()
{
}
template <typename T, typename Device>
PGemmCN<T, Device>::~PGemmCN()
{
}

template <typename T, typename Device>
void PGemmCN<T, Device>::set_dimension(
#ifdef __MPI
    MPI_Comm comm_col,
    MPI_Comm comm_row,
#endif
    const int ncolA_in,
    const int LDA_in,
    const int ncolB_in,
    const int LDB_in,
    const int nrow_in,
    const int LDC_in,
    const bool gatherC_in)
{
#ifdef __MPI
    MPI_Comm_rank(comm_col, &col_rank);
    MPI_Comm_size(comm_col, &col_nproc);
    if (comm_row != MPI_COMM_NULL)
    {
        MPI_Comm_rank(comm_row, &row_rank);
        MPI_Comm_size(comm_row, &row_nproc);
    }
    col_world = comm_col;
    row_world = comm_row;
#endif
    this->LDA = LDA_in;
    this->LDB = LDB_in;
    this->LDC = LDC_in;
    this->ncolA = ncolA_in;
    this->ncolB = ncolB_in;
    this->nrow = nrow_in;
#ifdef __MPI
    this->gatherC = gatherC_in;
    requests.resize(col_nproc);
    colA_loc.resize(col_nproc);
    MPI_Allgather(&ncolA, 1, MPI_INT, colA_loc.data(), 1, MPI_INT, col_world);
    for (int ip = 0; ip < col_nproc; ip++)
    {
        max_colA = std::max(max_colA, colA_loc[ip]);
    }

    if (this->gatherC)
    {
        colB_loc.resize(col_nproc);
        recv_counts.resize(col_nproc);
        displs.resize(col_nproc);
        MPI_Allgather(&ncolB, 1, MPI_INT, colB_loc.data(), 1, MPI_INT, col_world);
        for (int ip = 0; ip < col_nproc; ip++)
        {
            recv_counts[ip] = LDC * colB_loc[ip];
        }
        displs[0] = 0;
        for (int ip = 1; ip < col_nproc; ip++)
        {
            displs[ip] = displs[ip - 1] + recv_counts[ip - 1];
        }
        size_C_global = displs[col_nproc - 1] + recv_counts[col_nproc - 1];
    }
    size_C_local = ncolB * LDC;
#endif
}

template <typename T, typename Device>
void PGemmCN<T, Device>::multiply(const T alpha, const T* A, const T* B, const T beta, T* C)
{
    const Device* ctx = {};
#ifdef __MPI
    if (col_nproc > 1)
    {
        std::vector<T> A_tmp(max_colA * LDA);
        for (int ip = 0; ip < col_nproc; ip++)
        {
            if (col_rank != ip)
            {
                int size = ncolA * LDA;
                Parallel_Common::isend_dev<T, Device>(A, size, ip, 0, col_world, &requests[ip], A_tmp.data());
            }
        }

        T* C_local = C;
        std::vector<T> C_tmp;
        if (this->gatherC)
        {
            C_tmp.resize(size_C_local);
            if (std::is_same<Device, base_device::DEVICE_GPU>::value)
            {
                C_local = nullptr;
                resmem_dev_op()(C_local, size_C_local);
            }
            else
            {
                C_local = C_tmp.data();
            }
            syncmem_dev_op()(C_local, C + displs[col_rank], size_C_local);
        }

        T* Atmp_device = nullptr;
        if (std::is_same<Device, base_device::DEVICE_GPU>::value)
        {
            resmem_dev_op()(Atmp_device, max_colA * LDA);
        }
        else
        {
            Atmp_device = A_tmp.data();
        }

        int shift = 0;
        T real_beta = row_rank == 0 ? beta : 0;
        for (int ip = 0; ip < col_nproc; ip++)
        {
            T* C_start = C_local + shift;
            if (col_rank == ip)
            {
                ModuleBase::gemm_op<T, Device>()(ctx,
                                                 'C',
                                                 'N',
                                                 ncolA,
                                                 ncolB,
                                                 nrow,
                                                 &alpha,
                                                 A,
                                                 LDA,
                                                 B,
                                                 LDB,
                                                 &real_beta,
                                                 C_start,
                                                 LDC);
                shift += ncolA;
            }
            else
            {
                int m = colA_loc[ip];
                int size = m * LDA;
                MPI_Status status;
                Parallel_Common::recv_dev<T, Device>(Atmp_device, size, ip, 0, col_world, &status, A_tmp.data());
                MPI_Wait(&requests[ip], &status);
                ModuleBase::gemm_op<T, Device>()(ctx,
                                                 'C',
                                                 'N',
                                                 m,
                                                 ncolB,
                                                 nrow,
                                                 &alpha,
                                                 Atmp_device,
                                                 LDA,
                                                 B,
                                                 LDB,
                                                 &real_beta,
                                                 C_start,
                                                 LDC);
                shift += m;
            }
        }

        if (this->gatherC)
        {
            T* Cglobal_cpu = nullptr;
            T* Clocal_cpu = C_tmp.data();;
            if (std::is_same<Device, base_device::DEVICE_GPU>::value)
            {
                delmem_dev_op()(Atmp_device);

                syncmem_d2h_op()(Clocal_cpu, C_local, size_C_local);
                delmem_dev_op()(C_local);
                
                resmem_dev_op()(Cglobal_cpu, size_C_global);
            }
            else
            {
                Cglobal_cpu = C;
            }
            if (this->row_nproc > 1)
            {
                Parallel_Common::reduce_data(Clocal_cpu, size_C_local, row_world);
            }
            Parallel_Common::gatherv_data(Clocal_cpu,
                                          size_C_local,
                                          Cglobal_cpu,
                                          recv_counts.data(),
                                          displs.data(),
                                          col_world);

            if (std::is_same<Device, base_device::DEVICE_GPU>::value)
            {
                syncmem_h2d_op()(C, Cglobal_cpu, size_C_global);
                delmem_dev_op()(Cglobal_cpu);
            }
        }
        else
        {
            if (this->row_nproc > 1)
            {
                Parallel_Common::reduce_dev<T, Device>(C, size_C_local, row_world);
            }
        }
    }
    else
    {
        T real_beta = row_rank == 0 ? beta : 0;
#else
    T real_beta = beta;
#endif
        ModuleBase::gemm_op<T, Device>()(ctx, 'C', 'N', ncolA, ncolB, nrow, &alpha, A, LDA, B, LDB, &real_beta, C, LDC);
#ifdef __MPI
        if (this->row_nproc > 1)
        {
            Parallel_Common::reduce_dev<T, Device>(C, size_C_local, row_world);
        }
    }
#endif
}

template class PGemmCN<double, base_device::DEVICE_CPU>;
template class PGemmCN<float, base_device::DEVICE_CPU>;
template class PGemmCN<std::complex<double>, base_device::DEVICE_CPU>;
template class PGemmCN<std::complex<float>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class PGemmCN<double, base_device::DEVICE_GPU>;
template class PGemmCN<float, base_device::DEVICE_GPU>;
template class PGemmCN<std::complex<double>, base_device::DEVICE_GPU>;
template class PGemmCN<std::complex<float>, base_device::DEVICE_GPU>;
#endif

} // namespace ModuleBase