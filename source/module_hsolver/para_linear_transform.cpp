#include "para_linear_transform.h"
#include <vector>
#include <algorithm>
namespace hsolver
{
template <typename T, typename Device>
void para_linear_transform_op<T, Device>::operator()(T* A,
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
                                                     const int nproc_col

)
{
    const Device* ctx = {};
#ifdef __MPI
    if (nproc_col > 1)
    {
        std::vector<int> colA_loc(nproc_col);
        MPI_Allgather(&ncol_loc, 1, MPI_INT, colA_loc.data(), 1, MPI_INT, col_world);
        std::vector<int> start_col(nproc_col);
        start_col[0] = 0;
        for (int ip = 1; ip < nproc_col; ++ip)
        {
            start_col[ip] = start_col[ip - 1] + colA_loc[ip - 1];
        }
        int max_col = *std::max_element(colA_loc.begin(), colA_loc.end());
        std::vector<MPI_Request> requests(nproc_col);

        std::vector<T> A_tmp(max_col * LDA);
        T* A_tmp_device = A_tmp.data();
        if (std::is_same<Device, base_device::DEVICE_GPU>::value)
        {
            A_tmp_device = nullptr;
            resmem_dev_op()(A_tmp_device, max_col * LDA);
        }
        T* A_tmp2 = nullptr;
        resmem_dev_op()(A_tmp2, ncol_loc * LDA);
        syncmem_dev_op()(A_tmp2, A, ncol_loc * LDA);
        T* A_sum = nullptr;
        resmem_dev_op()(A_sum, ncol_loc * LDA);
        setmem_dev_op()(A_sum, 0.0, ncol_loc * LDA);

        // Send
        for (int ip = 0; ip < nproc_col; ++ip)
        {
            if (rank_col != ip)
            {
                int size = LDA * ncol_loc;
                Parallel_Common::isend_dev<T, Device>(A, size, ip, 0, col_world, &requests[ip], A_tmp.data());
            }
        }

        // Receive
        T* U_local = nullptr;
        resmem_dev_op()(U_local, max_col * ncol_loc);
        const int start = start_col[rank_col];
        for (int ip = 0; ip < nproc_col; ++ip)
        {
            T real_beta = ip == 0 ? beta : 0;
            const int start_row = start_col[ip];
            const int ncol_ip = colA_loc[ip];
            // get U_local
            for (int i = 0; i < ncol_loc; ++i)
            {
                const T* U_glo_tmp = U_global + start_row + (i + start) * ncol_glo;
                syncmem_dev_op()(U_local + i * ncol_ip, U_glo_tmp, ncol_ip);
            }

            if (ip == rank_col)
            {
                ModuleBase::gemm_op<T, Device>()(ctx,
                                                 'N',
                                                 'N',
                                                 nrow,
                                                 ncol_loc,
                                                 ncol_ip,
                                                 &alpha,
                                                 A,
                                                 LDA,
                                                 U_local,
                                                 ncol_ip,
                                                 &real_beta,
                                                 A_tmp2,
                                                 LDA);
            }
            else
            {
                int size = LDA * ncol_ip;
                MPI_Status status;
                Parallel_Common::recv_dev<T, Device>(A_tmp_device, size, ip, 0, col_world, &status, A_tmp.data());
                MPI_Wait(&requests[ip], &status);
                ModuleBase::gemm_op<T, Device>()(ctx,
                                                 'N',
                                                 'N',
                                                 nrow,
                                                 ncol_loc,
                                                 ncol_ip,
                                                 &alpha,
                                                 A_tmp_device,
                                                 LDA,
                                                 U_local,
                                                 ncol_ip,
                                                 &real_beta,
                                                 A_tmp2,
                                                 LDA);
            }
            // sum all the results
            T one = 1.0;
            ModuleBase::axpy_op<T, Device>()(ctx, ncol_loc * LDA, &one, A_tmp2, 1, A_sum, 1);
        }
        syncmem_dev_op()(A, A_sum, ncol_loc * LDA);
        delmem_dev_op()(U_local);
        delmem_dev_op()(A_tmp2);
        delmem_dev_op()(A_sum);
        if (std::is_same<Device, base_device::DEVICE_GPU>::value)
        {
            delmem_dev_op()(A_tmp_device);
        }
    }
    else
#endif
    {
        T* A_tmp = nullptr;
        resmem_dev_op()(A_tmp, LDA * ncol_glo);
        syncmem_dev_op()(A_tmp, A, LDA * ncol_loc);
        ModuleBase::gemm_op<T, Device>()(ctx,
                                         'N',
                                         'N',
                                         nrow,
                                         ncol_glo,
                                         ncol_glo,
                                         &alpha,
                                         A_tmp,
                                         LDA,
                                         U_global,
                                         ncol_glo,
                                         &beta,
                                         A,
                                         LDA);
        delmem_dev_op()(A_tmp);
    }
};

template struct para_linear_transform_op<double, base_device::DEVICE_CPU>;
template struct para_linear_transform_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct para_linear_transform_op<std::complex<float>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template struct para_linear_transform_op<double, base_device::DEVICE_GPU>;
template struct para_linear_transform_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct para_linear_transform_op<std::complex<float>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver