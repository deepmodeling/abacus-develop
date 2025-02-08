#include "para_linear_transform.h"

#include <algorithm>
#include <vector>
namespace hsolver
{
template <typename T, typename Device>
void PLinearTransform<T, Device>::set_dimension(const int nrowA,
                                                const int ncolA,
                                                const int ncolB,
                                                const int LDA,
#ifdef __MPI
                                                MPI_Comm col_world,
#endif
                                                const bool localU)
{
    this->nrowA = nrowA;
    this->ncolA = ncolA;
    this->ncolB = ncolB;
    this->LDA = LDA;
#ifdef __MPI
    this->col_world = col_world;
    MPI_Comm_rank(col_world, &rank_col);
    MPI_Comm_size(col_world, &nproc_col);
    if (nproc_col > 1)
    {
        this->localU = localU;
        colA_loc.resize(nproc_col);
        MPI_Allgather(&ncolA, 1, MPI_INT, colA_loc.data(), 1, MPI_INT, col_world);
        start_colA.resize(nproc_col);
        start_colA[0] = 0;
        for (int ip = 1; ip < nproc_col; ++ip)
        {
            start_colA[ip] = start_colA[ip - 1] + colA_loc[ip - 1];
        }
        this->ncolA_glo = start_colA[nproc_col - 1] + colA_loc[nproc_col - 1];
        this->max_colA = *std::max_element(colA_loc.begin(), colA_loc.end());

        std::vector<int> colB_loc(nproc_col);
        MPI_Allgather(&ncolB, 1, MPI_INT, colB_loc.data(), 1, MPI_INT, col_world);
        start_colB.resize(nproc_col);
        start_colB[0] = 0;
        for (int ip = 1; ip < nproc_col; ++ip)
        {
            start_colB[ip] = start_colB[ip - 1] + colB_loc[ip - 1];
        }
        this->max_colB = *std::max_element(colB_loc.begin(), colB_loc.end());
    }
#else
    nproc_col = 1;
    rank_col = 0;
#endif
}
template <typename T, typename Device>
void PLinearTransform<T, Device>::act(const T alpha, const T* A, const T* U, const T beta, T* B)
{
    const Device* ctx = {};
#ifdef __MPI
    if (nproc_col > 1)
    {
        std::vector<MPI_Request> requests(nproc_col);
        std::vector<T> A_tmp(max_colA * LDA);
        std::vector<T> isend_tmp;
        T* A_tmp_device = A_tmp.data();
        if (std::is_same<Device, base_device::DEVICE_GPU>::value)
        {
            A_tmp_device = nullptr;
            isend_tmp.resize(max_colA * LDA);
            resmem_dev_op()(A_tmp_device, max_colA * LDA);
        }
        T* B_tmp = nullptr;
        resmem_dev_op()(B_tmp, ncolB * LDA);
        syncmem_dev_op()(B_tmp, B, ncolB * LDA);
        setmem_dev_op()(B, 0.0, ncolB * LDA);

        T* U_tmp = nullptr;
        resmem_dev_op()(U_tmp, max_colA * max_colB);

        // Send
        for (int ip = 0; ip < nproc_col; ++ip)
        {
            if (rank_col != ip)
            {
                int size = LDA * ncolA;
                Parallel_Common::isend_dev<T, Device>(A, size, ip, 0, col_world, &requests[ip], isend_tmp.data());
            }
        }

        // Receive
        const int start = this->localU ? 0 : start_colB[rank_col];
        for (int ip = 0; ip < nproc_col; ++ip)
        {
            T real_beta = ip == 0 ? beta : 0;
            const int ncolA_ip = colA_loc[ip];
            // get U_tmp

            const int start_row = start_colA[ip];
            for (int i = 0; i < ncolB; ++i)
            {
                const T* U_part = U + start_row + (i + start) * ncolA_glo;
                syncmem_dev_op()(U_tmp + i * ncolA_ip, U_part, ncolA_ip);
            }

            if (ip == rank_col)
            {
                ModuleBase::gemm_op<T, Device>()(ctx,
                                                 'N',
                                                 'N',
                                                 nrowA,
                                                 ncolB,
                                                 ncolA_ip,
                                                 &alpha,
                                                 A,
                                                 LDA,
                                                 U_tmp,
                                                 ncolA_ip,
                                                 &real_beta,
                                                 B_tmp,
                                                 LDA);
            }
            else
            {
                int size = LDA * ncolA_ip;
                MPI_Status status;
                Parallel_Common::recv_dev<T, Device>(A_tmp_device, size, ip, 0, col_world, &status, A_tmp.data());
                MPI_Wait(&requests[ip], &status);
                ModuleBase::gemm_op<T, Device>()(ctx,
                                                 'N',
                                                 'N',
                                                 nrowA,
                                                 ncolB,
                                                 ncolA_ip,
                                                 &alpha,
                                                 A_tmp_device,
                                                 LDA,
                                                 U_tmp,
                                                 ncolA_ip,
                                                 &real_beta,
                                                 B_tmp,
                                                 LDA);
            }
            // sum all the results
            T one = 1.0;
            ModuleBase::axpy_op<T, Device>()(ctx, ncolB * LDA, &one, B_tmp, 1, B, 1);
        }
        delmem_dev_op()(U_tmp);
        delmem_dev_op()(B_tmp);
        if (std::is_same<Device, base_device::DEVICE_GPU>::value)
        {
            delmem_dev_op()(A_tmp_device);
        }
    }
    else
#endif
    {
        ModuleBase::gemm_op<T, Device>()(ctx,
                                         'N',
                                         'N',
                                         nrowA,
                                         ncolB,
                                         ncolA,
                                         &alpha,
                                         A,
                                         LDA,
                                         U,
                                         ncolA,
                                         &beta,
                                         B,
                                         LDA);
    }
};

template struct PLinearTransform<double, base_device::DEVICE_CPU>;
template struct PLinearTransform<std::complex<double>, base_device::DEVICE_CPU>;
template struct PLinearTransform<std::complex<float>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template struct PLinearTransform<double, base_device::DEVICE_GPU>;
template struct PLinearTransform<std::complex<double>, base_device::DEVICE_GPU>;
template struct PLinearTransform<std::complex<float>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver