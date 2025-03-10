#include "para_linear_transform.h"

#include "module_base/timer.h"

#include <algorithm>
#include <vector>
namespace hsolver
{
template <typename T, typename Device>
PLinearTransform<T, Device>::~PLinearTransform()
{
#ifdef __MPI
    delmem_dev_op()(U_tmp_);
    delmem_dev_op()(B_tmp_);
    delmem_dev_op()(A_tmp_device_);
#endif
}
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

        // allocate temperory memory
        resmem_dev_op()(B_tmp_, ncolB * LDA);
        resmem_dev_op()(U_tmp_, max_colA * max_colB);
        if (std::is_same<Device, base_device::DEVICE_GPU>::value)
        {
            resmem_dev_op()(A_tmp_device_, max_colA * LDA);
#ifndef __CUDA_MPI
            isend_tmp_.resize(max_colA * LDA);
#endif
        }
        A_tmp_.resize(max_colA * LDA);
    }
#else
    nproc_col = 1;
    rank_col = 0;
#endif
}
template <typename T, typename Device>
void PLinearTransform<T, Device>::act(const T alpha, const T* A, const T* U, const T beta, T* B)
{
    ModuleBase::timer::tick("PLinearTransform", "act");
#ifdef __MPI
    if (nproc_col > 1)
    {
        syncmem_dev_op()(B_tmp_, B, ncolB * LDA);
        std::vector<MPI_Request> requests(nproc_col);
        // Send
        for (int ip = 0; ip < nproc_col; ++ip)
        {
            if (rank_col != ip)
            {
                int size = LDA * ncolA;
                Parallel_Common::isend_dev<T, Device>(A, size, ip, 0, col_world, &requests[ip], isend_tmp_.data());
            }
        }

        // local part
        const int start = this->localU ? 0 : start_colB[rank_col];
        const T* U_part = U + start_colA[rank_col] + start * ncolA_glo;
        ModuleBase::matrixCopy<T, Device>()(ncolB, ncolA, U_part, ncolA_glo, U_tmp_, ncolA);
        ModuleBase::gemm_op<T, Device>()('N', 'N', nrowA, ncolB, ncolA, &alpha, A, LDA, U_tmp_, ncolA, &beta, B, LDA);

        // Receive
        T* Atmp_device = nullptr;
        if (std::is_same<Device, base_device::DEVICE_GPU>::value)
        {
            Atmp_device = A_tmp_device_;
        }
        else
        {
            Atmp_device = A_tmp_.data();
        }
        for (int ip = 0; ip < nproc_col; ++ip)
        {
            if (ip != rank_col)
            {
                T zero = 0.0;
                const int ncolA_ip = colA_loc[ip];
                const T* U_part = U + start_colA[ip] + start * ncolA_glo;
                ModuleBase::matrixCopy<T, Device>()(ncolB, ncolA_ip, U_part, ncolA_glo, U_tmp_, ncolA_ip);

                int size = LDA * ncolA_ip;
                MPI_Status status;
#ifdef __CUDA_MPI
                // If the memory is not set to zero, it may cause the result to be wrong when using CUDA Aware MPI
                // I am not sure if it is due to CUDA Aware MPI or not
                base_device::memory::set_memory_op<T, Device>()(Atmp_device, 0, size);
#endif
                Parallel_Common::recv_dev<T, Device>(Atmp_device, size, ip, 0, col_world, &status, A_tmp_.data());
                ModuleBase::gemm_op<T, Device>()('N',
                                                 'N',
                                                 nrowA,
                                                 ncolB,
                                                 ncolA_ip,
                                                 &alpha,
                                                 Atmp_device,
                                                 LDA,
                                                 U_tmp_,
                                                 ncolA_ip,
                                                 &zero,
                                                 B_tmp_,
                                                 LDA);
                // sum all the results
                T one = 1.0;
                ModuleBase::axpy_op<T, Device>()(ncolB * LDA, &one, B_tmp_, 1, B, 1);
            }
        }

        for (int ip = 0; ip < nproc_col; ++ip)
        {
            if (rank_col != ip)
            {
                MPI_Status status;
                MPI_Wait(&requests[ip], &status);
            }
        }
    }
    else
#endif
    {
        ModuleBase::gemm_op<T, Device>()('N',
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
    ModuleBase::timer::tick("PLinearTransform", "act");
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