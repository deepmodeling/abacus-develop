#include "para_linear_transform.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_hsolver/mpi_comm_helper.h"

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
    ModuleBase::timer::start("PLinearTransform", "act");
#ifdef __MPI
    if (nproc_col > 1)
    {
        syncmem_dev_op()(B_tmp_, B, ncolB * LDA);

        // Phase 1: Post all non-blocking sends
        MPIRequestTracker send_tracker;
        std::vector<MPI_Request> send_requests(nproc_col, MPI_REQUEST_NULL);
        for (int ip = 0; ip < nproc_col; ++ip)
        {
            if (rank_col != ip)
            {
                int size = LDA * ncolA;
                Parallel_Common::isend_dev<T, Device>(A, size, ip, 0, col_world,
                                                      &send_requests[ip], isend_tmp_.data());
            }
        }

        // Phase 2: Local computation (overlaps with sends in-flight)
        const int start = this->localU ? 0 : start_colB[rank_col];
        const T* U_part = U + start_colA[rank_col] + start * ncolA_glo;
        ModuleBase::matrixCopy<T, Device>()(ncolB, ncolA, U_part, ncolA_glo, U_tmp_, ncolA);
        ModuleBase::gemm_op<T, Device>()('N', 'N', nrowA, ncolB, ncolA,
                                         &alpha, A, LDA, U_tmp_, ncolA, &beta, B, LDA);

        // Phase 3: Post non-blocking receives and process remote data
        T* Atmp_device = nullptr;
        if (std::is_same<Device, base_device::DEVICE_GPU>::value)
        {
            Atmp_device = A_tmp_device_;
        }
        else
        {
            Atmp_device = A_tmp_.data();
        }

        MPIRequestTracker recv_tracker;
        for (int ip = 0; ip < nproc_col; ++ip)
        {
            if (ip != rank_col)
            {
                const int ncolA_ip = colA_loc[ip];
                const T* U_part_ip = U + start_colA[ip] + start * ncolA_glo;
                // Copy U partition (independent of recv, can be done while waiting)
                ModuleBase::matrixCopy<T, Device>()(ncolB, ncolA_ip, U_part_ip,
                                                    ncolA_glo, U_tmp_, ncolA_ip);

                int size = LDA * ncolA_ip;
                // Use non-blocking receive
                MPI_Request recv_req;
                MPI_Irecv(Atmp_device, size,
                          (std::is_same<T, std::complex<double>>::value) ? MPI_DOUBLE_COMPLEX
                          : (std::is_same<T, std::complex<float>>::value) ? MPI_C_FLOAT_COMPLEX
                          : MPI_DOUBLE,
                          ip, 0, col_world, &recv_req);
                recv_tracker.add(recv_req);

                // Wait for this receive before using the data
                MPI_Wait(&recv_req, MPI_STATUS_IGNORE);

                T zero = 0.0;
                ModuleBase::gemm_op<T, Device>()('N', 'N', nrowA, ncolB, ncolA_ip,
                                                 &alpha, Atmp_device, LDA,
                                                 U_tmp_, ncolA_ip, &zero, B_tmp_, LDA);
                // Accumulate into B
                T one = 1.0;
                ModuleBase::axpy_op<T, Device>()(ncolB * LDA, &one, B_tmp_, 1, B, 1);
            }
        }

        // Phase 4: Wait for all sends to complete
        for (int ip = 0; ip < nproc_col; ++ip)
        {
            if (rank_col != ip && send_requests[ip] != MPI_REQUEST_NULL)
            {
                MPI_Status status;
                MPI_Wait(&send_requests[ip], &status);
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
    ModuleBase::timer::end("PLinearTransform", "act");
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