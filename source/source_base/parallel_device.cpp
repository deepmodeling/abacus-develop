#include "parallel_device.h"

#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/types.h"

#if defined(__MPI) && defined(__NCCL_PARALLEL_DEVICE)
#include "source_base/module_device/device_check.h"

#include <cuda_runtime.h>

#include <map>
#include <mutex>

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#endif
#ifdef __MPI
namespace Parallel_Common
{
#if defined(__NCCL_PARALLEL_DEVICE)
namespace
{
struct NcclCommContext
{
    ncclComm_t comm = nullptr;
    cudaStream_t stream = nullptr;
    int size = 0;
};

class NcclCommRegistry
{
  public:
    ~NcclCommRegistry()
    {
        for (std::map<MPI_Fint, NcclCommContext>::iterator it = contexts_.begin(); it != contexts_.end(); ++it)
        {
            if (it->second.comm != nullptr)
            {
                ncclCommDestroy(it->second.comm);
            }
        }
    }

    NcclCommContext& get(MPI_Comm comm)
    {
        const MPI_Fint key = MPI_Comm_c2f(comm);
        std::lock_guard<std::mutex> lock(mutex_);
        std::map<MPI_Fint, NcclCommContext>::iterator found = contexts_.find(key);
        if (found != contexts_.end())
        {
            return found->second;
        }

        int rank = 0;
        int size = 0;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);

        NcclCommContext ctx;
        ctx.size = size;
        if (size > 1)
        {
            ncclUniqueId id;
            if (rank == 0)
            {
                CHECK_NCCL(ncclGetUniqueId(&id));
            }
            MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, comm);
            CHECK_NCCL(ncclCommInitRank(&ctx.comm, size, id, rank));
        }

        std::pair<std::map<MPI_Fint, NcclCommContext>::iterator, bool> inserted = contexts_.insert(std::make_pair(key, ctx));
        return inserted.first->second;
    }

  private:
    std::map<MPI_Fint, NcclCommContext> contexts_;
    std::mutex mutex_;
};

NcclCommRegistry& get_nccl_registry()
{
    static NcclCommRegistry registry;
    return registry;
}

template <typename T>
void nccl_bcast_impl(T* object, const int n, MPI_Comm& comm, ncclDataType_t datatype, const int count_scale = 1)
{
    NcclCommContext& ctx = get_nccl_registry().get(comm);
    if (ctx.size <= 1 || n <= 0)
    {
        return;
    }
    CHECK_NCCL(ncclBroadcast(object, object, static_cast<size_t>(n) * count_scale, datatype, 0, ctx.comm, ctx.stream));
    CHECK_CUDA(cudaStreamSynchronize(ctx.stream));
}

template <typename T>
void nccl_reduce_impl(T* object, const int n, MPI_Comm& comm, ncclDataType_t datatype, const int count_scale = 1)
{
    NcclCommContext& ctx = get_nccl_registry().get(comm);
    if (ctx.size <= 1 || n <= 0)
    {
        return;
    }
    CHECK_NCCL(ncclAllReduce(object, object, static_cast<size_t>(n) * count_scale, datatype, ncclSum, ctx.comm, ctx.stream));
    CHECK_CUDA(cudaStreamSynchronize(ctx.stream));
}

template <typename T>
void nccl_gatherv_impl(const T* sendbuf,
                       const int sendcount,
                       T* recvbuf,
                       const int* recvcounts,
                       const int* displs,
                       MPI_Comm& comm)
{
    NcclCommContext& ctx = get_nccl_registry().get(comm);
    if (ctx.size <= 1)
    {
        if (sendbuf != recvbuf && sendcount > 0)
        {
            CHECK_CUDA(cudaMemcpy(recvbuf, sendbuf, static_cast<size_t>(sendcount) * sizeof(T), cudaMemcpyDeviceToDevice));
        }
        return;
    }

    int chunk_count = 0;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    for (int i = 0; i < ctx.size; ++i)
    {
        if (recvcounts[i] > chunk_count)
        {
            chunk_count = recvcounts[i];
        }
    }
    if (recvcounts[rank] != sendcount)
    {
        throw std::runtime_error("nccl_gatherv_data: sendcount does not match recvcounts[rank]");
    }
    if (chunk_count <= 0)
    {
        return;
    }

    const size_t chunk_bytes = static_cast<size_t>(chunk_count) * sizeof(T);
    const size_t recv_bytes = chunk_bytes * ctx.size;
    unsigned char* staged_send = nullptr;
    unsigned char* staged_recv = nullptr;

    CHECK_CUDA(cudaMalloc(&staged_send, chunk_bytes));
    CHECK_CUDA(cudaMalloc(&staged_recv, recv_bytes));
    CHECK_CUDA(cudaMemsetAsync(staged_send, 0, chunk_bytes, ctx.stream));
    if (sendcount > 0)
    {
        CHECK_CUDA(cudaMemcpyAsync(staged_send,
                                   sendbuf,
                                   static_cast<size_t>(sendcount) * sizeof(T),
                                   cudaMemcpyDeviceToDevice,
                                   ctx.stream));
    }

    CHECK_NCCL(ncclAllGather(staged_send, staged_recv, chunk_bytes, ncclUint8, ctx.comm, ctx.stream));

    for (int i = 0; i < ctx.size; ++i)
    {
        if (recvcounts[i] > 0)
        {
            CHECK_CUDA(cudaMemcpyAsync(recvbuf + displs[i],
                                       staged_recv + static_cast<size_t>(i) * chunk_bytes,
                                       static_cast<size_t>(recvcounts[i]) * sizeof(T),
                                       cudaMemcpyDeviceToDevice,
                                       ctx.stream));
        }
    }

    CHECK_CUDA(cudaStreamSynchronize(ctx.stream));
    CHECK_CUDA(cudaFree(staged_send));
    CHECK_CUDA(cudaFree(staged_recv));
}
} // namespace

void nccl_bcast_data(double* object, const int& n, MPI_Comm& comm)
{
    nccl_bcast_impl(object, n, comm, ncclDouble);
}

void nccl_bcast_data(std::complex<double>* object, const int& n, MPI_Comm& comm)
{
    nccl_bcast_impl(reinterpret_cast<double*>(object), n, comm, ncclDouble, 2);
}

void nccl_bcast_data(float* object, const int& n, MPI_Comm& comm)
{
    nccl_bcast_impl(object, n, comm, ncclFloat);
}

void nccl_bcast_data(std::complex<float>* object, const int& n, MPI_Comm& comm)
{
    nccl_bcast_impl(reinterpret_cast<float*>(object), n, comm, ncclFloat, 2);
}

void nccl_reduce_data(double* object, const int& n, MPI_Comm& comm)
{
    nccl_reduce_impl(object, n, comm, ncclDouble);
}

void nccl_reduce_data(std::complex<double>* object, const int& n, MPI_Comm& comm)
{
    nccl_reduce_impl(reinterpret_cast<double*>(object), n, comm, ncclDouble, 2);
}

void nccl_reduce_data(float* object, const int& n, MPI_Comm& comm)
{
    nccl_reduce_impl(object, n, comm, ncclFloat);
}

void nccl_reduce_data(std::complex<float>* object, const int& n, MPI_Comm& comm)
{
    nccl_reduce_impl(reinterpret_cast<float*>(object), n, comm, ncclFloat, 2);
}

void nccl_gatherv_data(const double* sendbuf, int sendcount, double* recvbuf, const int* recvcounts, const int* displs, MPI_Comm& comm)
{
    nccl_gatherv_impl(sendbuf, sendcount, recvbuf, recvcounts, displs, comm);
}

void nccl_gatherv_data(const std::complex<double>* sendbuf,
                       int sendcount,
                       std::complex<double>* recvbuf,
                       const int* recvcounts,
                       const int* displs,
                       MPI_Comm& comm)
{
    nccl_gatherv_impl(sendbuf, sendcount, recvbuf, recvcounts, displs, comm);
}

void nccl_gatherv_data(const float* sendbuf, int sendcount, float* recvbuf, const int* recvcounts, const int* displs, MPI_Comm& comm)
{
    nccl_gatherv_impl(sendbuf, sendcount, recvbuf, recvcounts, displs, comm);
}

void nccl_gatherv_data(const std::complex<float>* sendbuf,
                       int sendcount,
                       std::complex<float>* recvbuf,
                       const int* recvcounts,
                       const int* displs,
                       MPI_Comm& comm)
{
    nccl_gatherv_impl(sendbuf, sendcount, recvbuf, recvcounts, displs, comm);
}
#endif

void isend_data(const double* buf, int count, int dest, int tag, MPI_Comm& comm, MPI_Request* request)
{
    MPI_Isend(buf, count, MPI_DOUBLE, dest, tag, comm, request);
}
void isend_data(const std::complex<double>* buf, int count, int dest, int tag, MPI_Comm& comm, MPI_Request* request)
{
    MPI_Isend(buf, count, MPI_DOUBLE_COMPLEX, dest, tag, comm, request);
}
void isend_data(const float* buf, int count, int dest, int tag, MPI_Comm& comm, MPI_Request* request)
{
    MPI_Isend(buf, count, MPI_FLOAT, dest, tag, comm, request);
}
void isend_data(const std::complex<float>* buf, int count, int dest, int tag, MPI_Comm& comm, MPI_Request* request)
{
    MPI_Isend(buf, count, MPI_COMPLEX, dest, tag, comm, request);
}
void send_data(const double* buf, int count, int dest, int tag, MPI_Comm& comm)
{
    MPI_Send(buf, count, MPI_DOUBLE, dest, tag, comm);
}
void send_data(const std::complex<double>* buf, int count, int dest, int tag, MPI_Comm& comm)
{
    MPI_Send(buf, count, MPI_DOUBLE_COMPLEX, dest, tag, comm);
}
void send_data(const float* buf, int count, int dest, int tag, MPI_Comm& comm)
{
    MPI_Send(buf, count, MPI_FLOAT, dest, tag, comm);
}
void send_data(const std::complex<float>* buf, int count, int dest, int tag, MPI_Comm& comm)
{
    MPI_Send(buf, count, MPI_COMPLEX, dest, tag, comm);
}
void recv_data(double* buf, int count, int source, int tag, MPI_Comm& comm, MPI_Status* status)
{
    MPI_Recv(buf, count, MPI_DOUBLE, source, tag, comm, status);
}
void recv_data(std::complex<double>* buf, int count, int source, int tag, MPI_Comm& comm, MPI_Status* status)
{
    MPI_Recv(buf, count, MPI_DOUBLE_COMPLEX, source, tag, comm, status);
}
void recv_data(float* buf, int count, int source, int tag, MPI_Comm& comm, MPI_Status* status)
{
    MPI_Recv(buf, count, MPI_FLOAT, source, tag, comm, status);
}
void recv_data(std::complex<float>* buf, int count, int source, int tag, MPI_Comm& comm, MPI_Status* status)
{
    MPI_Recv(buf, count, MPI_COMPLEX, source, tag, comm, status);
}
void bcast_data(std::complex<double>* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n * 2, MPI_DOUBLE, 0, comm);
}
void bcast_data(std::complex<float>* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n * 2, MPI_FLOAT, 0, comm);
}
void bcast_data(double* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n, MPI_DOUBLE, 0, comm);
}
void bcast_data(float* object, const int& n, const MPI_Comm& comm)
{
    MPI_Bcast(object, n, MPI_FLOAT, 0, comm);
}
void reduce_data(std::complex<double>* object, const int& n, const MPI_Comm& comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n * 2, MPI_DOUBLE, MPI_SUM, comm);
}
void reduce_data(std::complex<float>* object, const int& n, const MPI_Comm& comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n * 2, MPI_FLOAT, MPI_SUM, comm);
}
void reduce_data(double* object, const int& n, const MPI_Comm& comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_DOUBLE, MPI_SUM, comm);
}
void reduce_data(float* object, const int& n, const MPI_Comm& comm)
{
    MPI_Allreduce(MPI_IN_PLACE, object, n, MPI_FLOAT, MPI_SUM, comm);
}
void gatherv_data(const double* sendbuf, int sendcount, double* recvbuf, const int* recvcounts, const int* displs, MPI_Comm& comm)
{
    MPI_Allgatherv(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcounts, displs, MPI_DOUBLE, comm);
}
void gatherv_data(const std::complex<double>* sendbuf, int sendcount, std::complex<double>* recvbuf, const int* recvcounts, const int* displs, MPI_Comm& comm)
{
    MPI_Allgatherv(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, recvbuf, recvcounts, displs, MPI_DOUBLE_COMPLEX, comm);
}
void gatherv_data(const float* sendbuf, int sendcount, float* recvbuf, const int* recvcounts, const int* displs, MPI_Comm& comm)
{
    MPI_Allgatherv(sendbuf, sendcount, MPI_FLOAT, recvbuf, recvcounts, displs, MPI_FLOAT, comm);
}
void gatherv_data(const std::complex<float>* sendbuf, int sendcount, std::complex<float>* recvbuf, const int* recvcounts, const int* displs, MPI_Comm& comm)
{
    MPI_Allgatherv(sendbuf, sendcount, MPI_COMPLEX, recvbuf, recvcounts, displs, MPI_COMPLEX, comm);
}

#ifndef __CUDA_MPI
template <typename T>
struct object_cpu_point<T, base_device::DEVICE_GPU>
{
    bool alloc = false;
    T* get_buffer(const T* object, const int& n, T* tmp_space = nullptr)
    {
        T* object_cpu = nullptr;
        alloc = false;

        if (tmp_space == nullptr)
        {
            base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(object_cpu, n);
            alloc = true;
        }
        else
        {
            object_cpu = tmp_space;
        }
        return object_cpu;
    }
    T* get(const T* object, const int& n, T* tmp_space = nullptr)
    {
        T* object_cpu = get_buffer(object, n, tmp_space);
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(object_cpu,
                                                                                                          object,
                                                                                                          n);

        return object_cpu;
    }
    void sync_h2d(T* object, const T* object_cpu, const int& n)
    {
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(object,
                                                                                                          object_cpu,
                                                                                                          n);
    }
    void sync_d2h(T* object_cpu, const T* object, const int& n)
    {
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(object_cpu,
                                                                                                          object,
                                                                                                          n);
    }
    void del(T* object_cpu)
    {
        if (alloc)
        {
            base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(object_cpu);
        }
    }
};

template <typename T>
struct object_cpu_point<T, base_device::DEVICE_CPU>
{
    bool alloc = false;
    T* get_buffer(const T* object, const int& n, T* tmp_space = nullptr)
    {
        return const_cast<T*>(object);
    }
    T* get(const T* object, const int& n, T* tmp_space = nullptr)
    {
        return const_cast<T*>(object);
    }
    void sync_h2d(T* object, const T* object_cpu, const int& n)
    {
    }
    void sync_d2h(T* object_cpu, const T* object, const int& n)
    {
    }
    void del(T* object_cpu)
    {
    }
};

template struct object_cpu_point<double, base_device::DEVICE_CPU>;
template struct object_cpu_point<double, base_device::DEVICE_GPU>;
template struct object_cpu_point<std::complex<double>, base_device::DEVICE_CPU>;
template struct object_cpu_point<std::complex<double>, base_device::DEVICE_GPU>;
template struct object_cpu_point<float, base_device::DEVICE_CPU>;
template struct object_cpu_point<float, base_device::DEVICE_GPU>;
template struct object_cpu_point<std::complex<float>, base_device::DEVICE_CPU>;
template struct object_cpu_point<std::complex<float>, base_device::DEVICE_GPU>;
#endif

} // namespace Parallel_Common
#endif
