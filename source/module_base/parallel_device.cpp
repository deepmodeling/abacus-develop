#include "parallel_device.h"
#ifdef __MPI
namespace Parallel_Common
{
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

template <typename T>
struct object_cpu_point<T, base_device::DEVICE_GPU>
{
    bool alloc = false;
    T* get(const T* object, const int& n, T* tmp_space = nullptr)
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

} // namespace Parallel_Common
#endif