#ifndef __PARALLEL_DEVICE_H__
#define __PARALLEL_DEVICE_H__
#ifdef __MPI
#include "mpi.h"
#include "module_base/module_device/device.h"
#include "module_base/module_device/memory_op.h"
#include <complex>
namespace Parallel_Common
{
void isend_data(const double* buf, int count, int dest, int tag, MPI_Comm& comm, MPI_Request* request);
void isend_data(const std::complex<double>* buf, int count, int dest, int tag, MPI_Comm& comm, MPI_Request* request);
void isend_data(const float* buf, int count, int dest, int tag, MPI_Comm& comm, MPI_Request* request);
void isend_data(const std::complex<float>* buf, int count, int dest, int tag, MPI_Comm& comm, MPI_Request* request);
void send_data(const double* buf, int count, int dest, int tag, MPI_Comm& comm);
void send_data(const std::complex<double>* buf, int count, int dest, int tag, MPI_Comm& comm);
void send_data(const float* buf, int count, int dest, int tag, MPI_Comm& comm);
void send_data(const std::complex<float>* buf, int count, int dest, int tag, MPI_Comm& comm);
void recv_data(double* buf, int count, int source, int tag, MPI_Comm& comm, MPI_Status* status);
void recv_data(std::complex<double>* buf, int count, int source, int tag, MPI_Comm& comm, MPI_Status* status);
void recv_data(float* buf, int count, int source, int tag, MPI_Comm& comm, MPI_Status* status);
void recv_data(std::complex<float>* buf, int count, int source, int tag, MPI_Comm& comm, MPI_Status* status);
void bcast_data(std::complex<double>* object, const int& n, const MPI_Comm& comm);
void bcast_data(std::complex<float>* object, const int& n, const MPI_Comm& comm);
void bcast_data(double* object, const int& n, const MPI_Comm& comm);
void bcast_data(float* object, const int& n, const MPI_Comm& comm);
void reduce_data(std::complex<double>* object, const int& n, const MPI_Comm& comm);
void reduce_data(std::complex<float>* object, const int& n, const MPI_Comm& comm);
void reduce_data(double* object, const int& n, const MPI_Comm& comm);
void reduce_data(float* object, const int& n, const MPI_Comm& comm);
void gatherv_data(const double* sendbuf, int sendcount, double* recvbuf, const int* recvcounts, const int* displs, MPI_Comm& comm);
void gatherv_data(const std::complex<double>* sendbuf, int sendcount, std::complex<double>* recvbuf, const int* recvcounts, const int* displs, MPI_Comm& comm);
void gatherv_data(const float* sendbuf, int sendcount, float* recvbuf, const int* recvcounts, const int* displs, MPI_Comm& comm);
void gatherv_data(const std::complex<float>* sendbuf, int sendcount, std::complex<float>* recvbuf, const int* recvcounts, const int* displs, MPI_Comm& comm);

template<typename T, typename Device>
struct object_cpu_point
{
    bool alloc = false;
    T* get(const T* object, const int& n, T* tmp_space = nullptr);
    void del(T* object);
    void sync_d2h(T* object_cpu, const T* object, const int& n);
    void sync_h2d(T* object, const T* object_cpu, const int& n);
};

/**
 * @brief send data in Device
 * 
 */
template <typename T, typename Device>
void send_dev(const T* object, int count, int dest, int tag, MPI_Comm& comm, T* tmp_space = nullptr)
{
    object_cpu_point<T,Device> o;
    T* object_cpu = o.get(object, count, tmp_space);
    o.sync_d2h(object_cpu, object, count);
    send_data(object_cpu, count, dest, tag, comm);
    o.del(object_cpu);
    return;
}

/**
 * @brief isend data in Device
 * @note before the date in send_space is recieved, it should not be modified
 * 
 */
template <typename T, typename Device>
void isend_dev(const T* object, int count, int dest, int tag, MPI_Comm& comm, MPI_Request* request, T* send_space)
{
    object_cpu_point<T,Device> o;
    T* object_cpu = o.get(object, count, send_space);
    o.sync_d2h(object_cpu, object, count);
    isend_data(object_cpu, count, dest, tag, comm, request);
    o.del(object_cpu);
    return;
}

/**
 * @brief recv data in Device
 * 
 */
template <typename T, typename Device>
void recv_dev(T* object, int count, int source, int tag, MPI_Comm& comm, MPI_Status* status, T* tmp_space = nullptr)
{
    object_cpu_point<T,Device> o;
    T* object_cpu = o.get(object, count, tmp_space);
    recv_data(object_cpu, count, source, tag, comm, status);
    o.sync_h2d(object, object_cpu, count);
    o.del(object_cpu);
    return;
}

/**
 * @brief bcast data in Device
 * 
 * @tparam T: float, double, std::complex<float>, std::complex<double>
 * @tparam Device 
 * @param ctx Device ctx
 * @param object complex arrays in Device
 * @param n the size of complex arrays
 * @param comm MPI_Comm
 * @param tmp_space tmp space in CPU
 */
template <typename T, typename Device>
void bcast_dev(T* object, const int& n, const MPI_Comm& comm, T* tmp_space = nullptr)
{
    object_cpu_point<T,Device> o;
    T* object_cpu = o.get(object, n, tmp_space);
    o.sync_d2h(object_cpu, object, n);
    bcast_data(object_cpu, n, comm);
    o.sync_h2d(object, object_cpu, n);
    o.del(object_cpu);
    return;
}

template <typename T, typename Device>
void reduce_dev(T* object, const int& n, const MPI_Comm& comm, T* tmp_space = nullptr)
{
    object_cpu_point<T,Device> o;
    T* object_cpu = o.get(object, n, tmp_space);
    o.sync_d2h(object_cpu, object, n);
    reduce_data(object_cpu, n, comm);
    o.sync_h2d(object, object_cpu, n);
    o.del(object_cpu);
    return;
}
}
    

#endif
#endif