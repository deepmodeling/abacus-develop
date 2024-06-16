#include <iostream>
#include <cstring>

#include "cuda_tools.cuh"

cudaError_t check(cudaError_t result, const char *const func, const char *const file, const int line)
{
    if (result != cudaSuccess)
    {
        fprintf(stderr, "CUDA Runtime Error at %s:%d code=%s \"%s\" \n", file, line, cudaGetErrorString(result), func);
        exit(EXIT_FAILURE);
    }
    return result;
}
cudaError_t __checkCudaLastError(const char *file, const int line)
{
    cudaError_t result = cudaGetLastError();
    if (result != cudaSuccess)
    {
        fprintf(stderr, "%s(%i) : getLastCudaError():%s\n", file, line, cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
    return result;
}

void dump_cuda_array_to_file(double* cuda_array,
                             int width,
                             int hight,
                             const std::string& filename)
{
    double* h_data = new double[width * hight];
    cudaMemcpy(h_data,
               cuda_array,
               width * hight * sizeof(double),
               cudaMemcpyDeviceToHost);

    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
        std::cerr << "Failed to open file for writing." << std::endl;
    }
    for (int j = 0; j < hight; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            outFile << "hight" << j << "   width:" << i << "   "
                    << h_data[j * width + i] << std::endl;
        }
    }
    outFile.close();
    delete[] h_data;
}

template <typename T>
Cuda_Mem_Wrapper<T>::Cuda_Mem_Wrapper()
{
    this->device_pointer = nullptr;
    this->host_pointer = nullptr;
    this->one_stream_size = 0;
    this->one_stream_size_aligned = 0;
    this->stream_number = 1;
    this->total_size_aligned = 0;
}

template <typename T>
Cuda_Mem_Wrapper<T>::Cuda_Mem_Wrapper(int one_stream_size_in,
                                      int one_stream_size_aligned_in,
                                      int stream_number_in,
                                      bool malloc_host_in)
{
    this->stream_number = stream_number_in;
    this->one_stream_size = one_stream_size_in;
    this->one_stream_size_aligned = one_stream_size_aligned_in;
    this->total_size_aligned
        = this->one_stream_size_aligned * this->stream_number;

    checkCuda(cudaMalloc((void**)&this->device_pointer,
                         this->total_size_aligned * sizeof(T)));
    checkCuda(cudaMemset(this->device_pointer,
                         0,
                         this->total_size_aligned * sizeof(T)));
    this->host_pointer = nullptr;

    if (malloc_host_in)
    {
        checkCuda(cudaMallocHost((void**)&this->host_pointer,
                                 this->total_size_aligned * sizeof(T)));
        memset(this->host_pointer, 0, this->total_size_aligned * sizeof(T));
    }
}

template <typename T>
Cuda_Mem_Wrapper<T>::Cuda_Mem_Wrapper(int one_stream_size_in,
                                      int stream_number_in,
                                      bool malloc_host_in)
    : Cuda_Mem_Wrapper(one_stream_size_in,
                       one_stream_size_in,
                       stream_number_in,
                       malloc_host_in)
{
}

template <typename T>
Cuda_Mem_Wrapper<T>::Cuda_Mem_Wrapper(Cuda_Mem_Wrapper&& other) noexcept
{
    this->device_pointer = other.device_pointer;
    this->host_pointer = other.host_pointer;
    this->one_stream_size = other.one_stream_size;
    this->one_stream_size_aligned = other.one_stream_size_aligned;
    this->stream_number = other.stream_number;
    this->total_size_aligned = other.total_size_aligned;

    other.device_pointer = nullptr;
    other.host_pointer = nullptr;
    other.one_stream_size = 0;
    other.one_stream_size_aligned = 0;
    other.stream_number = 0;
    other.total_size_aligned = 0;
}

template <typename T>
Cuda_Mem_Wrapper<T>& Cuda_Mem_Wrapper<T>::operator=(Cuda_Mem_Wrapper&& other) noexcept
{
    if (this != &other)
    {
        this->free_all();
        this->device_pointer = other.device_pointer;
        this->host_pointer = other.host_pointer;
        this->one_stream_size = other.one_stream_size;
        this->one_stream_size_aligned = other.one_stream_size_aligned;
        this->stream_number = other.stream_number;
        this->total_size_aligned = other.total_size_aligned;

        other.device_pointer = nullptr;
        other.host_pointer = nullptr;
        other.one_stream_size = 0;
        other.one_stream_size_aligned = 0;
        other.stream_number = 0;
        other.total_size_aligned = 0;
    }
    return *this;
}

template <typename T>
void Cuda_Mem_Wrapper<T>::free_all()
{
    checkCuda(cudaFree(this->device_pointer));
    if (this->host_pointer != nullptr)
    {
        checkCuda(cudaFreeHost(this->host_pointer));
    }
}

template <typename T>
Cuda_Mem_Wrapper<T>::~Cuda_Mem_Wrapper()
{
    this->free_all();
}

template <typename T>
inline void Cuda_Mem_Wrapper<T>::copy_host_to_device_sync(int stream_id)
{
    if (this->host_pointer == nullptr || this->device_pointer == nullptr)
    {
        std::cerr << "host_pointer is nullptr, can not copy host to device"
                  << std::endl;
        exit(1);
    }
    checkCuda(cudaMemcpy(
        this->device_pointer + stream_id * this->one_stream_size_aligned,
        this->host_pointer + stream_id * this->one_stream_size_aligned,
        this->one_stream_size * sizeof(T),
        cudaMemcpyHostToDevice));
}

template <typename T>
inline void Cuda_Mem_Wrapper<T>::copy_host_to_device_async(cudaStream_t stream,
                                                    int stream_id)
{
    if (this->host_pointer == nullptr || this->device_pointer == nullptr)
    {
        std::cerr << "host_pointer is nullptr, can not copy host to device"
                  << std::endl;
        exit(1);
    }
    checkCuda(cudaMemcpyAsync(
        this->device_pointer + stream_id * this->one_stream_size_aligned,
        this->host_pointer + stream_id * this->one_stream_size_aligned,
        this->one_stream_size * sizeof(T),
        cudaMemcpyHostToDevice,
        stream));
}

template <typename T>
inline void Cuda_Mem_Wrapper<T>::copy_device_to_host_sync(int stream_id)
{
    if (this->host_pointer == nullptr || this->device_pointer == nullptr)
    {
        std::cerr << "host_pointer is nullptr, can not copy device to host"
                  << std::endl;
        exit(1);
    }
    checkCuda(cudaMemcpy(
        this->host_pointer + stream_id * this->one_stream_size_aligned,
        this->device_pointer + stream_id * this->one_stream_size_aligned,
        this->one_stream_size * sizeof(T),
        cudaMemcpyDeviceToHost));
}

template <typename T>
inline void Cuda_Mem_Wrapper<T>::copy_device_to_host_async(cudaStream_t stream,
                                                    int stream_id)
{
    if (this->host_pointer == nullptr || this->device_pointer == nullptr)
    {
        std::cerr << "host_pointer is nullptr, can not copy device to host"
                  << std::endl;
        exit(1);
    }
    checkCuda(cudaMemcpyAsync(
        this->host_pointer + stream_id * this->one_stream_size_aligned,
        this->device_pointer + stream_id * this->one_stream_size_aligned,
        this->one_stream_size * sizeof(T),
        cudaMemcpyDeviceToHost,
        stream));
}

template <typename T>
inline void Cuda_Mem_Wrapper<T>::memset_device_sync(int stream_id, int value)
{
    checkCuda(cudaMemset(this->device_pointer
                             + stream_id * this->one_stream_size_aligned,
                         value,
                         this->one_stream_size * sizeof(T)));
}

template <typename T>
inline void Cuda_Mem_Wrapper<T>::memset_device_async(cudaStream_t stream,
                                              int stream_id,
                                              int value)
{
    checkCuda(cudaMemsetAsync(this->device_pointer
                                  + stream_id * this->one_stream_size_aligned,
                              value,
                              this->one_stream_size * sizeof(T),
                              stream));
}

template <typename T>
inline void Cuda_Mem_Wrapper<T>::memset_host(int stream_id, int value)
{
    memset(this->host_pointer + stream_id * this->one_stream_size_aligned,
           value,
           this->one_stream_size * sizeof(T));
}

template <typename T>
inline T* Cuda_Mem_Wrapper<T>::get_device_pointer(int stream_id)
{
    return this->device_pointer + stream_id * this->one_stream_size_aligned;
}

template <typename T>
inline T* Cuda_Mem_Wrapper<T>::get_host_pointer(int stream_id)
{
    return this->host_pointer + stream_id * this->one_stream_size_aligned;
}
template class Cuda_Mem_Wrapper<double>;
template class Cuda_Mem_Wrapper<double*>;
template class Cuda_Mem_Wrapper<int>;
