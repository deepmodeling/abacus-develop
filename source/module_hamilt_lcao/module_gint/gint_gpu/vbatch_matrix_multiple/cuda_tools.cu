#include "module_hamilt_lcao/module_gint/gint_gpu/vbatch_matrix_multiple/cuda_tools.cuh"
cudaError_t checkCuda(cudaError_t result)
{
#if defined(__DEBUG)
    if (result != cudaSuccess)
    {
        fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
#endif
    return result;
}

cudaError_t checkCudaLastError()
{
#if defined(__DEBUG)
    cudaError_t result = cudaGetLastError();
    if (result != cudaSuccess)
    {
        fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
    return result;
#else
    return cudaSuccess;
#endif
}

void dump_cuda_array_to_file(double *cuda_array, int width, int hight, const std::string &filename)
{
    double *h_data = new double[width * hight];
    cudaMemcpy(h_data, cuda_array, width * hight * sizeof(double), cudaMemcpyDeviceToHost);

    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
        std::cerr << "Failed to open file for writing." << std::endl;
    }
    for (int j = 0; j < hight; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            outFile << "hight" << j << "   width:" << i << "   " << h_data[j * width + i] << std::endl;
        }
    }
    outFile.close();
    delete[] h_data;
}

template <typename T>
Cuda_Mem_Wrapper<T>::Cuda_Mem_Wrapper(int one_stream_size, int stream_number, bool malloc_host)
{
    this->stream_number = stream_number;
    this->one_stream_size = one_stream_size;
    this->total_size = stream_number * one_stream_size;
    checkCuda(cudaMalloc((void **)&this->device_pointer, this->total_size * sizeof(T)));
    checkCuda(cudaMemset(this->device_pointer, 0, this->total_size * sizeof(T)));
    this->host_pointer = nullptr;

    if (malloc_host)
    {
        checkCuda(cudaMallocHost((void **)&this->host_pointer, this->total_size * sizeof(T)));
        memset(this->host_pointer, 0, this->total_size * sizeof(T));
    }
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
void Cuda_Mem_Wrapper<T>::copy_host_to_device_sync(int stream_id)
{
    #if defined(__DEBUG)
    if (this->host_pointer == nullptr || this->device_pointer == nullptr)
    {
        std::cerr << "host_pointer is nullptr, can not copy host to device" << std::endl;
        exit(1);
    }
    #endif
    checkCuda(cudaMemcpy(this->device_pointer + stream_id * this->one_stream_size, this->host_pointer + stream_id * this->one_stream_size,
               this->one_stream_size * sizeof(T), cudaMemcpyHostToDevice));
}

template <typename T>
void Cuda_Mem_Wrapper<T>::copy_host_to_device_async(int stream_id, cudaStream_t stream)
{
    #if defined(__DEBUG)
    if (this->host_pointer == nullptr || this->device_pointer == nullptr)
    {
        std::cerr << "host_pointer is nullptr, can not copy host to device" << std::endl;
        exit(1);
    }
    #endif
    checkCuda(cudaMemcpyAsync(this->device_pointer + stream_id * this->one_stream_size,
                              this->host_pointer + stream_id * this->one_stream_size,
                              this->one_stream_size * sizeof(T), cudaMemcpyHostToDevice, stream));
}

template <typename T>
void Cuda_Mem_Wrapper<T>::copy_device_to_host_sync(int stream_id)
{
    #if defined(__DEBUG)
    if (this->host_pointer == nullptr || this->device_pointer == nullptr)
    {
        std::cerr << "host_pointer is nullptr, can not copy device to host" << std::endl;
        exit(1);
    }
    #endif
    checkCuda(cudaMemcpy(this->host_pointer + stream_id * this->one_stream_size, this->device_pointer + stream_id * this->one_stream_size,
               this->one_stream_size * sizeof(T), cudaMemcpyDeviceToHost));
}

template <typename T>
void Cuda_Mem_Wrapper<T>::copy_device_to_host_async(int stream_id, cudaStream_t stream)
{
    #if defined(__DEBUG)
    if (this->host_pointer == nullptr || this->device_pointer == nullptr)
    {
        std::cerr << "host_pointer is nullptr, can not copy device to host" << std::endl;
        exit(1);
    }
    #endif
    checkCuda(cudaMemcpyAsync(this->host_pointer + stream_id * this->one_stream_size,
                              this->device_pointer + stream_id * this->one_stream_size,
                              this->one_stream_size * sizeof(T), cudaMemcpyDeviceToHost, stream));
}

template <typename T>
T* Cuda_Mem_Wrapper<T>::get_device_pointer(int stream_id)
{
    return this->device_pointer + stream_id * this->one_stream_size;
}

template <typename T>
T* Cuda_Mem_Wrapper<T>::get_host_pointer(int stream_id)
{
    return this->host_pointer + stream_id * this->one_stream_size;
}
template class Cuda_Mem_Wrapper<double>;
template class Cuda_Mem_Wrapper<double*>;
template class Cuda_Mem_Wrapper<int>;