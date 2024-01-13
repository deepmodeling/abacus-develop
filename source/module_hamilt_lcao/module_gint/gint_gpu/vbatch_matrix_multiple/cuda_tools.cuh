#ifndef CUDA_TOOLS_CUH
#define CUDA_TOOLS_CUH
#include <cuda.h>    // for CUDA_VERSION
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <assert.h>   // for assert

cudaError_t checkCuda(cudaError_t result);
cudaError_t checkCudaLastError();

void dump_cuda_array_to_file(double *cuda_array, int width, int hight, const std::string &filename);

template <typename T>
class Cuda_Mem_Wrapper
{
public:
    Cuda_Mem_Wrapper(int one_stream_size, int stream_number = 1, bool malloc_host = true);
    ~Cuda_Mem_Wrapper();
    void copy_host_to_device_sync(int stream_id = 0);
    void copy_host_to_device_async(int stream_id, cudaStream_t stream);
    void copy_device_to_host_sync(int stream_id = 0);
    void copy_device_to_host_async(int stream_id, cudaStream_t stream);
    T *get_device_pointer(int stream_id = 0);
    T *get_host_pointer(int stream_id = 0);
    void free_all();
private:
    T *device_pointer;
    T *host_pointer;
    int one_stream_size;
    int stream_number;
    int total_size;
};

#endif // CUDA_TOOLS_CUH#ifndef CUDA_TOOLS_CUH
