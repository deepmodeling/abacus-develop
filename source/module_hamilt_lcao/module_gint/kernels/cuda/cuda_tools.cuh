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

/*
 * @brief: A simple wrapper for cudaMalloc and cudaFree, sync and async CUDA memory copy
 * @param: T: the type of the data
 * 
 * @note: 
 * Manual management of CUDA memory is a very delicate task; complex pointers and 
 * malloc/free operations make it easy for us to encounter memory bugs. 
 * The severity of the issues increases significantly when introducing multi-node, 
 * multi-GPU, and multi-stream parallelism. 
 * Debugging after encountering bugs is also very difficult, finding the leaking 
 * pointer from dozens of variables can be quite a headache. 
 * Therefore, considering that our use and management of memory have some 
 * homogeneity, we have abstracted these needs into the following encapsulations 
 * to reduce the cost of maintenance and development. The memory is allocated in 
 * the constructor and freed in the destructor.
 * 
 * The following interface is primarily designed for the following requirements: 
 * 1. We need to split a large task into multiple subtasks to run on multiple 
 *    streams across multiple GPUs on multiple nodes. 
 * 2. It is necessary to allocate memory of the same shape on both host and device. 
 * 3. Data copying between host and device sync or async is required.
 */

template <typename T>
class Cuda_Mem_Wrapper
{
public:
    Cuda_Mem_Wrapper(int one_stream_size, int one_stream_size_aligned, int stream_number = 1, bool malloc_host = true);
    Cuda_Mem_Wrapper(int one_stream_size, int stream_number = 1, bool malloc_host = true);
    ~Cuda_Mem_Wrapper();
    void copy_host_to_device_sync(int stream_id = 0);
    void copy_host_to_device_async(cudaStream_t stream, int stream_id);
    void copy_device_to_host_sync(int stream_id = 0);
    void copy_device_to_host_async(cudaStream_t stream, int stream_id);
    void memset_device_sync(int stream_id = 0, int value = 0);
    void memset_device_async(cudaStream_t stream, int stream_id = 0, int value = 0);
    void memset_host(int stream_id = 0, int value = 0);
    T *get_device_pointer(int stream_id = 0);
    T *get_host_pointer(int stream_id = 0);
    void free_all();
private:
    T *device_pointer;
    T *host_pointer;
    int one_stream_size;
    int one_stream_size_aligned;
    int stream_number;
    int total_size_aligned;
};

#endif // CUDA_TOOLS_CUH#ifndef CUDA_TOOLS_CUH