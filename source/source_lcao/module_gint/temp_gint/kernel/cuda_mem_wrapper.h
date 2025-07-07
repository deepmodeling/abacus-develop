#pragma once
#include <cuda_runtime.h>
#include "source_base/tool_quit.h"
#include "gint_helper.cuh"

template <typename T>
class CudaMemWrapper
{
  public:

    CudaMemWrapper() = default;
    CudaMemWrapper(const CudaMemWrapper& other) = delete;
    CudaMemWrapper& operator=(const CudaMemWrapper& other) = delete;
    CudaMemWrapper(CudaMemWrapper&& other) noexcept
    {
      this->device_ptr_ = other.device_ptr_;
      this->host_ptr_ = other.host_ptr_;
      this->size_ = other.size_;
      this->malloc_host_ = other.malloc_host_;
      this->stream_ = other.stream_;

      other.device_ptr_ = nullptr;
      other.host_ptr_ = nullptr;
      other.size_ = 0;
      other.malloc_host_ = false;
      other.stream_ = 0;
    }
    
    CudaMemWrapper& operator=(CudaMemWrapper&& other) noexcept
    {
      if (this != &other)
      {
        this->device_ptr_ = other.device_ptr_;
        this->host_ptr_ = other.host_ptr_;
        this->size_ = other.size_;
        this->malloc_host_ = other.malloc_host_;
        this->stream_ = other.stream_;

        other.device_ptr_ = nullptr;
        other.host_ptr_ = nullptr;
        other.size_ = 0;
        other.malloc_host_ = false;
        other.stream_ = 0;
      }
      return *this;
    }

    CudaMemWrapper(size_t size,
                   cudaStream_t stream = 0,
                   bool malloc_host = true)
    {
      size_ = size;
      malloc_host_ = malloc_host;
      stream_ = stream;

      if (malloc_host)
      { 
        checkCuda(cudaMallocHost((void**)&host_ptr_, size_* sizeof(T)));
        memset(host_ptr_, 0, size_ * sizeof(T));
      }
      else
      { host_ptr_ = nullptr; }

      checkCuda(cudaMalloc((void**)&device_ptr_, size_ * sizeof(T)));
      checkCuda(cudaMemset(device_ptr_, 0, size_ * sizeof(T)));
    }
    
    ~CudaMemWrapper()
    {
      free();
    }

    void copy_host_to_device_sync(size_t size)
    {
      if (host_ptr_ == nullptr)
        { ModuleBase::WARNING_QUIT("cuda_mem_wrapper", "Host pointer is null, cannot copy to device."); } 
      checkCuda(cudaMemcpy(device_ptr_, host_ptr_, size * sizeof(T), cudaMemcpyHostToDevice));
    }

    void copy_host_to_device_sync()
    {
      copy_host_to_device_sync(size_);
    }

    void copy_host_to_device_async(size_t size)
    {
      if (host_ptr_ == nullptr)
        { ModuleBase::WARNING_QUIT("cuda_mem_wrapper", "Host pointer is null, cannot copy to device."); } 
      checkCuda(cudaMemcpyAsync(device_ptr_, host_ptr_, size * sizeof(T), cudaMemcpyHostToDevice, stream_));
    }

    void copy_host_to_device_async()
    {
      copy_host_to_device_async(size_);
    }

    void copy_device_to_host_sync(size_t size)
    {
      if (host_ptr_ == nullptr)
        { ModuleBase::WARNING_QUIT("cuda_mem_wrapper", "Host pointer is null, cannot copy to host."); } 
      checkCuda(cudaMemcpy(host_ptr_, device_ptr_, size * sizeof(T), cudaMemcpyDeviceToHost));
    }

    void copy_device_to_host_sync()
    {
      copy_device_to_host_sync(size_);
    }

    void copy_device_to_host_async(size_t size)
    {
      if (host_ptr_ == nullptr)
        { ModuleBase::WARNING_QUIT("cuda_mem_wrapper", "Host pointer is null, cannot copy to host."); } 
      checkCuda(cudaMemcpyAsync(host_ptr_, device_ptr_, size * sizeof(T), cudaMemcpyDeviceToHost, stream_));
    }

    void copy_device_to_host_async()
    {
      copy_device_to_host_async(size_);
    }
    
    void memset_device_sync(const size_t size, const int value = 0)
    {
      checkCuda(cudaMemset(device_ptr_, value, size * sizeof(T)));
    }

    void memset_device_sync(const int value = 0)
    {
      memset_device_sync(size_, value);
    }

    void memset_device_async(const size_t size, const int value = 0)
    {
      checkCuda(cudaMemsetAsync(device_ptr_, value, size * sizeof(T), stream_));
    }

    void memset_device_async(const int value = 0)
    {
      memset_device_async(size_, value);
    }

    void memset_host(const size_t size, const int value = 0)
    {
      if (host_ptr_ == nullptr)
        { ModuleBase::WARNING_QUIT("cuda_mem_wrapper", "Host pointer is null, cannot memset host."); } 
      checkCuda(cudaMemset(host_ptr_, value, size * sizeof(T)));
    }

    void memset_host(const int value = 0)
    {
      memset_host(size_, value);
    }

    void free()
    {
      checkCuda(cudaFree(device_ptr_));
      checkCuda(cudaFreeHost(host_ptr_));
    }

    T* get_device_ptr() { return device_ptr_; }
    T* get_host_ptr() { return host_ptr_; }
    const T* get_device_ptr() const { return device_ptr_; }
    const T* get_host_ptr() const { return host_ptr_; }
    size_t get_size() const { return size_; }

  private:
    T* device_ptr_ = nullptr;
    T* host_ptr_ = nullptr;
    size_t size_ = 0;
    bool malloc_host_ = false;
    cudaStream_t stream_ = 0;
};