#include "device.h"

#include "source_base/tool_quit.h"

#include <base/macros/macros.h>
#include <cstring>
#include <iostream>
#ifdef __MPI
#include "mpi.h"
#endif

#if defined(__CUDA)
#include <cuda_runtime.h>
#endif

#if defined(__ROCM)
#include <hip/hip_runtime.h>
#endif

namespace base_device {

// for device
template <>
base_device::AbacusDevice_t
get_device_type<base_device::DEVICE_CPU>(const base_device::DEVICE_CPU *dev) {
  return base_device::CpuDevice;
}
template <>
base_device::AbacusDevice_t
get_device_type<base_device::DEVICE_GPU>(const base_device::DEVICE_GPU *dev) {
  return base_device::GpuDevice;
}

// for precision
template <> std::string get_current_precision(const float *var) {
  return "single";
}
template <> std::string get_current_precision(const double *var) {
  return "double";
}
template <> std::string get_current_precision(const std::complex<float> *var) {
  return "single";
}
template <> std::string get_current_precision(const std::complex<double> *var) {
  return "double";
}

namespace information {

#if __MPI
int stringCmp(const void *a, const void *b) {
  char *m = (char *)a;
  char *n = (char *)b;
  int i, sum = 0;

  for (i = 0; i < MPI_MAX_PROCESSOR_NAME; i++) {
    if (m[i] == n[i]) {
      continue;
    } else {
      sum = m[i] - n[i];
      break;
    }
  }
  return sum;
}
int get_node_rank() {
  char host_name[MPI_MAX_PROCESSOR_NAME];
  memset(host_name, '\0', sizeof(char) * MPI_MAX_PROCESSOR_NAME);
  char(*host_names)[MPI_MAX_PROCESSOR_NAME];
  int n, namelen, color, rank, nprocs, myrank;
  size_t bytes;
  MPI_Comm nodeComm;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Get_processor_name(host_name, &namelen);

  bytes = nprocs * sizeof(char[MPI_MAX_PROCESSOR_NAME]);
  host_names = (char(*)[MPI_MAX_PROCESSOR_NAME])malloc(bytes);
  for (int ii = 0; ii < nprocs; ii++) {
    memset(host_names[ii], '\0', sizeof(char) * MPI_MAX_PROCESSOR_NAME);
  }

  strcpy(host_names[rank], host_name);

  for (n = 0; n < nprocs; n++) {
    MPI_Bcast(&(host_names[n]), MPI_MAX_PROCESSOR_NAME, MPI_CHAR, n,
              MPI_COMM_WORLD);
  }
  qsort(host_names, nprocs, sizeof(char[MPI_MAX_PROCESSOR_NAME]), stringCmp);

  color = 0;
  for (n = 0; n < nprocs - 1; n++) {
    if (strcmp(host_name, host_names[n]) == 0) {
      break;
    }
    if (strcmp(host_names[n], host_names[n + 1])) {
      color++;
    }
  }

  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &nodeComm);
  MPI_Comm_rank(nodeComm, &myrank);

  MPI_Barrier(MPI_COMM_WORLD);
  int looprank = myrank;
  // printf (" Assigning device %d  to process on node %s rank %d,
  // OK\n",looprank,  host_name, rank );
  free(host_names);
  return looprank;
}

int get_node_rank_with_mpi_shared(const MPI_Comm mpi_comm) {
  // 20240530 zhanghaochong
  // The main difference between this function and the above is that it does not
  // use hostname, but uses MPI's built-in function to achieve similar
  // functions.
  MPI_Comm localComm;
  int localMpiRank;
  MPI_Comm_split_type(mpi_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &localComm);
  MPI_Comm_rank(localComm, &localMpiRank);
  MPI_Comm_free(&localComm);
  return localMpiRank;
}
#if defined(__CUDA)

int set_device_by_rank(const MPI_Comm mpi_comm) {
  int localMpiRank = get_node_rank_with_mpi_shared(mpi_comm);
  int device_num = -1;

  cudaGetDeviceCount(&device_num);
  if (device_num <= 0) {
    ModuleBase::WARNING_QUIT("device", "can not find gpu device!");
  }
  // warning: this is not a good way to assign devices, user should assign One
  // process per GPU
  int local_device_id = localMpiRank % device_num;
  int ret = cudaSetDevice(local_device_id);
  if (ret != cudaSuccess) {
    ModuleBase::WARNING_QUIT("device", "cudaSetDevice failed!");
  }
  return local_device_id;
}
#endif

#endif

bool probe_gpu_availability() {
#if defined(__CUDA)
    int device_count = 0;
    // Directly call cudaGetDeviceCount without cudaErrcheck to prevent program exit
    cudaError_t error_id = cudaGetDeviceCount(&device_count);
    if (error_id == cudaSuccess && device_count > 0) {
        return true;
    }
    return false;
#elif defined(__ROCM)
    int device_count = 0;
    hipError_t error_id = hipGetDeviceCount(&device_count);
    if (error_id == hipSuccess && device_count > 0) {
        return true;
    }
    return false;
#else
    // If not compiled with GPU support, GPU is not available
    return false;
#endif
}

std::string get_device_flag(const std::string &device,
                            const std::string &basis_type) {
    // 1. Validate input string
    if (device != "cpu" && device != "gpu" && device != "auto") {
        ModuleBase::WARNING_QUIT("device", "Parameter \"device\" can only be set to \"cpu\", \"gpu\", or \"auto\"!");
    }
    
    // NOTE: This function is called only on rank 0 during input parsing.
    // The result will be broadcast to other ranks via the standard bcast mechanism.
    // DO NOT use MPI_Bcast here as other ranks are not in this code path.
    
    std::string result = "cpu";
    
    if (device == "gpu") {
        if (probe_gpu_availability()) {
            result = "gpu";
            // std::cout << " INFO: 'device=gpu' specified. GPU will be used." << std::endl;
        } else {
            ModuleBase::WARNING_QUIT("device", "Device is set to 'gpu', but no available GPU was found. Please check your hardware/drivers or set 'device=cpu'.");
        }
    } else if (device == "auto") {
        if (probe_gpu_availability()) {
            result = "gpu";
            // std::cout << " INFO: 'device=auto' specified. GPU detected and will be used." << std::endl;
        } else {
            result = "cpu";
            // std::cout << " WARNING: 'device=auto' specified, but no GPU was found. Falling back to CPU." << std::endl;
            // std::cout << "          To suppress this warning, please explicitly set 'device=cpu' in your input." << std::endl;
        }
    } else { // device == "cpu"
        result = "cpu";
        // std::cout << " INFO: 'device=cpu' specified. CPU will be used." << std::endl;
    }

    // 2. Final check for incompatible basis type
    if (result == "gpu" && basis_type == "lcao_in_pw") {
        ModuleBase::WARNING_QUIT("device", "The GPU currently does not support the basis type \"lcao_in_pw\"!");
    }

    // 3. Return the final decision
    return result;
}

int get_device_kpar(const int& kpar, const int& bndpar)
{
#if __MPI && (__CUDA || __ROCM)
    // This function should only be called when device mode is GPU
    // The device decision has already been made by get_device_flag()
    int temp_nproc = 0;
    int new_kpar = kpar;
    MPI_Comm_size(MPI_COMM_WORLD, &temp_nproc);
    if (temp_nproc != kpar * bndpar)
    {
        new_kpar = temp_nproc / bndpar;
        ModuleBase::WARNING("Input_conv", "kpar is not compatible with the number of processors, auto set kpar value.");
    }
    
    // get the CPU rank of current node
    int node_rank = base_device::information::get_node_rank();

    int device_num = -1;
#if defined(__CUDA)
    cudaErrcheck(cudaGetDeviceCount(&device_num)); // get the number of GPU devices of current node
    cudaErrcheck(cudaSetDevice(node_rank % device_num)); // bind the CPU processor to the devices
#elif defined(__ROCM)
    hipErrcheck(hipGetDeviceCount(&device_num));
    hipErrcheck(hipSetDevice(node_rank % device_num));
#endif
    return new_kpar;
#endif
    return kpar;
}

} // end of namespace information
} // end of namespace base_device
