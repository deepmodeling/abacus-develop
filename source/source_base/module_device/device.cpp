
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

std::string get_device_flag(const std::string &device,
                            const std::string &basis_type) {
if (device == "cpu") {
  return "cpu"; // no extra checks required
}
std::string error_message;
if (device != "auto" and device != "gpu")
{
  error_message += "Parameter \"device\" can only be set to \"cpu\" or \"gpu\"!";
  ModuleBase::WARNING_QUIT("device", error_message);
}

// Get available GPU count
int device_count = -1;
#if ((defined __CUDA) || (defined __ROCM))
#if defined(__CUDA)
cudaGetDeviceCount(&device_count);
#elif defined(__ROCM)
hipGetDeviceCount(&device_count);
/***auto start_time = std::chrono::high_resolution_clock::now();
std::cout << "Starting hipGetDeviceCount.." << std::endl;
auto end_time = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
std::cout << "hipGetDeviceCount took " << duration.count() << "seconds" << std::endl;***/
#endif
if (device_count <= 0)
{
  error_message += "Cannot find GPU on this computer!\n";
}
#else // CPU only
error_message += "ABACUS is built with CPU support only. Please rebuild with GPU support.\n";
#endif

if (basis_type == "lcao_in_pw") {
  error_message +=
      "The GPU currently does not support the basis type \"lcao_in_pw\"!";
}
if(error_message.empty())
{
  return "gpu"; // possibly automatically set to GPU
}
else if (device == "gpu")
{
  ModuleBase::WARNING_QUIT("device", error_message);
}
else { return "cpu";
}
}

int get_device_kpar(const int& kpar, const int& bndpar)
{
#if __MPI && (__CUDA || __ROCM)
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
  cudaGetDeviceCount(&device_num); // get the number of GPU devices of current node
  cudaSetDevice(node_rank % device_num); // band the CPU processor to the devices
#elif defined(__ROCM)
  hipGetDeviceCount(&device_num);
  hipSetDevice(node_rank % device_num);
#endif
  return new_kpar;
#endif
  return kpar;
}

} // end of namespace information
} // end of namespace base_device
