
#include <iostream>
#include <cstring>
#include "module_psi/kernels/types.h"
#include "module_psi/kernels/device.h"
#include "module_base/tool_quit.h"

#if defined(__CUDA)
#include <cuda_runtime.h>
#endif

#if defined(__ROCM)
#include <hip/hip_runtime.h>
#endif

#ifdef __MPI
#include "mpi.h"
#endif

namespace psi{

namespace device{

static bool is_init = false;

// functions used in custom ops
template<> AbacusDevice_t get_device_type <DEVICE_CPU> (const DEVICE_CPU* dev) {
    return CpuDevice;
}
template<> std::string get_current_precision(const double* var) {
    return "double";
}
template<> std::string get_current_precision (const std::complex<float> * var) {
    return "single";
}
template<> std::string get_current_precision (const std::complex<double> * var) {
    return "double";
}

#if ((defined __CUDA) || (defined __ROCM))
template<> AbacusDevice_t get_device_type <DEVICE_GPU> (const DEVICE_GPU* dev) {
    return GpuDevice;
}

void set_device(const int rank) {
#if defined (__CUDA)
    cudaSetDevice(rank);
#elif defined (__ROCM)
    hipSetDevice(rank);
#endif
}

int get_device_num() {
    int device_num = -1;
#if defined (__CUDA)
    cudaGetDeviceCount(&device_num);
#elif defined (__ROCM)
    hipGetDeviceCount(&device_num);
#endif
    return device_num;
}
#endif

#if defined(__CUDA)
template<> void print_device_info <DEVICE_GPU> (const DEVICE_GPU* ctx, std::ofstream& ofs_device) {
  if (is_init) {return;}
  int deviceCount = 0;
  cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
  if (error_id != cudaSuccess) {
    ofs_device << "cudaGetDeviceCount returned "
               << static_cast<int>(error_id) << "\n-> "
               << cudaGetErrorString(error_id) << std::endl;
    ModuleBase::WARNING_QUIT("device", "GPU returned is without cudaSuccess");
  }
  // This function call returns 0 if there are no CUDA capable devices.
  if (deviceCount == 0) {
    ofs_device << "There are no available device(s) that support CUDA\n";
  } else {
    ofs_device << "Detected " << deviceCount << " CUDA Capable device(s)\n";
  }
  int dev = 0, driverVersion = 0, runtimeVersion = 0;
  cudaSetDevice(dev);
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, dev);
  ofs_device << "\nDevice " << dev << ":\t " << deviceProp.name << std::endl;
  // Console log
  cudaDriverGetVersion(&driverVersion);
  cudaRuntimeGetVersion(&runtimeVersion);
  char msg[1024];
  sprintf(msg,
          "  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
          driverVersion / 1000, (driverVersion % 100) / 10,
          runtimeVersion / 1000, (runtimeVersion % 100) / 10);
  ofs_device << msg << std::endl;
  sprintf(msg,
          "  CUDA Capability Major/Minor version number:    %d.%d\n",
          deviceProp.major, deviceProp.minor);
  ofs_device << msg << std::endl;
  sprintf(msg,
          "  GPU Max Clock rate:                            %.0f MHz (%0.2f "
          "GHz)\n",
          deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);
  ofs_device << msg << std::endl;
  // This is supported in CUDA 5.0 (runtime API device properties)
  sprintf(msg, 
          "  Memory Clock rate:                             %.0f Mhz\n",
          deviceProp.memoryClockRate * 1e-3f);
  ofs_device << msg << std::endl;   
   
  sprintf(msg, 
          "  Memory Bus Width:                              %d-bit\n",
          deviceProp.memoryBusWidth);
  ofs_device << msg << std::endl;   
  sprintf(msg,
      "  Maximum Texture Dimension Size (x,y,z)         1D=(%d), 2D=(%d, "
      "%d), 3D=(%d, %d, %d)\n",
      deviceProp.maxTexture1D, deviceProp.maxTexture2D[0],
      deviceProp.maxTexture2D[1], deviceProp.maxTexture3D[0],
      deviceProp.maxTexture3D[1], deviceProp.maxTexture3D[2]);
  ofs_device << msg << std::endl;   
  
  sprintf(msg, 
      "  Maximum Layered 1D Texture Size, (num) layers  1D=(%d), %d layers\n",
      deviceProp.maxTexture1DLayered[0], deviceProp.maxTexture1DLayered[1]);
  ofs_device << msg << std::endl;   
  sprintf(msg, 
      "  Maximum Layered 2D Texture Size, (num) layers  2D=(%d, %d), %d "
      "layers\n",
      deviceProp.maxTexture2DLayered[0], deviceProp.maxTexture2DLayered[1],
      deviceProp.maxTexture2DLayered[2]);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Total amount of constant memory:               %zu bytes\n",
         deviceProp.totalConstMem);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Total amount of shared memory per block:       %zu bytes\n",
         deviceProp.sharedMemPerBlock);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Total shared memory per multiprocessor:        %zu bytes\n",
         deviceProp.sharedMemPerMultiprocessor);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Total number of registers available per block: %d\n",
         deviceProp.regsPerBlock);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Warp size:                                     %d\n",
         deviceProp.warpSize);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Maximum number of threads per multiprocessor:  %d\n",
         deviceProp.maxThreadsPerMultiProcessor);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Maximum number of threads per block:           %d\n",
         deviceProp.maxThreadsPerBlock);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Max dimension size of a thread block (x,y,z): (%d, %d, %d)\n",
         deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1],
         deviceProp.maxThreadsDim[2]);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Max dimension size of a grid size    (x,y,z): (%d, %d, %d)\n",
         deviceProp.maxGridSize[0], deviceProp.maxGridSize[1],
         deviceProp.maxGridSize[2]);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Maximum memory pitch:                          %zu bytes\n",
         deviceProp.memPitch);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Texture alignment:                             %zu bytes\n",
         deviceProp.textureAlignment);
  ofs_device << msg << std::endl;   
  sprintf(msg, 
      "  Concurrent copy and kernel execution:          %s with %d copy "
      "engine(s)\n",
      (deviceProp.deviceOverlap ? "Yes" : "No"), deviceProp.asyncEngineCount);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Run time limit on kernels:                     %s\n",
         deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Integrated GPU sharing Host Memory:            %s\n",
         deviceProp.integrated ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Support host page-locked memory mapping:       %s\n",
         deviceProp.canMapHostMemory ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Alignment requirement for Surfaces:            %s\n",
         deviceProp.surfaceAlignment ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Device has ECC support:                        %s\n",
         deviceProp.ECCEnabled ? "Enabled" : "Disabled");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Device supports Unified Addressing (UVA):      %s\n",
         deviceProp.unifiedAddressing ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Device supports Managed Memory:                %s\n",
         deviceProp.managedMemory ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Device supports Compute Preemption:            %s\n",
         deviceProp.computePreemptionSupported ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Supports Cooperative Kernel Launch:            %s\n",
         deviceProp.cooperativeLaunch ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Supports MultiDevice Co-op Kernel Launch:      %s\n",
         deviceProp.cooperativeMultiDeviceLaunch ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Device PCI Domain ID / Bus ID / location ID:   %d / %d / %d\n",
         deviceProp.pciDomainID, deviceProp.pciBusID, deviceProp.pciDeviceID);
  ofs_device << msg << std::endl;   
  const char *sComputeMode[] = {
      "Default (multiple host threads can use ::cudaSetDevice() with device "
      "simultaneously)",
      "Exclusive (only one host thread in one process is able to use "
      "::cudaSetDevice() with this device)",
      "Prohibited (no host thread can use ::cudaSetDevice() with this "
      "device)",
      "Exclusive Process (many threads in one process is able to use "
      "::cudaSetDevice() with this device)",
      "Unknown", NULL};
  sprintf(msg, "  Compute Mode:\n");
  ofs_device << msg << std::endl;   
  ofs_device << "  " << sComputeMode[deviceProp.computeMode] << std::endl << std::endl;

  // If there are 2 or more GPUs, query to determine whether RDMA is supported
  if (deviceCount >= 2) {
    cudaDeviceProp prop[64];
    int gpuid[64];  // we want to find the first two GPUs that can support P2P
    int gpu_p2p_count = 0;

    for (int i = 0; i < deviceCount; i++) {
      cudaGetDeviceProperties(&prop[i], i);

      // Only boards based on Fermi or later can support P2P
      if (prop[i].major >= 2) {
        // This is an array of P2P capable GPUs
        gpuid[gpu_p2p_count++] = i;
      }
    }

    // Show all the combinations of support P2P GPUs
    int can_access_peer;

    if (gpu_p2p_count >= 2) {
      for (int i = 0; i < gpu_p2p_count; i++) {
        for (int j = 0; j < gpu_p2p_count; j++) {
          if (gpuid[i] == gpuid[j]) {
            continue;
          }
          cudaDeviceCanAccessPeer(&can_access_peer, gpuid[i], gpuid[j]);
          sprintf(msg, "> Peer access from %s (GPU%d) -> %s (GPU%d) : %s\n",
                 prop[gpuid[i]].name, gpuid[i], prop[gpuid[j]].name, gpuid[j],
                 can_access_peer ? "Yes" : "No");
          ofs_device << msg << std::endl;
        }
      }
    }
  }

  // csv masterlog info
  // *****************************
  // exe and CUDA driver name
  std::string sProfileString = "deviceQuery, CUDA Driver = CUDART";
  char cTemp[16];

  // driver version
  sProfileString += ", CUDA Driver Version = ";

  snprintf(cTemp, sizeof(cTemp), "%d.%d", driverVersion / 1000,
           (driverVersion % 100) / 10);
  sProfileString += cTemp;

  // Runtime version
  sProfileString += ", CUDA Runtime Version = ";
  snprintf(cTemp, sizeof(cTemp), "%d.%d", runtimeVersion / 1000,
           (runtimeVersion % 100) / 10);
  sProfileString += cTemp;

  // Device count
  sProfileString += ", NumDevs = ";
  snprintf(cTemp, sizeof(cTemp), "%d", deviceCount);
  sProfileString += cTemp;
  sProfileString += "\n";

  ofs_device << sProfileString.c_str() << std::endl;
  is_init = true;
  ofs_device << "End of device informations." << std::endl << std::endl;
}

template<> void record_device_memory<DEVICE_GPU> (const DEVICE_GPU* ctx, std::ofstream& ofs_device, std::string str, size_t size) {
  ofs_device << "Allocate " << static_cast<double>(size) / 8 / 1024 / 1024 << " \tMB device memory\t" 
             << "from " << str
             << std::endl << std::endl;
}

#elif defined(__ROCM)
template<> void print_device_info <DEVICE_GPU> (const DEVICE_GPU* ctx, std::ofstream& ofs_device) {
  if (is_init) {return;}
  int deviceCount = 0;
  hipError_t error_id = hipGetDeviceCount(&deviceCount);
  if (error_id != hipSuccess) {
    ofs_device << "hipGetDeviceCount returned "
               << static_cast<int>(error_id) << "\n-> "
               << hipGetErrorString(error_id) << std::endl;
    ModuleBase::WARNING_QUIT("device", "GPU returned is without hipSuccess");
  }
  // This function call returns 0 if there are no CUDA capable devices.
  if (deviceCount == 0) {
    ofs_device << "There are no available device(s) that support CUDA\n";
  } else {
    ofs_device << "Detected " << deviceCount << " CUDA Capable device(s)\n";
  }
  int dev = 0, driverVersion = 0, runtimeVersion = 0;
  hipSetDevice(dev);
  hipDeviceProp_t deviceProp;
  hipGetDeviceProperties(&deviceProp, dev);
  ofs_device << "\nDevice " << dev << ":\t " << deviceProp.name << std::endl;
  // Console log
  hipDriverGetVersion(&driverVersion);
  hipRuntimeGetVersion(&runtimeVersion);
  char msg[1024];
  sprintf(msg,
          "  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
          driverVersion / 1000, (driverVersion % 100) / 10,
          runtimeVersion / 1000, (runtimeVersion % 100) / 10);
  ofs_device << msg << std::endl;
  sprintf(msg,
          "  CUDA Capability Major/Minor version number:    %d.%d\n",
          deviceProp.major, deviceProp.minor);
  ofs_device << msg << std::endl;
  sprintf(msg,
          "  GPU Max Clock rate:                            %.0f MHz (%0.2f "
          "GHz)\n",
          deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);
  ofs_device << msg << std::endl;
  // This is supported in CUDA 5.0 (runtime API device properties)
  sprintf(msg, 
          "  Memory Clock rate:                             %.0f Mhz\n",
          deviceProp.memoryClockRate * 1e-3f);
  ofs_device << msg << std::endl;   
   
  sprintf(msg, 
          "  Memory Bus Width:                              %d-bit\n",
          deviceProp.memoryBusWidth);
  ofs_device << msg << std::endl;   
  sprintf(msg,
      "  Maximum Texture Dimension Size (x,y,z)         1D=(%d), 2D=(%d, "
      "%d), 3D=(%d, %d, %d)\n",
      deviceProp.maxTexture1D, deviceProp.maxTexture2D[0],
      deviceProp.maxTexture2D[1], deviceProp.maxTexture3D[0],
      deviceProp.maxTexture3D[1], deviceProp.maxTexture3D[2]);
  ofs_device << msg << std::endl;   
  
  sprintf(msg, "  Total amount of constant memory:               %zu bytes\n",
         deviceProp.totalConstMem);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Total amount of shared memory per block:       %zu bytes\n",
         deviceProp.sharedMemPerBlock);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Total number of registers available per block: %d\n",
         deviceProp.regsPerBlock);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Warp size:                                     %d\n",
         deviceProp.warpSize);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Maximum number of threads per multiprocessor:  %d\n",
         deviceProp.maxThreadsPerMultiProcessor);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Maximum number of threads per block:           %d\n",
         deviceProp.maxThreadsPerBlock);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Max dimension size of a thread block (x,y,z): (%d, %d, %d)\n",
         deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1],
         deviceProp.maxThreadsDim[2]);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Max dimension size of a grid size    (x,y,z): (%d, %d, %d)\n",
         deviceProp.maxGridSize[0], deviceProp.maxGridSize[1],
         deviceProp.maxGridSize[2]);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Maximum memory pitch:                          %zu bytes\n",
         deviceProp.memPitch);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Texture alignment:                             %zu bytes\n",
         deviceProp.textureAlignment);
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Run time limit on kernels:                     %s\n",
         deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Integrated GPU sharing Host Memory:            %s\n",
         deviceProp.integrated ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Support host page-locked memory mapping:       %s\n",
         deviceProp.canMapHostMemory ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Device has ECC support:                        %s\n",
         deviceProp.ECCEnabled ? "Enabled" : "Disabled");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Device supports Managed Memory:                %s\n",
         deviceProp.managedMemory ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Supports Cooperative Kernel Launch:            %s\n",
         deviceProp.cooperativeLaunch ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Supports MultiDevice Co-op Kernel Launch:      %s\n",
         deviceProp.cooperativeMultiDeviceLaunch ? "Yes" : "No");
  ofs_device << msg << std::endl;   
  sprintf(msg, "  Device PCI Domain ID / Bus ID / location ID:   %d / %d / %d\n",
         deviceProp.pciDomainID, deviceProp.pciBusID, deviceProp.pciDeviceID);
  ofs_device << msg << std::endl;   
  const char *sComputeMode[] = {
      "Default (multiple host threads can use ::hipSetDevice() with device "
      "simultaneously)",
      "Exclusive (only one host thread in one process is able to use "
      "::hipSetDevice() with this device)",
      "Prohibited (no host thread can use ::hipSetDevice() with this "
      "device)",
      "Exclusive Process (many threads in one process is able to use "
      "::hipSetDevice() with this device)",
      "Unknown", NULL};
  sprintf(msg, "  Compute Mode:\n");
  ofs_device << msg << std::endl;   
  ofs_device << "  " << sComputeMode[deviceProp.computeMode] << std::endl << std::endl;

  // If there are 2 or more GPUs, query to determine whether RDMA is supported
  if (deviceCount >= 2) {
    hipDeviceProp_t prop[64];
    int gpuid[64];  // we want to find the first two GPUs that can support P2P
    int gpu_p2p_count = 0;

    for (int i = 0; i < deviceCount; i++) {
      hipGetDeviceProperties(&prop[i], i);

      // Only boards based on Fermi or later can support P2P
      if (prop[i].major >= 2) {
        // This is an array of P2P capable GPUs
        gpuid[gpu_p2p_count++] = i;
      }
    }

    // Show all the combinations of support P2P GPUs
    int can_access_peer;

    if (gpu_p2p_count >= 2) {
      for (int i = 0; i < gpu_p2p_count; i++) {
        for (int j = 0; j < gpu_p2p_count; j++) {
          if (gpuid[i] == gpuid[j]) {
            continue;
          }
          hipDeviceCanAccessPeer(&can_access_peer, gpuid[i], gpuid[j]);
          sprintf(msg, "> Peer access from %s (GPU%d) -> %s (GPU%d) : %s\n",
                 prop[gpuid[i]].name, gpuid[i], prop[gpuid[j]].name, gpuid[j],
                 can_access_peer ? "Yes" : "No");
          ofs_device << msg << std::endl;
        }
      }
    }
  }

  // csv masterlog info
  // *****************************
  // exe and CUDA driver name
  std::string sProfileString = "deviceQuery, CUDA Driver = CUDART";
  char cTemp[16];

  // driver version
  sProfileString += ", CUDA Driver Version = ";

  snprintf(cTemp, sizeof(cTemp), "%d.%d", driverVersion / 1000,
           (driverVersion % 100) / 10);
  sProfileString += cTemp;

  // Runtime version
  sProfileString += ", CUDA Runtime Version = ";
  snprintf(cTemp, sizeof(cTemp), "%d.%d", runtimeVersion / 1000,
           (runtimeVersion % 100) / 10);
  sProfileString += cTemp;

  // Device count
  sProfileString += ", NumDevs = ";
  snprintf(cTemp, sizeof(cTemp), "%d", deviceCount);
  sProfileString += cTemp;
  sProfileString += "\n";

  ofs_device << sProfileString.c_str() << std::endl;
  is_init = true;
  ofs_device << "End of device informations." << std::endl << std::endl;
}

template<> void record_device_memory<DEVICE_GPU> (const DEVICE_GPU* ctx, std::ofstream& ofs_device, std::string str, size_t size) {
  ofs_device << "Allocate " << static_cast<double>(size) / 8 / 1024 / 1024 << " \tMB device memory\t" 
             << "from " << str
             << std::endl << std::endl;
}

#endif

#if __MPI
int stringCmp(const void *a, const void* b)
{
    char* m = (char*)a;
    char* n = (char*)b;
    int i, sum = 0;

    for(i = 0; i < MPI_MAX_PROCESSOR_NAME; i++)
        if (m[i] == n[i])
            continue;
        else
        {
            sum = m[i] - n[i];
            break;
        }
    return sum;
}

int get_node_rank() {
    char host_name[MPI_MAX_PROCESSOR_NAME];
    memset(host_name, '\0', sizeof(char) * MPI_MAX_PROCESSOR_NAME);
    char (*host_names)[MPI_MAX_PROCESSOR_NAME];
    int n, namelen, color, rank, nprocs, myrank;
    size_t bytes;
    MPI_Comm nodeComm;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Get_processor_name(host_name,&namelen);

    bytes = nprocs * sizeof(char[MPI_MAX_PROCESSOR_NAME]);
    host_names = (char (*)[MPI_MAX_PROCESSOR_NAME]) malloc(bytes);
    for (int ii = 0; ii < nprocs; ii++) {
        memset(host_names[ii], '\0', sizeof(char) * MPI_MAX_PROCESSOR_NAME);
    }

    strcpy(host_names[rank], host_name);

    for (n=0; n<nprocs; n++)
        MPI_Bcast(&(host_names[n]),MPI_MAX_PROCESSOR_NAME, MPI_CHAR, n, MPI_COMM_WORLD);
    qsort(host_names, nprocs,  sizeof(char[MPI_MAX_PROCESSOR_NAME]), stringCmp);

    color = 0;
    for (n=0; n<nprocs-1; n++)
    {
        if(strcmp(host_name, host_names[n]) == 0)
        {
            break;
        }
        if(strcmp(host_names[n], host_names[n+1]))
        {
            color++;
        }
    }

    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &nodeComm);
    MPI_Comm_rank(nodeComm, &myrank);

    MPI_Barrier(MPI_COMM_WORLD);
    int looprank=myrank;
    // printf (" Assigning device %d  to process on node %s rank %d, OK\n",looprank,  host_name, rank );
    free(host_names);
    return looprank;
}
#endif

std::string get_device_flag(const std::string& device, const std::string& ks_solver, const std::string& basis_type) {
    std::string str = "gpu";
    std::string env = "host";
    int device_num = -1;
#if ((defined __CUDA) || (defined __ROCM))
    device_num = device::get_device_num();
    if (device_num <= 0) {
        str = "cpu";
    }
    env = "device";
#else
    str = "cpu";
#endif
    if (ks_solver != "cg" && ks_solver != "dav" && ks_solver != "bpcg") {
        str = "cpu";
    }
    if (basis_type != "pw") {
        str = "cpu";
    }
    if (device == "cpu") {
        str = "cpu";
    }
    if (str == device) {
        return str;
    }
    else {
        std::string msg = "INPUT device setting does not match the request!";
        msg += "\n Input device = " + device;
        msg += "\n Input basis_type = " + basis_type;
        msg += "\n Input ks_solver = " + ks_solver;
        msg += "\n Compile setting = " + env;
        msg += "\n Environment device_num = " + std::to_string(device_num) + "\n";
        ModuleBase::WARNING_QUIT("device", msg);
        return "unknown";
    }
}

int get_device_kpar(const int& kpar) {
#if __MPI && (__CUDA || __ROCM)
    int temp_nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &temp_nproc);
    if (temp_nproc != kpar)
    {
        ModuleBase::WARNING("Input_conv", "None kpar set in INPUT file, auto set kpar value.");
    }
    // GlobalV::KPAR = temp_nproc;
    // band the CPU processor to the devices
    int node_rank = psi::device::get_node_rank();
    int device_num = psi::device::get_device_num();
    psi::device::set_device(node_rank % device_num);
    return temp_nproc;
#endif
    return kpar;
}

} // end of namespace device
} // end of namespace psi