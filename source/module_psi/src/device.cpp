
#include <stdio.h>
#include <complex>
#include <fstream>
#include <iostream>
#include "module_psi/psi.h"
#include "module_psi/include/types.h"
#include "module_psi/include/device.h"
#include "module_base/tool_quit.h"

#if defined(__CUDA)
#include <cuda_runtime.h>
#endif

namespace psi{

namespace device{

static bool is_init = false;

// functions used in custom ops
template<> AbacusDevice_t get_device_type <DEVICE_CPU> (const DEVICE_CPU* dev) {
    return CpuDevice;
}

#if ((defined __CUDA) || (defined __ROCM))
template<> AbacusDevice_t get_device_type <DEVICE_GPU> (const DEVICE_GPU* dev) {
    return GpuDevice;
}
#endif

#if defined(__CUDA) || defined(__ROCM)
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
  char msg[256];
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

#endif

} // end of namespace device
} // end of namespace psi