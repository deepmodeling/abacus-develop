#include "cuda_compat.h"

namespace ModuleBase {
namespace cuda_compat {

//---------------------------------------------------------------------------
// Implementation of printDeprecatedDeviceInfo and printComputeModeInfo
//---------------------------------------------------------------------------
void printDeprecatedDeviceInfo(std::ostream& ofs_device, const cudaDeviceProp& deviceProp) 
{
#if defined(CUDA_VERSION) && CUDA_VERSION < 13000
  char msg[1024];
  sprintf(msg,
          "  GPU Max Clock rate:                            %.0f MHz (%0.2f "
          "GHz)\n",
          deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);
  ofs_device << msg << std::endl;
  // This is supported in CUDA 5.0 (runtime API device properties)
  sprintf(msg, "  Memory Clock rate:                             %.0f Mhz\n",
          deviceProp.memoryClockRate * 1e-3f);
  ofs_device << msg << std::endl;

  sprintf(msg, "  Memory Bus Width:                              %d-bit\n",
          deviceProp.memoryBusWidth);
  ofs_device << msg << std::endl;
  
  sprintf(msg,
          "  Concurrent copy and kernel execution:          %s with %d copy "
          "engine(s)\n",
          (deviceProp.deviceOverlap ? "Yes" : "No"),
          deviceProp.asyncEngineCount);
  ofs_device << msg << std::endl;
  sprintf(msg, "  Run time limit on kernels:                     %s\n",
          deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
  ofs_device << msg << std::endl;
#endif
}

void printComputeModeInfo(std::ostream& ofs_device, const cudaDeviceProp& deviceProp) 
{
#if defined(CUDA_VERSION) && CUDA_VERSION < 13000
  char msg[1024];
  sprintf(msg, "  Supports MultiDevice Co-op Kernel Launch:      %s\n",
          deviceProp.cooperativeMultiDeviceLaunch ? "Yes" : "No");
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
      "Unknown",
      NULL};
  sprintf(msg, "  Compute Mode:\n");
  ofs_device << msg << std::endl;
  ofs_device << "  " << sComputeMode[deviceProp.computeMode] << std::endl
             << std::endl;
#endif
}

} // namespace cuda_compat
} // namespace ModuleBase
