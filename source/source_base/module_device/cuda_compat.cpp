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

//-------------------------------------------------------------------------------------------------
// Implementation of cufftGetErrorStringCompat
//-------------------------------------------------------------------------------------------------
const char* cufftGetErrorStringCompat(cufftResult_t error) 
{
    switch (error)
    {
    case CUFFT_SUCCESS:
        return "CUFFT_SUCCESS";
    case CUFFT_INVALID_PLAN:
        return "CUFFT_INVALID_PLAN";
    case CUFFT_ALLOC_FAILED:
        return "CUFFT_ALLOC_FAILED";
    case CUFFT_INVALID_TYPE:
        return "CUFFT_INVALID_TYPE";
    case CUFFT_INVALID_VALUE:
        return "CUFFT_INVALID_VALUE";
    case CUFFT_INTERNAL_ERROR:
        return "CUFFT_INTERNAL_ERROR";
    case CUFFT_EXEC_FAILED:
        return "CUFFT_EXEC_FAILED";
    case CUFFT_SETUP_FAILED:
        return "CUFFT_SETUP_FAILED";
    case CUFFT_INVALID_SIZE:
        return "CUFFT_INVALID_SIZE";
    case CUFFT_UNALIGNED_DATA:
        return "CUFFT_UNALIGNED_DATA";
    case CUFFT_INVALID_DEVICE:
        return "CUFFT_INVALID_DEVICE";
    case CUFFT_NO_WORKSPACE:
        return "CUFFT_NO_WORKSPACE";
    case CUFFT_NOT_IMPLEMENTED:
        return "CUFFT_NOT_IMPLEMENTED";
    case CUFFT_NOT_SUPPORTED:
        return "CUFFT_NOT_SUPPORTED";

#if defined(CUDA_VERSION) && CUDA_VERSION < 13000
    case CUFFT_INCOMPLETE_PARAMETER_LIST:
        return "CUFFT_INCOMPLETE_PARAMETER_LIST";
    case CUFFT_PARSE_ERROR:
        return "CUFFT_PARSE_ERROR";
    case CUFFT_LICENSE_ERROR:
        return "CUFFT_LICENSE_ERROR";
#endif
    
    default:
        return "<unknown>";
    }
}

} // namespace cuda_compat
} // namespace ModuleBase
