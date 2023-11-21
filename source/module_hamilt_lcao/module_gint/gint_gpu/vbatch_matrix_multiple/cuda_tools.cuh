#ifndef CUDA_TOOLS_CUH
#define CUDA_TOOLS_CUH

#include <cuda_runtime.h>
#include <fstream>
#include <sstream>
#include <iostream>

cudaError_t checkCuda(cudaError_t result);

void dump_cuda_array_to_file(double *cuda_array, int width, int hight, const std::string &filename);


#endif // CUDA_TOOLS_CUH#ifndef CUDA_TOOLS_CUH
