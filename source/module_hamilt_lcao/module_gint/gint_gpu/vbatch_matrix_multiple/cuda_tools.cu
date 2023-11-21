#include "module_hamilt_lcao/module_gint/gint_gpu/vbatch_matrix_multiple/cuda_tools.cuh"

cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(__DEBUG)
    if (result != cudaSuccess)
    {
        fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
#endif
    return result;
}

void dump_cuda_array_to_file(double *cuda_array, int width, int hight, const std::string &filename)
{
    double *h_data = new double[width * hight];
    cudaMemcpy(h_data, cuda_array, width * hight * sizeof(double), cudaMemcpyDeviceToHost);

    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
        std::cerr << "Failed to open file for writing." << std::endl;
    }
    for (int j = 0; j < hight; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            outFile << "hight" << j << "   width:" << i << "   " << h_data[j * width + i] << std::endl;
        }
    }
    outFile.close();
    delete[] h_data;
}