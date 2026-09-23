#ifndef ESOLVER_NEP_GPU_KERNELS_H
#define ESOLVER_NEP_GPU_KERNELS_H

class NEP;
struct NEP_GPU_Model;

struct NEP_GPU_Result
{
    double* potential;
    double* force;
    double* virial;
};

int nep_gpu_select_device(int local_rank, int* device_id, char* message, int message_size);
NEP_GPU_Model* nep_gpu_create(const NEP* nep, char* message, int message_size);
void nep_gpu_destroy(NEP_GPU_Model* model);
int nep_gpu_compute(NEP_GPU_Model* model,
                    int nlocal,
                    int nall,
                    int max_radial,
                    int max_angular,
                    const int* type,
                    const double* x,
                    const double* y,
                    const double* z,
                    const int* nn_radial,
                    const int* nl_radial,
                    const int* nn_angular,
                    const int* nl_angular,
                    NEP_GPU_Result result,
                    char* message,
                    int message_size);

int nep_gpu_compute_device_neighbors(NEP_GPU_Model* model,
                                     int nlocal,
                                     int nall,
                                     int max_radial,
                                     int max_angular,
                                     const int* type,
                                     const double* x,
                                     const double* y,
                                     const double* z,
                                     const void* device_nn_radial,
                                     const void* device_nl_radial,
                                     const void* device_nn_angular,
                                     const void* device_nl_angular,
                                     NEP_GPU_Result result,
                                     char* message,
                                     int message_size);

#endif
