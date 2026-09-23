/* ----------------------------------------------------------------------
   CUDA backend for ABACUS NEP molecular dynamics
------------------------------------------------------------------------- */

#include "nep_gpu_kernels.h"

#include "nep.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <vector>

#define NUM_ELEMENTS 94
#include "nep_gpu_math.cuh"

namespace {

struct GPUParaMB {
  int version = 4;
  int model_type = 0;
  float rc_radial[NUM_ELEMENTS];
  float rc_angular[NUM_ELEMENTS];
  int n_max_radial = 0;
  int n_max_angular = 0;
  int L_max = 0;
  int dim_angular = 0;
  int num_L = 0;
  int basis_size_radial = 8;
  int basis_size_angular = 8;
  int num_types_sq = 0;
  int num_c_radial = 0;
  int num_types = 0;
};

struct GPUANN {
  int dim = 0;
  int num_neurons1 = 0;
  int num_para = 0;
  int num_para_ann = 0;
  const float *w0[NUM_ELEMENTS];
  const float *b0[NUM_ELEMENTS];
  const float *w1[NUM_ELEMENTS];
  const float *b1;
  const float *c;
  const float *q_scaler;
};

}    // namespace

struct NEP_GPU_Model {
  GPUParaMB paramb;
  GPUANN annmb;
  float *d_parameters = nullptr;
  int *d_type = nullptr;
  int *d_nn_radial = nullptr;
  int *d_nl_radial = nullptr;
  int *d_nn_angular = nullptr;
  int *d_nl_angular = nullptr;
  double *d_x = nullptr;
  double *d_y = nullptr;
  double *d_z = nullptr;
  double *d_pe = nullptr;
  double *d_fx = nullptr;
  double *d_fy = nullptr;
  double *d_fz = nullptr;
  double *d_virial = nullptr;
  float *d_Fp = nullptr;
  float *d_sum_fxyz = nullptr;
  float *d_f12x = nullptr;
  float *d_f12y = nullptr;
  float *d_f12z = nullptr;
  int nall_capacity = 0;
  int max_radial_capacity = 0;
  int max_angular_capacity = 0;
};

namespace {

int set_message(char *message, int message_size, const char *text)
{
  if (message && message_size > 0) snprintf(message, message_size, "%s", text);
  return 1;
}

int check_cuda(cudaError_t status, char *message, int message_size, const char *where)
{
  if (status == cudaSuccess) return 0;
  if (message && message_size > 0)
    snprintf(message, message_size, "%s: %s", where, cudaGetErrorString(status));
  return 1;
}

void update_ann_pointers(NEP_GPU_Model *model)
{
  float *pointer = model->d_parameters;
  for (int t = 0; t < model->paramb.num_types; ++t) {
    if (t > 0 && model->paramb.version == 3) {
      // NEP3 uses one shared neural-network parameter block for all types.
      pointer -= (model->annmb.dim + 2) * model->annmb.num_neurons1;
    }
    model->annmb.w0[t] = pointer;
    pointer += model->annmb.num_neurons1 * model->annmb.dim;
    model->annmb.b0[t] = pointer;
    pointer += model->annmb.num_neurons1;
    model->annmb.w1[t] = pointer;
    pointer += model->annmb.num_neurons1;
    if (model->paramb.version == 5) pointer += 1;
  }
  model->annmb.b1 = pointer;
  pointer += 1;
  model->annmb.c = pointer;
  model->annmb.q_scaler = model->d_parameters + model->annmb.num_para;
}

template <typename T> void free_device(T *&ptr)
{
  if (ptr) cudaFree(ptr);
  ptr = nullptr;
}

template <typename T>
int ensure_capacity(T *&ptr, int current, int requested, char *message, int message_size,
                    const char *name)
{
  if (requested <= current) return 0;
  free_device(ptr);
  return check_cuda(cudaMalloc(reinterpret_cast<void **>(&ptr), sizeof(T) * requested), message,
                    message_size, name);
}

int ensure_workspace(NEP_GPU_Model *model, int nall, int max_radial, int max_angular,
                     char *message, int message_size)
{
  if (ensure_capacity(model->d_type, model->nall_capacity, nall, message, message_size, "type"))
    return 1;
  if (ensure_capacity(model->d_x, model->nall_capacity, nall, message, message_size, "x"))
    return 1;
  if (ensure_capacity(model->d_y, model->nall_capacity, nall, message, message_size, "y"))
    return 1;
  if (ensure_capacity(model->d_z, model->nall_capacity, nall, message, message_size, "z"))
    return 1;
  if (ensure_capacity(model->d_pe, model->nall_capacity, nall, message, message_size, "pe"))
    return 1;
  if (ensure_capacity(model->d_fx, model->nall_capacity, nall, message, message_size, "fx"))
    return 1;
  if (ensure_capacity(model->d_fy, model->nall_capacity, nall, message, message_size, "fy"))
    return 1;
  if (ensure_capacity(model->d_fz, model->nall_capacity, nall, message, message_size, "fz"))
    return 1;
  if (ensure_capacity(model->d_virial, model->nall_capacity * 9, nall * 9, message, message_size,
                      "virial"))
    return 1;
  if (ensure_capacity(model->d_Fp, model->nall_capacity * model->annmb.dim,
                      nall * model->annmb.dim, message, message_size, "Fp"))
    return 1;

  const int num_abc = (model->paramb.L_max + 1) * (model->paramb.L_max + 1) - 1;
  const int sum_size = nall * (model->paramb.n_max_angular + 1) * num_abc;
  const int old_sum_size =
      model->nall_capacity * (model->paramb.n_max_angular + 1) * num_abc;
  if (ensure_capacity(model->d_sum_fxyz, old_sum_size, sum_size, message, message_size,
                      "sum_fxyz"))
    return 1;

  if (ensure_capacity(model->d_nn_radial, model->nall_capacity, nall, message, message_size,
                      "nn_radial"))
    return 1;
  if (ensure_capacity(model->d_nn_angular, model->nall_capacity, nall, message, message_size,
                      "nn_angular"))
    return 1;
  if (ensure_capacity(model->d_nl_radial, model->nall_capacity * model->max_radial_capacity,
                      nall * max_radial, message, message_size, "nl_radial"))
    return 1;
  if (ensure_capacity(model->d_nl_angular, model->nall_capacity * model->max_angular_capacity,
                      nall * max_angular, message, message_size, "nl_angular"))
    return 1;
  if (ensure_capacity(model->d_f12x, model->nall_capacity * model->max_angular_capacity,
                      nall * max_angular, message, message_size, "f12x"))
    return 1;
  if (ensure_capacity(model->d_f12y, model->nall_capacity * model->max_angular_capacity,
                      nall * max_angular, message, message_size, "f12y"))
    return 1;
  if (ensure_capacity(model->d_f12z, model->nall_capacity * model->max_angular_capacity,
                      nall * max_angular, message, message_size, "f12z"))
    return 1;

  model->nall_capacity = std::max(model->nall_capacity, nall);
  model->max_radial_capacity = std::max(model->max_radial_capacity, max_radial);
  model->max_angular_capacity = std::max(model->max_angular_capacity, max_angular);
  return 0;
}

__global__ void find_descriptor_lammps(GPUParaMB paramb, GPUANN annmb, int N, const int *g_NN,
                                       const int *g_NL, const int *g_NN_angular,
                                       const int *g_NL_angular, const int *__restrict__ g_type,
                                       const double *__restrict__ g_x,
                                       const double *__restrict__ g_y,
                                       const double *__restrict__ g_z, double *g_pe, float *g_Fp,
                                       float *g_sum_fxyz)
{
  int n1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (n1 >= N) return;

  int t1 = g_type[n1];
  double x1 = g_x[n1];
  double y1 = g_y[n1];
  double z1 = g_z[n1];
  float q[MAX_DIM] = {0.0f};

  for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
    int n2 = g_NL[n1 + N * i1];
    float x12 = g_x[n2] - x1;
    float y12 = g_y[n2] - y1;
    float z12 = g_z[n2] - z1;
    float d12 = sqrtf(x12 * x12 + y12 * y12 + z12 * z12);
    float fc12;
    int t2 = g_type[n2];
    float rc = (paramb.rc_radial[t1] + paramb.rc_radial[t2]) * 0.5f;
    float rcinv = 1.0f / rc;
    find_fc(rc, rcinv, d12, fc12);
    float fn12[MAX_NUM_N];
    find_fn(paramb.basis_size_radial, rcinv, d12, fc12, fn12);
    for (int n = 0; n <= paramb.n_max_radial; ++n) {
      float gn12 = 0.0f;
      for (int k = 0; k <= paramb.basis_size_radial; ++k) {
        int c_index = (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq;
        c_index += t1 * paramb.num_types + t2;
        gn12 += fn12[k] * annmb.c[c_index];
      }
      q[n] += gn12;
    }
  }

  const int num_abc = (paramb.L_max + 1) * (paramb.L_max + 1) - 1;
  for (int n = 0; n <= paramb.n_max_angular; ++n) {
    float s[NUM_OF_ABC] = {0.0f};
    for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
      int n2 = g_NL_angular[n1 + N * i1];
      float x12 = g_x[n2] - x1;
      float y12 = g_y[n2] - y1;
      float z12 = g_z[n2] - z1;
      float d12 = sqrtf(x12 * x12 + y12 * y12 + z12 * z12);
      float fc12;
      int t2 = g_type[n2];
      float rc = (paramb.rc_angular[t1] + paramb.rc_angular[t2]) * 0.5f;
      float rcinv = 1.0f / rc;
      find_fc(rc, rcinv, d12, fc12);
      float fn12[MAX_NUM_N];
      find_fn(paramb.basis_size_angular, rcinv, d12, fc12, fn12);
      float gn12 = 0.0f;
      for (int k = 0; k <= paramb.basis_size_angular; ++k) {
        int c_index = (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq;
        c_index += t1 * paramb.num_types + t2 + paramb.num_c_radial;
        gn12 += fn12[k] * annmb.c[c_index];
      }
      accumulate_s(paramb.L_max, d12, x12, y12, z12, gn12, s);
    }
    find_q(paramb.L_max, paramb.num_L, paramb.n_max_angular + 1, n, s,
           q + (paramb.n_max_radial + 1));
    for (int abc = 0; abc < num_abc; ++abc)
      g_sum_fxyz[(n * num_abc + abc) * N + n1] = s[abc];
  }

  for (int d = 0; d < annmb.dim; ++d) q[d] *= annmb.q_scaler[d];

  float F = 0.0f;
  float Fp[MAX_DIM] = {0.0f};
  if (paramb.version == 5) {
    apply_ann_one_layer_nep5(annmb.dim, annmb.num_neurons1, annmb.w0[t1], annmb.b0[t1],
                             annmb.w1[t1], annmb.b1, q, F, Fp);
  } else {
    apply_ann_one_layer(annmb.dim, annmb.num_neurons1, annmb.w0[t1], annmb.b0[t1],
                        annmb.w1[t1], annmb.b1, q, F, Fp);
  }
  g_pe[n1] += F;
  for (int d = 0; d < annmb.dim; ++d) g_Fp[d * N + n1] = Fp[d] * annmb.q_scaler[d];
}

__global__ void find_force_radial_lammps(GPUParaMB paramb, GPUANN annmb, int N, int nlocal,
                                         const int *g_NN, const int *g_NL,
                                         const int *__restrict__ g_type,
                                         const double *__restrict__ g_x,
                                         const double *__restrict__ g_y,
                                         const double *__restrict__ g_z,
                                         const float *__restrict__ g_Fp, double *g_fx,
                                         double *g_fy, double *g_fz, double *g_virial)
{
  int n1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (n1 >= nlocal) return;

  int t1 = g_type[n1];
  float s_fx = 0.0f, s_fy = 0.0f, s_fz = 0.0f;
  float s_sxx = 0.0f, s_sxy = 0.0f, s_sxz = 0.0f;
  float s_syx = 0.0f, s_syy = 0.0f, s_syz = 0.0f;
  float s_szx = 0.0f, s_szy = 0.0f, s_szz = 0.0f;
  double x1 = g_x[n1], y1 = g_y[n1], z1 = g_z[n1];
  for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
    int n2 = g_NL[n1 + N * i1];
    int t2 = g_type[n2];
    float x12 = g_x[n2] - x1;
    float y12 = g_y[n2] - y1;
    float z12 = g_z[n2] - z1;
    float r12[3] = {x12, y12, z12};
    float d12 = sqrtf(x12 * x12 + y12 * y12 + z12 * z12);
    float d12inv = 1.0f / d12;
    float f12[3] = {0.0f};
    float f21[3] = {0.0f};
    float fc12, fcp12;
    float rc = (paramb.rc_radial[t1] + paramb.rc_radial[t2]) * 0.5f;
    float rcinv = 1.0f / rc;
    find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);
    float fn12[MAX_NUM_N], fnp12[MAX_NUM_N];
    find_fn_and_fnp(paramb.basis_size_radial, rcinv, d12, fc12, fcp12, fn12, fnp12);
    for (int n = 0; n <= paramb.n_max_radial; ++n) {
      float gnp12 = 0.0f, gnp21 = 0.0f;
      for (int k = 0; k <= paramb.basis_size_radial; ++k) {
        int c_index = (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq;
        gnp12 += fnp12[k] * annmb.c[c_index + t1 * paramb.num_types + t2];
        gnp21 += fnp12[k] * annmb.c[c_index + t2 * paramb.num_types + t1];
      }
      float tmp12 = g_Fp[n1 + n * N] * gnp12 * d12inv;
      float tmp21 = g_Fp[n2 + n * N] * gnp21 * d12inv;
      for (int d = 0; d < 3; ++d) {
        f12[d] += tmp12 * r12[d];
        f21[d] -= tmp21 * r12[d];
      }
    }
    s_fx += f12[0] - f21[0];
    s_fy += f12[1] - f21[1];
    s_fz += f12[2] - f21[2];
    s_sxx += r12[0] * f21[0];
    s_syy += r12[1] * f21[1];
    s_szz += r12[2] * f21[2];
    s_sxy += r12[0] * f21[1];
    s_sxz += r12[0] * f21[2];
    s_syx += r12[1] * f21[0];
    s_syz += r12[1] * f21[2];
    s_szx += r12[2] * f21[0];
    s_szy += r12[2] * f21[1];
  }
  g_fx[n1] += s_fx;
  g_fy[n1] += s_fy;
  g_fz[n1] += s_fz;
  g_virial[n1 + 0 * N] += s_sxx;
  g_virial[n1 + 1 * N] += s_syy;
  g_virial[n1 + 2 * N] += s_szz;
  g_virial[n1 + 3 * N] += s_sxy;
  g_virial[n1 + 4 * N] += s_sxz;
  g_virial[n1 + 5 * N] += s_syz;
  g_virial[n1 + 6 * N] += s_syx;
  g_virial[n1 + 7 * N] += s_szx;
  g_virial[n1 + 8 * N] += s_szy;
}

__global__ void find_partial_force_angular_lammps(GPUParaMB paramb, GPUANN annmb, int N,
                                                  const int *g_NN_angular,
                                                  const int *g_NL_angular,
                                                  const int *__restrict__ g_type,
                                                  const double *__restrict__ g_x,
                                                  const double *__restrict__ g_y,
                                                  const double *__restrict__ g_z,
                                                  const float *__restrict__ g_Fp,
                                                  const float *__restrict__ g_sum_fxyz,
                                                  float *g_f12x, float *g_f12y, float *g_f12z)
{
  int n1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (n1 >= N) return;

  float Fp[MAX_DIM_ANGULAR] = {0.0f};
  float sum_fxyz[NUM_OF_ABC * MAX_NUM_N];
  const int num_abc = (paramb.L_max + 1) * (paramb.L_max + 1) - 1;
  for (int d = 0; d < paramb.dim_angular; ++d)
    Fp[d] = g_Fp[(paramb.n_max_radial + 1 + d) * N + n1];
  for (int n = 0; n < paramb.n_max_angular + 1; ++n) {
    for (int abc = 0; abc < num_abc; ++abc)
      sum_fxyz[n * NUM_OF_ABC + abc] = g_sum_fxyz[(n * num_abc + abc) * N + n1];
  }

  int t1 = g_type[n1];
  double x1 = g_x[n1], y1 = g_y[n1], z1 = g_z[n1];
  for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
    int index = i1 * N + n1;
    int n2 = g_NL_angular[n1 + N * i1];
    float x12 = g_x[n2] - x1;
    float y12 = g_y[n2] - y1;
    float z12 = g_z[n2] - z1;
    float r12[3] = {x12, y12, z12};
    float d12 = sqrtf(x12 * x12 + y12 * y12 + z12 * z12);
    float f12[3] = {0.0f};
    float fc12, fcp12;
    int t2 = g_type[n2];
    float rc = (paramb.rc_angular[t1] + paramb.rc_angular[t2]) * 0.5f;
    float rcinv = 1.0f / rc;
    find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);
    float fn12[MAX_NUM_N], fnp12[MAX_NUM_N];
    find_fn_and_fnp(paramb.basis_size_angular, rcinv, d12, fc12, fcp12, fn12, fnp12);
    for (int n = 0; n <= paramb.n_max_angular; ++n) {
      float gn12 = 0.0f, gnp12 = 0.0f;
      for (int k = 0; k <= paramb.basis_size_angular; ++k) {
        int c_index = (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq;
        c_index += t1 * paramb.num_types + t2 + paramb.num_c_radial;
        gn12 += fn12[k] * annmb.c[c_index];
        gnp12 += fnp12[k] * annmb.c[c_index];
      }
      accumulate_f12(paramb.L_max, paramb.num_L, n, paramb.n_max_angular + 1, d12, r12, gn12,
                     gnp12, Fp, sum_fxyz, f12);
    }
    g_f12x[index] = f12[0];
    g_f12y[index] = f12[1];
    g_f12z[index] = f12[2];
  }
}

__global__ void find_force_many_body_lammps(int N, int nlocal, const int *g_NN,
                                            const int *g_NL,
                                            const float *__restrict__ g_f12x,
                                            const float *__restrict__ g_f12y,
                                            const float *__restrict__ g_f12z,
                                            const double *__restrict__ g_x,
                                            const double *__restrict__ g_y,
                                            const double *__restrict__ g_z, double *g_fx,
                                            double *g_fy, double *g_fz, double *g_virial)
{
  int n1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (n1 >= nlocal) return;

  float s_fx = 0.0f, s_fy = 0.0f, s_fz = 0.0f;
  float s_sxx = 0.0f, s_sxy = 0.0f, s_sxz = 0.0f;
  float s_syx = 0.0f, s_syy = 0.0f, s_syz = 0.0f;
  float s_szx = 0.0f, s_szy = 0.0f, s_szz = 0.0f;
  double x1 = g_x[n1], y1 = g_y[n1], z1 = g_z[n1];

  for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
    int index = i1 * N + n1;
    int n2 = g_NL[index];
    float x12 = g_x[n2] - x1;
    float y12 = g_y[n2] - y1;
    float z12 = g_z[n2] - z1;
    float f12x = g_f12x[index];
    float f12y = g_f12y[index];
    float f12z = g_f12z[index];

    float f21x = 0.0f, f21y = 0.0f, f21z = 0.0f;
    for (int k = 0; k < g_NN[n2]; ++k) {
      if (g_NL[n2 + k * N] == n1) {
        int reverse = n2 + k * N;
        f21x = g_f12x[reverse];
        f21y = g_f12y[reverse];
        f21z = g_f12z[reverse];
        break;
      }
    }

    s_fx += f12x - f21x;
    s_fy += f12y - f21y;
    s_fz += f12z - f21z;
    s_sxx += x12 * f21x;
    s_syy += y12 * f21y;
    s_szz += z12 * f21z;
    s_sxy += x12 * f21y;
    s_sxz += x12 * f21z;
    s_syx += y12 * f21x;
    s_syz += y12 * f21z;
    s_szx += z12 * f21x;
    s_szy += z12 * f21y;
  }

  g_fx[n1] += s_fx;
  g_fy[n1] += s_fy;
  g_fz[n1] += s_fz;
  g_virial[n1 + 0 * N] += s_sxx;
  g_virial[n1 + 1 * N] += s_syy;
  g_virial[n1 + 2 * N] += s_szz;
  g_virial[n1 + 3 * N] += s_sxy;
  g_virial[n1 + 4 * N] += s_sxz;
  g_virial[n1 + 5 * N] += s_syz;
  g_virial[n1 + 6 * N] += s_syx;
  g_virial[n1 + 7 * N] += s_szx;
  g_virial[n1 + 8 * N] += s_szy;
}

}    // namespace

int nep_gpu_select_device(int local_rank, int *device_id, char *message, int message_size)
{
  int count = 0;
  if (check_cuda(cudaGetDeviceCount(&count), message, message_size, "cudaGetDeviceCount")) return 1;
  if (count <= 0) return set_message(message, message_size, "no visible CUDA devices");
  const int selected = local_rank % count;
  if (check_cuda(cudaSetDevice(selected), message, message_size, "cudaSetDevice")) return 1;
  if (device_id) *device_id = selected;
  return 0;
}

NEP_GPU_Model *nep_gpu_create(const NEP *nep, char *message, int message_size)
{
  NEP_GPU_Model *model = new NEP_GPU_Model();
  model->paramb.version = nep->paramb.version;
  model->paramb.model_type = nep->paramb.model_type;
  model->paramb.n_max_radial = nep->paramb.n_max_radial;
  model->paramb.n_max_angular = nep->paramb.n_max_angular;
  model->paramb.L_max = nep->paramb.L_max;
  model->paramb.dim_angular = nep->paramb.dim_angular;
  model->paramb.num_L = nep->paramb.num_L;
  model->paramb.basis_size_radial = nep->paramb.basis_size_radial;
  model->paramb.basis_size_angular = nep->paramb.basis_size_angular;
  model->paramb.num_types_sq = nep->paramb.num_types_sq;
  model->paramb.num_c_radial = nep->paramb.num_c_radial;
  model->paramb.num_types = static_cast<int>(nep->paramb.num_types);
  for (int t = 0; t < NUM_ELEMENTS; ++t) {
    model->paramb.rc_radial[t] = static_cast<float>(nep->paramb.rc_radial[t]);
    model->paramb.rc_angular[t] = static_cast<float>(nep->paramb.rc_angular[t]);
  }

  model->annmb.dim = nep->annmb.dim;
  model->annmb.num_neurons1 = nep->annmb.num_neurons1;
  model->annmb.num_para = nep->annmb.num_para;
  model->annmb.num_para_ann = nep->annmb.num_para_ann;

  std::vector<float> parameters(nep->annmb.num_para + nep->annmb.dim);
  for (int i = 0; i < nep->annmb.num_para; ++i)
    parameters[i] = static_cast<float>(nep->parameters[i]);
  for (int i = 0; i < nep->annmb.dim; ++i)
    parameters[nep->annmb.num_para + i] = static_cast<float>(nep->paramb.q_scaler[i]);

  if (check_cuda(cudaMalloc(reinterpret_cast<void **>(&model->d_parameters),
                            sizeof(float) * parameters.size()),
                 message, message_size, "model parameters")) {
    delete model;
    return nullptr;
  }
  if (check_cuda(cudaMemcpy(model->d_parameters, parameters.data(), sizeof(float) * parameters.size(),
                            cudaMemcpyHostToDevice),
                 message, message_size, "copy model parameters")) {
    nep_gpu_destroy(model);
    return nullptr;
  }
  update_ann_pointers(model);
  return model;
}

void nep_gpu_destroy(NEP_GPU_Model *model)
{
  if (!model) return;
  free_device(model->d_parameters);
  free_device(model->d_type);
  free_device(model->d_nn_radial);
  free_device(model->d_nl_radial);
  free_device(model->d_nn_angular);
  free_device(model->d_nl_angular);
  free_device(model->d_x);
  free_device(model->d_y);
  free_device(model->d_z);
  free_device(model->d_pe);
  free_device(model->d_fx);
  free_device(model->d_fy);
  free_device(model->d_fz);
  free_device(model->d_virial);
  free_device(model->d_Fp);
  free_device(model->d_sum_fxyz);
  free_device(model->d_f12x);
  free_device(model->d_f12y);
  free_device(model->d_f12z);
  delete model;
}

static int nep_gpu_compute_impl(NEP_GPU_Model *model, int nlocal, int nall, int max_radial,
                                int max_angular, const int *type, const double *x, const double *y,
                                const double *z, const int *nn_radial, const int *nl_radial,
                                const int *nn_angular, const int *nl_angular,
                                bool device_neighbors,
                                const int *device_nn_radial, const int *device_nl_radial,
                                const int *device_nn_angular, const int *device_nl_angular,
                                NEP_GPU_Result result, char *message, int message_size)
{
  if (!model) return set_message(message, message_size, "GPU model is null");
  if (ensure_workspace(model, nall, max_radial, max_angular, message, message_size)) return 1;

  const int *d_nn_radial = model->d_nn_radial;
  const int *d_nl_radial = model->d_nl_radial;
  const int *d_nn_angular = model->d_nn_angular;
  const int *d_nl_angular = model->d_nl_angular;

  if (check_cuda(cudaMemcpy(model->d_type, type, sizeof(int) * nall, cudaMemcpyHostToDevice),
                 message, message_size, "copy type"))
    return 1;
  if (check_cuda(cudaMemcpy(model->d_x, x, sizeof(double) * nall, cudaMemcpyHostToDevice), message,
                 message_size, "copy x"))
    return 1;
  if (check_cuda(cudaMemcpy(model->d_y, y, sizeof(double) * nall, cudaMemcpyHostToDevice), message,
                 message_size, "copy y"))
    return 1;
  if (check_cuda(cudaMemcpy(model->d_z, z, sizeof(double) * nall, cudaMemcpyHostToDevice), message,
                 message_size, "copy z"))
    return 1;
  if (device_neighbors)
  {
    if (device_nn_radial == nullptr || device_nn_angular == nullptr
        || (max_radial > 0 && device_nl_radial == nullptr)
        || (max_angular > 0 && device_nl_angular == nullptr))
      return set_message(message, message_size, "device neighbor-list pointer is null");
    d_nn_radial = device_nn_radial;
    d_nl_radial = device_nl_radial;
    d_nn_angular = device_nn_angular;
    d_nl_angular = device_nl_angular;
  }
  else
  {
    if (check_cuda(cudaMemcpy(model->d_nn_radial, nn_radial, sizeof(int) * nall,
                              cudaMemcpyHostToDevice),
                   message, message_size, "copy nn_radial"))
      return 1;
    if (check_cuda(cudaMemcpy(model->d_nn_angular, nn_angular, sizeof(int) * nall,
                              cudaMemcpyHostToDevice),
                   message, message_size, "copy nn_angular"))
      return 1;
    if (check_cuda(cudaMemcpy(model->d_nl_radial, nl_radial, sizeof(int) * nall * max_radial,
                              cudaMemcpyHostToDevice),
                   message, message_size, "copy nl_radial"))
      return 1;
    if (check_cuda(cudaMemcpy(model->d_nl_angular, nl_angular, sizeof(int) * nall * max_angular,
                              cudaMemcpyHostToDevice),
                   message, message_size, "copy nl_angular"))
      return 1;
  }

  cudaMemset(model->d_pe, 0, sizeof(double) * nall);
  cudaMemset(model->d_fx, 0, sizeof(double) * nall);
  cudaMemset(model->d_fy, 0, sizeof(double) * nall);
  cudaMemset(model->d_fz, 0, sizeof(double) * nall);
  cudaMemset(model->d_virial, 0, sizeof(double) * nall * 9);
  cudaMemset(model->d_f12x, 0, sizeof(float) * nall * max_angular);
  cudaMemset(model->d_f12y, 0, sizeof(float) * nall * max_angular);
  cudaMemset(model->d_f12z, 0, sizeof(float) * nall * max_angular);

  constexpr int block_size = 64;
  const int grid_all = (nall + block_size - 1) / block_size;
  const int grid_local = (nlocal + block_size - 1) / block_size;
  find_descriptor_lammps<<<grid_all, block_size>>>(
      model->paramb, model->annmb, nall, d_nn_radial, d_nl_radial,
      d_nn_angular, d_nl_angular, model->d_type, model->d_x, model->d_y, model->d_z,
      model->d_pe, model->d_Fp, model->d_sum_fxyz);
  if (check_cuda(cudaGetLastError(), message, message_size, "find_descriptor_lammps")) return 1;

  find_force_radial_lammps<<<grid_local, block_size>>>(
      model->paramb, model->annmb, nall, nlocal, d_nn_radial, d_nl_radial,
      model->d_type, model->d_x, model->d_y, model->d_z, model->d_Fp, model->d_fx, model->d_fy,
      model->d_fz, model->d_virial);
  if (check_cuda(cudaGetLastError(), message, message_size, "find_force_radial_lammps")) return 1;

  find_partial_force_angular_lammps<<<grid_all, block_size>>>(
      model->paramb, model->annmb, nall, d_nn_angular, d_nl_angular, model->d_type,
      model->d_x, model->d_y, model->d_z, model->d_Fp, model->d_sum_fxyz, model->d_f12x,
      model->d_f12y, model->d_f12z);
  if (check_cuda(cudaGetLastError(), message, message_size, "find_partial_force_angular_lammps"))
    return 1;

  find_force_many_body_lammps<<<grid_local, block_size>>>(
      nall, nlocal, d_nn_angular, d_nl_angular, model->d_f12x, model->d_f12y,
      model->d_f12z, model->d_x, model->d_y, model->d_z, model->d_fx, model->d_fy, model->d_fz,
      model->d_virial);
  if (check_cuda(cudaGetLastError(), message, message_size, "find_force_many_body_lammps"))
    return 1;
  if (check_cuda(cudaDeviceSynchronize(), message, message_size, "cudaDeviceSynchronize")) return 1;

  if (check_cuda(cudaMemcpy(result.potential, model->d_pe, sizeof(double) * nlocal,
                            cudaMemcpyDeviceToHost),
                 message, message_size, "copy potential"))
    return 1;
  if (check_cuda(cudaMemcpy(result.force, model->d_fx, sizeof(double) * nlocal,
                            cudaMemcpyDeviceToHost),
                 message, message_size, "copy fx"))
    return 1;
  if (check_cuda(cudaMemcpy(result.force + nlocal, model->d_fy, sizeof(double) * nlocal,
                            cudaMemcpyDeviceToHost),
                 message, message_size, "copy fy"))
    return 1;
  if (check_cuda(cudaMemcpy(result.force + 2 * nlocal, model->d_fz, sizeof(double) * nlocal,
                            cudaMemcpyDeviceToHost),
                 message, message_size, "copy fz"))
    return 1;
  for (int k = 0; k < 9; ++k) {
    if (check_cuda(cudaMemcpy(result.virial + k * nlocal, model->d_virial + k * nall,
                              sizeof(double) * nlocal, cudaMemcpyDeviceToHost),
                   message, message_size, "copy virial"))
      return 1;
  }
  return 0;
}

int nep_gpu_compute(NEP_GPU_Model *model, int nlocal, int nall, int max_radial,
                    int max_angular, const int *type, const double *x, const double *y,
                    const double *z, const int *nn_radial, const int *nl_radial,
                    const int *nn_angular, const int *nl_angular,
                    NEP_GPU_Result result, char *message, int message_size)
{
  return nep_gpu_compute_impl(model, nlocal, nall, max_radial, max_angular, type, x, y, z,
                              nn_radial, nl_radial, nn_angular, nl_angular, false, nullptr,
                              nullptr, nullptr, nullptr, result, message, message_size);
}

int nep_gpu_compute_device_neighbors(NEP_GPU_Model *model, int nlocal, int nall, int max_radial,
                                     int max_angular, const int *type, const double *x,
                                     const double *y, const double *z,
                                     const void *device_nn_radial, const void *device_nl_radial,
                                     const void *device_nn_angular, const void *device_nl_angular,
                                     NEP_GPU_Result result, char *message, int message_size)
{
  return nep_gpu_compute_impl(model, nlocal, nall, max_radial, max_angular, type, x, y, z,
                              nullptr, nullptr, nullptr, nullptr, true,
                              static_cast<const int*>(device_nn_radial),
                              static_cast<const int*>(device_nl_radial),
                              static_cast<const int*>(device_nn_angular),
                              static_cast<const int*>(device_nl_angular), result, message,
                              message_size);
}
