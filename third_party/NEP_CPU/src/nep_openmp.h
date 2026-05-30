/*
 * OpenMP helpers for NEP force accumulation.
 * Each thread writes to private per-atom buffers; a parallel reduction merges into global arrays.
 */
#pragma once

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <vector>

namespace nep_openmp
{

struct ForceBuffers {
  std::vector<std::vector<double>> fx;
  std::vector<std::vector<double>> fy;
  std::vector<std::vector<double>> fz;
  std::vector<std::vector<double>> virial;
  std::vector<std::vector<double>> pe;
  int nthreads = 1;

  void init(const int N, const bool with_virial, const bool with_pe)
  {
#if defined(_OPENMP)
    nthreads = omp_get_max_threads();
#else
    nthreads = 1;
#endif
    fx.assign(nthreads, std::vector<double>(N, 0.0));
    fy.assign(nthreads, std::vector<double>(N, 0.0));
    fz.assign(nthreads, std::vector<double>(N, 0.0));
    if (with_virial) {
      virial.assign(nthreads, std::vector<double>(N * 9, 0.0));
    } else {
      virial.clear();
    }
    if (with_pe) {
      pe.assign(nthreads, std::vector<double>(N, 0.0));
    } else {
      pe.clear();
    }
  }

  double* fx_ptr(const int tid) { return fx[tid].data(); }
  double* fy_ptr(const int tid) { return fy[tid].data(); }
  double* fz_ptr(const int tid) { return fz[tid].data(); }
  double* virial_ptr(const int tid) { return virial.empty() ? nullptr : virial[tid].data(); }
  double* pe_ptr(const int tid) { return pe.empty() ? nullptr : pe[tid].data(); }
};

inline void accumulate_force(
  const int n1,
  const int n2,
  const double f12[3],
  double* fx,
  double* fy,
  double* fz)
{
  if (fx) {
    fx[n1] += f12[0];
    fx[n2] -= f12[0];
  }
  if (fy) {
    fy[n1] += f12[1];
    fy[n2] -= f12[1];
  }
  if (fz) {
    fz[n1] += f12[2];
    fz[n2] -= f12[2];
  }
}

inline void accumulate_virial(
  const int n2,
  const int N,
  const double r12[3],
  const double f12[3],
  double* virial)
{
  virial[n2 + 0 * N] -= r12[0] * f12[0];
  virial[n2 + 1 * N] -= r12[0] * f12[1];
  virial[n2 + 2 * N] -= r12[0] * f12[2];
  virial[n2 + 3 * N] -= r12[1] * f12[0];
  virial[n2 + 4 * N] -= r12[1] * f12[1];
  virial[n2 + 5 * N] -= r12[1] * f12[2];
  virial[n2 + 6 * N] -= r12[2] * f12[0];
  virial[n2 + 7 * N] -= r12[2] * f12[1];
  virial[n2 + 8 * N] -= r12[2] * f12[2];
}

inline void accumulate_dipole_virial(
  const int n2,
  const int N,
  const double r12[3],
  const double f12[3],
  double* virial)
{
  const double r12_square = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
  virial[n2 + 0 * N] -= r12_square * f12[0];
  virial[n2 + 1 * N] -= r12_square * f12[1];
  virial[n2 + 2 * N] -= r12_square * f12[2];
}

inline void reduce_forces(
  const ForceBuffers& buf,
  const int N,
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial,
  double* g_pe)
{
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < N; ++i) {
    if (g_fx) {
      double sx = 0.0;
      for (int t = 0; t < buf.nthreads; ++t) {
        sx += buf.fx[t][i];
      }
      g_fx[i] += sx;
    }
    if (g_fy) {
      double sy = 0.0;
      for (int t = 0; t < buf.nthreads; ++t) {
        sy += buf.fy[t][i];
      }
      g_fy[i] += sy;
    }
    if (g_fz) {
      double sz = 0.0;
      for (int t = 0; t < buf.nthreads; ++t) {
        sz += buf.fz[t][i];
      }
      g_fz[i] += sz;
    }
    if (g_pe && !buf.pe.empty()) {
      double sp = 0.0;
      for (int t = 0; t < buf.nthreads; ++t) {
        sp += buf.pe[t][i];
      }
      g_pe[i] += sp;
    }
  }

  if (g_virial && !buf.virial.empty()) {
    const int nvir = N * 9;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < nvir; ++i) {
      double sv = 0.0;
      for (int t = 0; t < buf.nthreads; ++t) {
        sv += buf.virial[t][i];
      }
      g_virial[i] += sv;
    }
  }
}

} // namespace nep_openmp
