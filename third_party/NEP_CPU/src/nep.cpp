/*
    Copyright 2022 Zheyong Fan, Junjie Wang, Eric Lindgren
    This file is part of NEP_CPU.
    NEP_CPU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    NEP_CPU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with NEP_CPU.  If not, see <http://www.gnu.org/licenses/>.
*/

/*----------------------------------------------------------------------------80
A CPU implementation of the neuroevolution potential (NEP)
Ref: Zheyong Fan et al., Neuroevolution machine learning potentials:
Combining high accuracy and low cost in atomistic simulations and application to
heat transport, Phys. Rev. B. 104, 104309 (2021).
------------------------------------------------------------------------------*/

#include "nep.h"
#include "dftd3para.h"
#include "nep_openmp.h"
#include "nep_utilities.h"
#include "neighbor_nep.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace
{

void find_descriptor_small_box(
  const bool calculating_potential,
  const bool calculating_descriptor,
  const bool calculating_latent_space,
  const bool calculating_polarizability,
  NEP::ParaMB& paramb,
  NEP::ANN& annmb,
  const int N,
  const int* g_NN_radial,
  const int* g_NL_radial,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12_radial,
  const double* g_y12_radial,
  const double* g_z12_radial,
  const double* g_x12_angular,
  const double* g_y12_angular,
  const double* g_z12_angular,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gn_radial,
  const double* g_gn_angular,
#endif
  double* g_Fp,
  double* g_sum_fxyz,
  double* g_potential,
  double* g_descriptor,
  double* g_latent_space,
  double* g_virial,
  bool calculating_B_projection,
  double* g_B_projection)
{
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int n1 = 0; n1 < N; ++n1) {
    const int t1 = g_type[n1];
    double q[MAX_DIM] = {0.0};

    for (int i1 = 0; i1 < g_NN_radial[n1]; ++i1) {
      const int index = i1 * N + n1;
      const int n2 = g_NL_radial[index];
      const int t2 = g_type[n2];
      double r12[3] = {g_x12_radial[index], g_y12_radial[index], g_z12_radial[index]};
      const double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);

#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
      int index_left, index_right;
      double weight_left, weight_right;
      find_index_and_weight(
        d12 * paramb.rcinv_radial, index_left, index_right, weight_left, weight_right);
      const int t12 = t1 * paramb.num_types + t2;
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        q[n] +=
          g_gn_radial[(index_left * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
            weight_left +
          g_gn_radial[(index_right * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
            weight_right;
      }
#else
      double fc12;
      double rc = (paramb.rc_radial[t1] + paramb.rc_radial[t2]) * 0.5;
      double rcinv = 1.0 / rc;
      find_fc(rc, rcinv, d12, fc12);
      double fn12[MAX_NUM_N];
      find_fn(paramb.basis_size_radial, rcinv, d12, fc12, fn12);
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        double gn12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_radial; ++k) {
          int c_index = (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2;
          gn12 += fn12[k] * annmb.c[c_index];
        }
        q[n] += gn12;
      }
#endif
    }

    for (int n = 0; n <= paramb.n_max_angular; ++n) {
      double s[NUM_OF_ABC] = {0.0};
      for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
        int index = i1 * N + n1;
        int n2 = g_NL_angular[index];
        double r12[3] = {g_x12_angular[index], g_y12_angular[index], g_z12_angular[index]};
        double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
        int index_left, index_right;
        double weight_left, weight_right;
        find_index_and_weight(
          d12 * paramb.rcinv_angular, index_left, index_right, weight_left, weight_right);
        int t12 = t1 * paramb.num_types + g_type[n2];
        double gn12 =
          g_gn_angular[(index_left * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n] *
            weight_left +
          g_gn_angular[(index_right * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n] *
            weight_right;
        accumulate_s(paramb.L_max, d12, r12[0], r12[1], r12[2], gn12, s);
#else
        int t2 = g_type[n2];
        double fc12;
        double rc = (paramb.rc_angular[t1] + paramb.rc_angular[t2]) * 0.5;
        double rcinv = 1.0 / rc;
        find_fc(rc, rcinv, d12, fc12);
        double fn12[MAX_NUM_N];
        find_fn(paramb.basis_size_angular, rcinv, d12, fc12, fn12);
        double gn12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_angular; ++k) {
          int c_index = (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2 + paramb.num_c_radial;
          gn12 += fn12[k] * annmb.c[c_index];
        }
        accumulate_s(paramb.L_max, d12, r12[0], r12[1], r12[2], gn12, s);
#endif
      }
      find_q(
        paramb.L_max, paramb.num_L, paramb.n_max_angular + 1, n, s, q + (paramb.n_max_radial + 1));
      for (int abc = 0; abc < NUM_OF_ABC; ++abc) {
        g_sum_fxyz[(n * NUM_OF_ABC + abc) * N + n1] = s[abc];
      }
    }

    if (calculating_descriptor) {
      for (int d = 0; d < annmb.dim; ++d) {
        g_descriptor[d * N + n1] = q[d] * paramb.q_scaler[d];
      }
    }

    if (
      calculating_potential || calculating_latent_space || calculating_polarizability ||
      calculating_B_projection) {
      for (int d = 0; d < annmb.dim; ++d) {
        q[d] = q[d] * paramb.q_scaler[d];
      }

      double F = 0.0, Fp[MAX_DIM] = {0.0}, latent_space[MAX_NEURON] = {0.0};

      if (calculating_polarizability) {
        apply_ann_one_layer(
          annmb.dim, annmb.num_neurons1, annmb.w0_pol[t1], annmb.b0_pol[t1], annmb.w1_pol[t1],
          annmb.b1_pol, q, F, Fp, latent_space, false, nullptr);
        g_virial[n1] = F;
        g_virial[n1 + N * 4] = F;
        g_virial[n1 + N * 8] = F;

        for (int d = 0; d < annmb.dim; ++d) {
          Fp[d] = 0.0;
        }
        for (int d = 0; d < annmb.num_neurons1; ++d) {
          latent_space[d] = 0.0;
        }
      }

      if (paramb.version == 5) {
        apply_ann_one_layer_nep5(
          annmb.dim, annmb.num_neurons1, annmb.w0[t1], annmb.b0[t1], annmb.w1[t1], annmb.b1, q, F,
          Fp, latent_space);
      } else {
        apply_ann_one_layer(
          annmb.dim, annmb.num_neurons1, annmb.w0[t1], annmb.b0[t1], annmb.w1[t1], annmb.b1, q, F,
          Fp, latent_space, calculating_B_projection,
          g_B_projection + n1 * (annmb.num_neurons1 * (annmb.dim + 2)));
      }

      if (calculating_latent_space) {
        for (int n = 0; n < annmb.num_neurons1; ++n) {
          g_latent_space[n * N + n1] = latent_space[n];
        }
      }

      if (calculating_potential) {
        g_potential[n1] += F;
      }

      for (int d = 0; d < annmb.dim; ++d) {
        g_Fp[d * N + n1] = Fp[d] * paramb.q_scaler[d];
      }
    }
  }
}

void find_force_radial_small_box(
  const bool is_dipole,
  NEP::ParaMB& paramb,
  NEP::ANN& annmb,
  const int N,
  const int* g_NN,
  const int* g_NL,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  const double* g_Fp,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gnp_radial,
#endif
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial)
{
  nep_openmp::ForceBuffers buf;
  buf.init(N, g_virial != nullptr, false);

#if defined(_OPENMP)
#pragma omp parallel
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    double* fx = buf.fx_ptr(tid);
    double* fy = buf.fy_ptr(tid);
    double* fz = buf.fz_ptr(tid);
    double* virial = buf.virial_ptr(tid);

#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
    for (int n1 = 0; n1 < N; ++n1) {
      const int t1 = g_type[n1];
      for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
        const int index = i1 * N + n1;
        const int n2 = g_NL[index];
        const int t2 = g_type[n2];
        double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
        const double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
        const double d12inv = 1.0 / d12;
        double f12[3] = {0.0};
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
        int index_left, index_right;
        double weight_left, weight_right;
        find_index_and_weight(
          d12 * paramb.rcinv_radial, index_left, index_right, weight_left, weight_right);
        const int t12 = t1 * paramb.num_types + t2;
        for (int n = 0; n <= paramb.n_max_radial; ++n) {
          const double gnp12 =
            g_gnp_radial[(index_left * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
              weight_left +
            g_gnp_radial[(index_right * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
              weight_right;
          const double tmp12 = g_Fp[n1 + n * N] * gnp12 * d12inv;
          for (int d = 0; d < 3; ++d) {
            f12[d] += tmp12 * r12[d];
          }
        }
#else
        double fc12, fcp12;
        const double rc = (paramb.rc_radial[t1] + paramb.rc_radial[t2]) * 0.5;
        const double rcinv = 1.0 / rc;
        find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);
        double fn12[MAX_NUM_N];
        double fnp12[MAX_NUM_N];
        find_fn_and_fnp(paramb.basis_size_radial, rcinv, d12, fc12, fcp12, fn12, fnp12);
        const int c_base = t1 * paramb.num_types + t2;
        for (int n = 0; n <= paramb.n_max_radial; ++n) {
          double gnp12 = 0.0;
          for (int k = 0; k <= paramb.basis_size_radial; ++k) {
            const int c_index =
              (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq + c_base;
            gnp12 += fnp12[k] * annmb.c[c_index];
          }
          const double tmp12 = g_Fp[n1 + n * N] * gnp12 * d12inv;
          for (int d = 0; d < 3; ++d) {
            f12[d] += tmp12 * r12[d];
          }
        }
#endif

        nep_openmp::accumulate_force(n1, n2, f12, fx, fy, fz);
        if (virial) {
          if (!is_dipole) {
            nep_openmp::accumulate_virial(n2, N, r12, f12, virial);
          } else {
            nep_openmp::accumulate_dipole_virial(n2, N, r12, f12, virial);
          }
        }
      }
    }
  }

  nep_openmp::reduce_forces(buf, N, g_fx, g_fy, g_fz, g_virial, nullptr);
}

void find_force_angular_small_box(
  const bool is_dipole,
  NEP::ParaMB& paramb,
  NEP::ANN& annmb,
  const int N,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  const double* g_Fp,
  const double* g_sum_fxyz,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gn_angular,
  const double* g_gnp_angular,
#endif
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial)
{
  nep_openmp::ForceBuffers buf;
  buf.init(N, g_virial != nullptr, false);

#if defined(_OPENMP)
#pragma omp parallel
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    double* fx = buf.fx_ptr(tid);
    double* fy = buf.fy_ptr(tid);
    double* fz = buf.fz_ptr(tid);
    double* virial = buf.virial_ptr(tid);

#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
    for (int n1 = 0; n1 < N; ++n1) {
      double Fp[MAX_DIM_ANGULAR] = {0.0};
      double sum_fxyz[NUM_OF_ABC * MAX_NUM_N];
      for (int d = 0; d < paramb.dim_angular; ++d) {
        Fp[d] = g_Fp[(paramb.n_max_radial + 1 + d) * N + n1];
      }
      for (int d = 0; d < (paramb.n_max_angular + 1) * NUM_OF_ABC; ++d) {
        sum_fxyz[d] = g_sum_fxyz[d * N + n1];
      }

      const int t1 = g_type[n1];

      for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
        const int index = i1 * N + n1;
        const int n2 = g_NL_angular[n1 + N * i1];
        double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
        const double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
        double f12[3] = {0.0};
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
        int index_left, index_right;
        double weight_left, weight_right;
        find_index_and_weight(
          d12 * paramb.rcinv_angular, index_left, index_right, weight_left, weight_right);
        const int t12 = t1 * paramb.num_types + g_type[n2];
        for (int n = 0; n <= paramb.n_max_angular; ++n) {
          const int index_left_all =
            (index_left * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n;
          const int index_right_all =
            (index_right * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n;
          const double gn12 =
            g_gn_angular[index_left_all] * weight_left + g_gn_angular[index_right_all] * weight_right;
          const double gnp12 = g_gnp_angular[index_left_all] * weight_left +
                               g_gnp_angular[index_right_all] * weight_right;
          accumulate_f12(
            paramb.L_max, paramb.num_L, n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp,
            sum_fxyz, f12);
        }
#else
        const int t2 = g_type[n2];
        double fc12, fcp12;
        const double rc = (paramb.rc_angular[t1] + paramb.rc_angular[t2]) * 0.5;
        const double rcinv = 1.0 / rc;
        find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);

        double fn12[MAX_NUM_N];
        double fnp12[MAX_NUM_N];
        find_fn_and_fnp(paramb.basis_size_angular, rcinv, d12, fc12, fcp12, fn12, fnp12);
        const int c_base = t1 * paramb.num_types + t2 + paramb.num_c_radial;
        for (int n = 0; n <= paramb.n_max_angular; ++n) {
          double gn12 = 0.0;
          double gnp12 = 0.0;
          for (int k = 0; k <= paramb.basis_size_angular; ++k) {
            const int c_index =
              (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq + c_base;
            gn12 += fn12[k] * annmb.c[c_index];
            gnp12 += fnp12[k] * annmb.c[c_index];
          }
          accumulate_f12(
            paramb.L_max, paramb.num_L, n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp,
            sum_fxyz, f12);
        }
#endif

        nep_openmp::accumulate_force(n1, n2, f12, fx, fy, fz);
        if (virial) {
          if (!is_dipole) {
            nep_openmp::accumulate_virial(n2, N, r12, f12, virial);
          } else {
            nep_openmp::accumulate_dipole_virial(n2, N, r12, f12, virial);
          }
        }
      }
    }
  }

  nep_openmp::reduce_forces(buf, N, g_fx, g_fy, g_fz, g_virial, nullptr);
}

void find_force_ZBL_small_box(
  const int N,
  NEP::ParaMB& paramb,
  const NEP::ZBL& zbl,
  const int* g_NN,
  const int* g_NL,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial,
  double* g_pe)
{
  nep_openmp::ForceBuffers buf;
  buf.init(N, g_virial != nullptr, g_pe != nullptr);

#if defined(_OPENMP)
#pragma omp parallel
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    double* fx = buf.fx_ptr(tid);
    double* fy = buf.fy_ptr(tid);
    double* fz = buf.fz_ptr(tid);
    double* virial = buf.virial_ptr(tid);
    double* pe = buf.pe_ptr(tid);

#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
    for (int n1 = 0; n1 < N; ++n1) {
      const int type1 = g_type[n1];
      const int zi = paramb.atomic_numbers[type1] + 1;
      const double pow_zi = pow(double(zi), 0.23);
      for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
        const int index = i1 * N + n1;
        const int n2 = g_NL[index];
        double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
        const double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
        const double d12inv = 1.0 / d12;
        double f, fp;
        const int type2 = g_type[n2];
        const int zj = paramb.atomic_numbers[type2] + 1;
        const double a_inv = (pow_zi + pow(double(zj), 0.23)) * 2.134563;
        const double zizj = K_C_SP * zi * zj;
        if (zbl.flexibled) {
          int t1, t2;
          if (type1 < type2) {
            t1 = type1;
            t2 = type2;
          } else {
            t1 = type2;
            t2 = type1;
          }
          const int zbl_index = t1 * zbl.num_types - (t1 * (t1 - 1)) / 2 + (t2 - t1);
          double ZBL_para[10];
          for (int i = 0; i < 10; ++i) {
            ZBL_para[i] = zbl.para[10 * zbl_index + i];
          }
          find_f_and_fp_zbl(ZBL_para, zizj, a_inv, d12, d12inv, f, fp);
        } else {
          double rc_inner = zbl.rc_inner;
          double rc_outer = zbl.rc_outer;
          if (paramb.use_typewise_cutoff_zbl) {
            rc_outer = std::min(
              (COVALENT_RADIUS[zi - 1] + COVALENT_RADIUS[zj - 1]) * paramb.typewise_cutoff_zbl_factor,
              rc_outer);
            rc_inner = 0.0;
          }
          find_f_and_fp_zbl(zizj, a_inv, rc_inner, rc_outer, d12, d12inv, f, fp);
        }
        const double f2 = fp * d12inv * 0.5;
        const double f12[3] = {r12[0] * f2, r12[1] * f2, r12[2] * f2};
        nep_openmp::accumulate_force(n1, n2, f12, fx, fy, fz);
        if (virial) {
          nep_openmp::accumulate_virial(n2, N, r12, f12, virial);
        }
        if (pe) {
          pe[n1] += f * 0.5;
        }
      }
    }
  }

  nep_openmp::reduce_forces(buf, N, g_fx, g_fy, g_fz, g_virial, g_pe);
}

void find_descriptor_small_box(
  const bool calculating_potential,
  const bool calculating_descriptor,
  NEP::ParaMB& paramb,
  NEP::ANN& annmb,
  const int N,
  const int* g_NN_radial,
  const int* g_NL_radial,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12_radial,
  const double* g_y12_radial,
  const double* g_z12_radial,
  const double* g_x12_angular,
  const double* g_y12_angular,
  const double* g_z12_angular,
  double* g_Fp,
  double* g_sum_fxyz,
  double* g_charge,
  double* g_charge_derivative,
  double* g_potential,
  double* g_descriptor)
{
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int n1 = 0; n1 < N; ++n1) {
    int t1 = g_type[n1];
    double q[MAX_DIM] = {0.0};

    for (int i1 = 0; i1 < g_NN_radial[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL_radial[index];
      double r12[3] = {g_x12_radial[index], g_y12_radial[index], g_z12_radial[index]};
      double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);

      double fc12;
      int t2 = g_type[n2];
      double rc = paramb.rc_radial_max;
      double rcinv = 1.0 / rc;
      find_fc(rc, rcinv, d12, fc12);
      double fn12[MAX_NUM_N];
      find_fn(paramb.basis_size_radial, rcinv, d12, fc12, fn12);
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        double gn12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_radial; ++k) {
          int c_index = (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2;
          gn12 += fn12[k] * annmb.c[c_index];
        }
        q[n] += gn12;
      }
    }

    for (int n = 0; n <= paramb.n_max_angular; ++n) {
      double s[NUM_OF_ABC] = {0.0};
      for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
        int index = i1 * N + n1;
        int n2 = g_NL_angular[index];
        double r12[3] = {g_x12_angular[index], g_y12_angular[index], g_z12_angular[index]};
        double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
        int t2 = g_type[n2];
        double fc12;
        double rc = paramb.rc_angular_max;
        double rcinv = 1.0 / rc;
        find_fc(rc, rcinv, d12, fc12);
        double fn12[MAX_NUM_N];
        find_fn(paramb.basis_size_angular, rcinv, d12, fc12, fn12);
        double gn12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_angular; ++k) {
          int c_index = (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2 + paramb.num_c_radial;
          gn12 += fn12[k] * annmb.c[c_index];
        }
        accumulate_s(paramb.L_max, d12, r12[0], r12[1], r12[2], gn12, s);
      }
      find_q(
        paramb.L_max, paramb.num_L, paramb.n_max_angular + 1, n, s, q + (paramb.n_max_radial + 1));
      for (int abc = 0; abc < NUM_OF_ABC; ++abc) {
        g_sum_fxyz[(n * NUM_OF_ABC + abc) * N + n1] = s[abc];
      }
    }

    if (calculating_descriptor) {
      for (int d = 0; d < annmb.dim; ++d) {
        g_descriptor[d * N + n1] = q[d] * paramb.q_scaler[d];
      }
    }

    if (calculating_potential) {
      for (int d = 0; d < annmb.dim; ++d) {
        q[d] = q[d] * paramb.q_scaler[d];
      }

      double F = 0.0, Fp[MAX_DIM] = {0.0};
      double charge = 0.0;
      double charge_derivative[MAX_DIM] = {0.0};

      apply_ann_one_layer_charge(
        annmb.dim,
        annmb.num_neurons1,
        annmb.w0[t1],
        annmb.b0[t1],
        annmb.w1[t1],
        annmb.b1,
        q,
        F,
        Fp,
        charge,
        charge_derivative);

      if (calculating_potential) {
        g_potential[n1] += F;
        g_charge[n1] = charge;
      }

      for (int d = 0; d < annmb.dim; ++d) {
        g_Fp[d * N + n1] = Fp[d] * paramb.q_scaler[d];
        g_charge_derivative[d * N + n1] = charge_derivative[d] * paramb.q_scaler[d];
      }
    }
  }
}

void zero_total_charge(const int N, double* g_charge)
{
  double mean_charge = 0.0;
  for (int n = 0; n < N; ++n) {
    mean_charge += g_charge[n];
  }
  mean_charge /= N;
  for (int n = 0; n < N; ++n) {
    g_charge[n] -= mean_charge;
  }
}

void find_force_radial_small_box(
  NEP::ParaMB& paramb,
  NEP::ANN& annmb,
  const int N,
  const int* g_NN,
  const int* g_NL,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  const double* g_Fp,
  const double* g_charge_derivative,
  const double* g_D_real, 
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial)
{
  for (int n1 = 0; n1 < N; ++n1) {
    int t1 = g_type[n1];
    for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL[index];
      int t2 = g_type[n2];
      double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      double d12inv = 1.0 / d12;
      double f12[3] = {0.0};
      double fc12, fcp12;
      double rc = paramb.rc_radial_max;
      double rcinv = 1.0 / rc;
      find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);
      double fn12[MAX_NUM_N];
      double fnp12[MAX_NUM_N];
      find_fn_and_fnp(paramb.basis_size_radial, rcinv, d12, fc12, fcp12, fn12, fnp12);
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        double gnp12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_radial; ++k) {
          int c_index = (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2;
          gnp12 += fnp12[k] * annmb.c[c_index];
        }
        double tmp12 = (g_Fp[n1 + n * N] + g_charge_derivative[n1 + n * N] * g_D_real[n1]) * gnp12 * d12inv;
        for (int d = 0; d < 3; ++d) {
          f12[d] += tmp12 * r12[d];
        }
      }

      if (g_fx) {
        g_fx[n1] += f12[0];
        g_fx[n2] -= f12[0];
      }

      if (g_fy) {
        g_fy[n1] += f12[1];
        g_fy[n2] -= f12[1];
      }

      if (g_fz) {
        g_fz[n1] += f12[2];
        g_fz[n2] -= f12[2];
      }

      g_virial[n2 + 0 * N] -= r12[0] * f12[0];
      g_virial[n2 + 1 * N] -= r12[0] * f12[1];
      g_virial[n2 + 2 * N] -= r12[0] * f12[2];
      g_virial[n2 + 3 * N] -= r12[1] * f12[0];
      g_virial[n2 + 4 * N] -= r12[1] * f12[1];
      g_virial[n2 + 5 * N] -= r12[1] * f12[2];
      g_virial[n2 + 6 * N] -= r12[2] * f12[0];
      g_virial[n2 + 7 * N] -= r12[2] * f12[1];
      g_virial[n2 + 8 * N] -= r12[2] * f12[2];
    }
  }
}

void find_force_angular_small_box(
  NEP::ParaMB& paramb,
  NEP::ANN& annmb,
  const int N,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  const double* g_Fp,
  const double* g_charge_derivative,
  const double* g_D_real, 
  const double* g_sum_fxyz,
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial)
{
  for (int n1 = 0; n1 < N; ++n1) {

    double Fp[MAX_DIM_ANGULAR] = {0.0};
    double sum_fxyz[NUM_OF_ABC * MAX_NUM_N];
    for (int d = 0; d < paramb.dim_angular; ++d) {
      Fp[d] = g_Fp[(paramb.n_max_radial + 1 + d) * N + n1] 
        + g_charge_derivative[(paramb.n_max_radial + 1 + d) * N + n1] * g_D_real[n1];
    }
    for (int d = 0; d < (paramb.n_max_angular + 1) * NUM_OF_ABC; ++d) {
      sum_fxyz[d] = g_sum_fxyz[d * N + n1];
    }

    int t1 = g_type[n1];

    for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL_angular[n1 + N * i1];
      double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      double f12[3] = {0.0};
      int t2 = g_type[n2];
      double fc12, fcp12;
      double rc = paramb.rc_angular_max;
      double rcinv = 1.0 / rc;
      find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);

      double fn12[MAX_NUM_N];
      double fnp12[MAX_NUM_N];
      find_fn_and_fnp(paramb.basis_size_angular, rcinv, d12, fc12, fcp12, fn12, fnp12);
      for (int n = 0; n <= paramb.n_max_angular; ++n) {
        double gn12 = 0.0;
        double gnp12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_angular; ++k) {
          int c_index = (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2 + paramb.num_c_radial;
          gn12 += fn12[k] * annmb.c[c_index];
          gnp12 += fnp12[k] * annmb.c[c_index];
        }
        accumulate_f12(
          paramb.L_max, paramb.num_L, n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp,
          sum_fxyz, f12);
      }

      if (g_fx) {
        g_fx[n1] += f12[0];
        g_fx[n2] -= f12[0];
      }

      if (g_fy) {
        g_fy[n1] += f12[1];
        g_fy[n2] -= f12[1];
      }

      if (g_fz) {
        g_fz[n1] += f12[2];
        g_fz[n2] -= f12[2];
      }

      g_virial[n2 + 0 * N] -= r12[0] * f12[0];
      g_virial[n2 + 1 * N] -= r12[0] * f12[1];
      g_virial[n2 + 2 * N] -= r12[0] * f12[2];
      g_virial[n2 + 3 * N] -= r12[1] * f12[0];
      g_virial[n2 + 4 * N] -= r12[1] * f12[1];
      g_virial[n2 + 5 * N] -= r12[1] * f12[2];
      g_virial[n2 + 6 * N] -= r12[2] * f12[0];
      g_virial[n2 + 7 * N] -= r12[2] * f12[1];
      g_virial[n2 + 8 * N] -= r12[2] * f12[2];
    }
  }
}

void find_bec_diagonal(const int N, const double* g_q, double* g_bec)
{
  for (int n1 = 0; n1 < N; ++n1) {
    g_bec[n1 + N * 0] = g_q[n1];
    g_bec[n1 + N * 1] = 0.0;
    g_bec[n1 + N * 2] = 0.0;
    g_bec[n1 + N * 3] = 0.0;
    g_bec[n1 + N * 4] = g_q[n1];
    g_bec[n1 + N * 5] = 0.0;
    g_bec[n1 + N * 6] = 0.0;
    g_bec[n1 + N * 7] = 0.0;
    g_bec[n1 + N * 8] = g_q[n1];
  }
}

void find_bec_radial_small_box(
  const NEP::ParaMB paramb,
  const NEP::ANN annmb,
  const int N,
  const int* g_NN,
  const int* g_NL,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  const double* g_charge_derivative,
  double* g_bec)
{
  for (int n1 = 0; n1 < N; ++n1) {
    int t1 = g_type[n1];
    for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL[index];
      int t2 = g_type[n2];
      double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      double d12inv = 1.0 / d12;
      double fc12, fcp12;
      double rc = paramb.rc_radial_max;
      double rcinv = 1.0 / rc;
      find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);
      double fn12[MAX_NUM_N];
      double fnp12[MAX_NUM_N];
      double f12[3] = {0.0};

      find_fn_and_fnp(paramb.basis_size_radial, rcinv, d12, fc12, fcp12, fn12, fnp12);
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        double gnp12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_radial; ++k) {
          int c_index = (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2;
          gnp12 += fnp12[k] * annmb.c[c_index];
        }
        const double tmp12 = g_charge_derivative[n1 + n * N] * gnp12 * d12inv;
        for (int d = 0; d < 3; ++d) {
          f12[d] += tmp12 * r12[d];
        }
      }

      double bec_xx = 0.5* (r12[0] * f12[0]);
      double bec_xy = 0.5* (r12[0] * f12[1]);
      double bec_xz = 0.5* (r12[0] * f12[2]);
      double bec_yx = 0.5* (r12[1] * f12[0]);
      double bec_yy = 0.5* (r12[1] * f12[1]);
      double bec_yz = 0.5* (r12[1] * f12[2]);
      double bec_zx = 0.5* (r12[2] * f12[0]);
      double bec_zy = 0.5* (r12[2] * f12[1]);
      double bec_zz = 0.5* (r12[2] * f12[2]);

      g_bec[n1] += bec_xx;
      g_bec[n1 + N] += bec_xy;
      g_bec[n1 + N * 2] += bec_xz;
      g_bec[n1 + N * 3] += bec_yx;
      g_bec[n1 + N * 4] += bec_yy;
      g_bec[n1 + N * 5] += bec_yz;
      g_bec[n1 + N * 6] += bec_zx;
      g_bec[n1 + N * 7] += bec_zy;
      g_bec[n1 + N * 8] += bec_zz;

      g_bec[n2] -= bec_xx;
      g_bec[n2 + N] -= bec_xy;
      g_bec[n2 + N * 2] -= bec_xz;
      g_bec[n2 + N * 3] -= bec_yx;
      g_bec[n2 + N * 4] -= bec_yy;
      g_bec[n2 + N * 5] -= bec_yz;
      g_bec[n2 + N * 6] -= bec_zx;
      g_bec[n2 + N * 7] -= bec_zy;
      g_bec[n2 + N * 8] -= bec_zz;
    }
  }
}

void find_bec_angular_small_box(
  NEP::ParaMB paramb,
  NEP::ANN annmb,
  const int N,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  const double* g_charge_derivative,
  const double* g_sum_fxyz,
  double* g_bec)
{
  for (int n1 = 0; n1 < N; ++n1) {
    double Fp[MAX_DIM_ANGULAR] = {0.0};
    double sum_fxyz[NUM_OF_ABC * MAX_NUM_N];
    for (int d = 0; d < paramb.dim_angular; ++d) {
      Fp[d] = g_charge_derivative[(paramb.n_max_radial + 1 + d) * N + n1];
    }
    for (int d = 0; d < (paramb.n_max_angular + 1) * NUM_OF_ABC; ++d) {
      sum_fxyz[d] = g_sum_fxyz[d * N + n1];
    }
    int t1 = g_type[n1];
    for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL_angular[index];
      double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      double f12[3] = {0.0};
      double fc12, fcp12;
      int t2 = g_type[n2];
      double rc = paramb.rc_angular_max;
      double rcinv = 1.0 / rc;
      find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);

      double fn12[MAX_NUM_N];
      double fnp12[MAX_NUM_N];
      find_fn_and_fnp(paramb.basis_size_angular, rcinv, d12, fc12, fcp12, fn12, fnp12);
      for (int n = 0; n <= paramb.n_max_angular; ++n) {
        double gn12 = 0.0;
        double gnp12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_angular; ++k) {
          int c_index = (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2 + paramb.num_c_radial;
          gn12 += fn12[k] * annmb.c[c_index];
          gnp12 += fnp12[k] * annmb.c[c_index];
        }
        accumulate_f12(
          paramb.L_max,
          paramb.num_L,
          n,
          paramb.n_max_angular + 1,
          d12,
          r12,
          gn12,
          gnp12,
          Fp,
          sum_fxyz,
          f12);
      }

      double bec_xx = 0.5* (r12[0] * f12[0]);
      double bec_xy = 0.5* (r12[0] * f12[1]);
      double bec_xz = 0.5* (r12[0] * f12[2]);
      double bec_yx = 0.5* (r12[1] * f12[0]);
      double bec_yy = 0.5* (r12[1] * f12[1]);
      double bec_yz = 0.5* (r12[1] * f12[2]);
      double bec_zx = 0.5* (r12[2] * f12[0]);
      double bec_zy = 0.5* (r12[2] * f12[1]);
      double bec_zz = 0.5* (r12[2] * f12[2]);

      g_bec[n1] += bec_xx;
      g_bec[n1 + N] += bec_xy;
      g_bec[n1 + N * 2] += bec_xz;
      g_bec[n1 + N * 3] += bec_yx;
      g_bec[n1 + N * 4] += bec_yy;
      g_bec[n1 + N * 5] += bec_yz;
      g_bec[n1 + N * 6] += bec_zx;
      g_bec[n1 + N * 7] += bec_zy;
      g_bec[n1 + N * 8] += bec_zz;

      g_bec[n2] -= bec_xx;
      g_bec[n2 + N] -= bec_xy;
      g_bec[n2 + N * 2] -= bec_xz;
      g_bec[n2 + N * 3] -= bec_yx;
      g_bec[n2 + N * 4] -= bec_yy;
      g_bec[n2 + N * 5] -= bec_yz;
      g_bec[n2 + N * 6] -= bec_zx;
      g_bec[n2 + N * 7] -= bec_zy;
      g_bec[n2 + N * 8] -= bec_zz;
    }
  }
}

void scale_bec(const int N, const double* sqrt_epsilon_inf, double* g_bec)
{
  for (int n1 = 0; n1 < N; ++n1) {
    for (int d = 0; d < 9; ++d) {
      g_bec[n1 + N * d] *= sqrt_epsilon_inf[0];
    }
  }
}

void find_force_charge_real_space_only_small_box(
  const int N,
  const NEP::Charge_Para charge_para,
  const int* g_NN,
  const int* g_NL,
  const double* g_charge,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial,
  double* g_pe,
  double* g_D_real)
{
  for (int n1 = 0; n1 < N; ++n1) {
    double s_fx = 0.0;
    double s_fy = 0.0;
    double s_fz = 0.0;
    double s_sxx = 0.0;
    double s_sxy = 0.0;
    double s_sxz = 0.0;
    double s_syx = 0.0;
    double s_syy = 0.0;
    double s_syz = 0.0;
    double s_szx = 0.0;
    double s_szy = 0.0;
    double s_szz = 0.0;
    double q1 = g_charge[n1];
    double s_pe = 0; // no self energy
    double D_real = 0; // no self energy

    for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL[index];
      double q2 = g_charge[n2];
      double qq = q1 * q2;
      double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      double d12inv = 1.0 / d12;

      double erfc_r = erfc(charge_para.alpha * d12) * d12inv;
      D_real += q2 * (erfc_r + charge_para.A * d12 + charge_para.B);
      s_pe += 0.5 * qq * (erfc_r + charge_para.A * d12 + charge_para.B);
      double f2 = erfc_r + charge_para.two_alpha_over_sqrt_pi * exp(-charge_para.alpha * charge_para.alpha * d12 * d12);
      f2 = -0.5 * K_C_SP * qq * (f2 * d12inv * d12inv - charge_para.A * d12inv);
      double f12[3] = {r12[0] * f2, r12[1] * f2, r12[2] * f2};
      double f21[3] = {-r12[0] * f2, -r12[1] * f2, -r12[2] * f2};

      s_fx += f12[0] - f21[0];
      s_fy += f12[1] - f21[1];
      s_fz += f12[2] - f21[2];
      s_sxx -= r12[0] * f12[0];
      s_sxy -= r12[0] * f12[1];
      s_sxz -= r12[0] * f12[2];
      s_syx -= r12[1] * f12[0];
      s_syy -= r12[1] * f12[1];
      s_syz -= r12[1] * f12[2];
      s_szx -= r12[2] * f12[0];
      s_szy -= r12[2] * f12[1];
      s_szz -= r12[2] * f12[2];
    }
    g_fx[n1] += s_fx;
    g_fy[n1] += s_fy;
    g_fz[n1] += s_fz;
    g_virial[n1 + 0 * N] += s_sxx;
    g_virial[n1 + 1 * N] += s_sxy;
    g_virial[n1 + 2 * N] += s_sxz;
    g_virial[n1 + 3 * N] += s_syx;
    g_virial[n1 + 4 * N] += s_syy;
    g_virial[n1 + 5 * N] += s_syz;
    g_virial[n1 + 6 * N] += s_szx;
    g_virial[n1 + 7 * N] += s_szy;
    g_virial[n1 + 8 * N] += s_szz;
    g_D_real[n1] = K_C_SP * D_real;
    g_pe[n1] += K_C_SP * s_pe;
  }
}

void find_force_charge_real_space_small_box(
  const int N,
  const NEP::Charge_Para charge_para,
  const int* g_NN,
  const int* g_NL,
  const double* g_charge,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial,
  double* g_pe,
  double* g_D_real)
{
  for (int n1 = 0; n1 < N; ++n1) {
    double s_fx = 0.0;
    double s_fy = 0.0;
    double s_fz = 0.0;
    double s_sxx = 0.0;
    double s_sxy = 0.0;
    double s_sxz = 0.0;
    double s_syx = 0.0;
    double s_syy = 0.0;
    double s_syz = 0.0;
    double s_szx = 0.0;
    double s_szy = 0.0;
    double s_szz = 0.0;
    double q1 = g_charge[n1];
    double s_pe = -charge_para.two_alpha_over_sqrt_pi * 0.5 * q1 * q1; // self energy part
    double D_real = -q1 * charge_para.two_alpha_over_sqrt_pi; // self energy part

    for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL[index];
      double q2 = g_charge[n2];
      double qq = q1 * q2;
      double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      double d12inv = 1.0 / d12;

      double erfc_r = erfc(charge_para.alpha * d12) * d12inv;
      D_real += q2 * erfc_r;
      s_pe += 0.5 * qq * erfc_r;
      double f2 = erfc_r + charge_para.two_alpha_over_sqrt_pi * exp(-charge_para.alpha * charge_para.alpha * d12 * d12);
      f2 *= -0.5 * K_C_SP * qq * d12inv * d12inv;
      double f12[3] = {r12[0] * f2, r12[1] * f2, r12[2] * f2};
      double f21[3] = {-r12[0] * f2, -r12[1] * f2, -r12[2] * f2};

      s_fx += f12[0] - f21[0];
      s_fy += f12[1] - f21[1];
      s_fz += f12[2] - f21[2];
      s_sxx -= r12[0] * f12[0];
      s_sxy -= r12[0] * f12[1];
      s_sxz -= r12[0] * f12[2];
      s_syx -= r12[1] * f12[0];
      s_syy -= r12[1] * f12[1];
      s_syz -= r12[1] * f12[2];
      s_szx -= r12[2] * f12[0];
      s_szy -= r12[2] * f12[1];
      s_szz -= r12[2] * f12[2];
    }
    g_fx[n1] += s_fx;
    g_fy[n1] += s_fy;
    g_fz[n1] += s_fz;
    g_virial[n1 + 0 * N] += s_sxx;
    g_virial[n1 + 1 * N] += s_sxy;
    g_virial[n1 + 2 * N] += s_sxz;
    g_virial[n1 + 3 * N] += s_syx;
    g_virial[n1 + 4 * N] += s_syy;
    g_virial[n1 + 5 * N] += s_syz;
    g_virial[n1 + 6 * N] += s_szx;
    g_virial[n1 + 7 * N] += s_szy;
    g_virial[n1 + 8 * N] += s_szz;
    g_D_real[n1] += K_C_SP * D_real;
    g_pe[n1] += K_C_SP * s_pe;
  }
}

void find_dftd3_coordination_number(
  NEP::DFTD3& dftd3,
  const int N,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12)
{
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int n1 = 0; n1 < N; ++n1) {
    int z1 = dftd3.atomic_number[g_type[n1]];
    double R_cov_1 = dftd3para::Bohr * dftd3para::covalent_radius[z1];
    double cn_temp = 0.0;
    for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL_angular[index];
      int z2 = dftd3.atomic_number[g_type[n2]];
      double R_cov_2 = dftd3para::Bohr * dftd3para::covalent_radius[z2];
      double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      cn_temp += 1.0 / (exp(-16.0 * ((R_cov_1 + R_cov_2) / d12 - 1.0)) + 1.0);
    }
    dftd3.cn[n1] = cn_temp;
  }
}

void add_dftd3_force(
  NEP::DFTD3& dftd3,
  const int N,
  const int* g_NN_radial,
  const int* g_NL_radial,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  double* g_potential,
  double* g_force,
  double* g_virial)
{
  for (int n1 = 0; n1 < N; ++n1) {
    int z1 = dftd3.atomic_number[g_type[n1]];
    int num_cn_1 = dftd3para::num_cn[z1];
    double dc6_sum = 0.0;
    double dc8_sum = 0.0;
    for (int i1 = 0; i1 < g_NN_radial[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL_radial[index];
      int z2 = dftd3.atomic_number[g_type[n2]];
      int z_small = z1, z_large = z2;
      if (z1 > z2) {
        z_small = z2;
        z_large = z1;
      }
      int z12 = z_small * dftd3para::max_elem - (z_small * (z_small - 1)) / 2 + (z_large - z_small);
      double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      double d12_2 = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
      double d12_4 = d12_2 * d12_2;
      double d12_6 = d12_4 * d12_2;
      double d12_8 = d12_6 * d12_2;
      double c6 = 0.0;
      double dc6 = 0.0;
      int num_cn_2 = dftd3para::num_cn[z2];
      if (num_cn_1 == 1 && num_cn_2 == 1) {
        c6 = dftd3para::c6_ref[z12 * dftd3para::max_cn2];
      } else {
        double W = 0.0;
        double dW = 0.0;
        double Z = 0.0;
        double dZ = 0.0;
        for (int i = 0; i < num_cn_1; ++i) {
          for (int j = 0; j < num_cn_2; ++j) {
            double diff_i = dftd3.cn[n1] - dftd3para::cn_ref[z1 * dftd3para::max_cn + i];
            double diff_j = dftd3.cn[n2] - dftd3para::cn_ref[z2 * dftd3para::max_cn + j];
            double L_ij = exp(-4.0 * (diff_i * diff_i + diff_j * diff_j));
            W += L_ij;
            dW += L_ij * (-8.0 * diff_i);
            double c6_ref_ij =
              (z1 < z2) ? dftd3para::c6_ref[z12 * dftd3para::max_cn2 + i * dftd3para::max_cn + j]
                        : dftd3para::c6_ref[z12 * dftd3para::max_cn2 + j * dftd3para::max_cn + i];
            Z += c6_ref_ij * L_ij;
            dZ += c6_ref_ij * L_ij * (-8.0 * diff_i);
          }
        }
        if (W < 1.0e-30) {
          int i = num_cn_1 - 1;
          int j = num_cn_2 - 1;
          c6 = (z1 < z2) ? dftd3para::c6_ref[z12 * dftd3para::max_cn2 + i * dftd3para::max_cn + j]
                         : dftd3para::c6_ref[z12 * dftd3para::max_cn2 + j * dftd3para::max_cn + i];
        } else {
          W = 1.0 / W;
          c6 = Z * W;
          dc6 = dZ * W - c6 * dW * W;
        }
      }

      c6 *= dftd3para::HartreeBohr6;
      dc6 *= dftd3para::HartreeBohr6;
      double c8_over_c6 = 3.0 * dftd3para::r2r4[z1] * dftd3para::r2r4[z2] * dftd3para::Bohr2;
      double c8 = c6 * c8_over_c6;
      double damp = dftd3.a1 * sqrt(c8_over_c6) + dftd3.a2;
      double damp_2 = damp * damp;
      double damp_4 = damp_2 * damp_2;
      double damp_6 = 1.0 / (d12_6 + damp_4 * damp_2);
      double damp_8 = 1.0 / (d12_8 + damp_4 * damp_4);
      g_potential[n1] -= (dftd3.s6 * c6 * damp_6 + dftd3.s8 * c8 * damp_8) * 0.5;
      double f2 = dftd3.s6 * c6 * 3.0 * d12_4 * (damp_6 * damp_6) +
                  dftd3.s8 * c8 * 4.0 * d12_6 * (damp_8 * damp_8);
      double f12[3] = {r12[0] * f2, r12[1] * f2, r12[2] * f2};
      g_force[n1 + 0 * N] += f12[0];
      g_force[n1 + 1 * N] += f12[1];
      g_force[n1 + 2 * N] += f12[2];
      g_force[n2 + 0 * N] -= f12[0];
      g_force[n2 + 1 * N] -= f12[1];
      g_force[n2 + 2 * N] -= f12[2];
      g_virial[n2 + 0 * N] -= r12[0] * f12[0];
      g_virial[n2 + 1 * N] -= r12[0] * f12[1];
      g_virial[n2 + 2 * N] -= r12[0] * f12[2];
      g_virial[n2 + 3 * N] -= r12[1] * f12[0];
      g_virial[n2 + 4 * N] -= r12[1] * f12[1];
      g_virial[n2 + 5 * N] -= r12[1] * f12[2];
      g_virial[n2 + 6 * N] -= r12[2] * f12[0];
      g_virial[n2 + 7 * N] -= r12[2] * f12[1];
      g_virial[n2 + 8 * N] -= r12[2] * f12[2];
      dc6_sum += dc6 * dftd3.s6 * damp_6;
      dc8_sum += dc6 * c8_over_c6 * dftd3.s8 * damp_8;
    }
    dftd3.dc6_sum[n1] = dc6_sum;
    dftd3.dc8_sum[n1] = dc8_sum;
  }
}

void add_dftd3_force_extra(
  const NEP::DFTD3& dftd3,
  const int N,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* g_type,
  const double* g_x12,
  const double* g_y12,
  const double* g_z12,
  double* g_force,
  double* g_virial)
{
  for (int n1 = 0; n1 < N; ++n1) {
    int z1 = dftd3.atomic_number[g_type[n1]];
    double R_cov_1 = dftd3para::Bohr * dftd3para::covalent_radius[z1];
    double dc6_sum = dftd3.dc6_sum[n1];
    double dc8_sum = dftd3.dc8_sum[n1];
    for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL_angular[index];
      int z2 = dftd3.atomic_number[g_type[n2]];
      double R_cov_2 = dftd3para::Bohr * dftd3para::covalent_radius[z2];
      double r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      double d12_2 = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
      double d12 = sqrt(d12_2);
      double cn_exp_factor = exp(-16.0 * ((R_cov_1 + R_cov_2) / d12 - 1.0));
      double f2 = cn_exp_factor * 16.0 * (R_cov_1 + R_cov_2) * (dc6_sum + dc8_sum); // not 8.0
      f2 /= (cn_exp_factor + 1.0) * (cn_exp_factor + 1.0) * d12 * d12_2;
      double f12[3] = {r12[0] * f2, r12[1] * f2, r12[2] * f2};
      g_force[n1 + 0 * N] += f12[0];
      g_force[n1 + 1 * N] += f12[1];
      g_force[n1 + 2 * N] += f12[2];
      g_force[n2 + 0 * N] -= f12[0];
      g_force[n2 + 1 * N] -= f12[1];
      g_force[n2 + 2 * N] -= f12[2];
      g_virial[n2 + 0 * N] -= r12[0] * f12[0];
      g_virial[n2 + 1 * N] -= r12[0] * f12[1];
      g_virial[n2 + 2 * N] -= r12[0] * f12[2];
      g_virial[n2 + 3 * N] -= r12[1] * f12[0];
      g_virial[n2 + 4 * N] -= r12[1] * f12[1];
      g_virial[n2 + 5 * N] -= r12[1] * f12[2];
      g_virial[n2 + 6 * N] -= r12[2] * f12[0];
      g_virial[n2 + 7 * N] -= r12[2] * f12[1];
      g_virial[n2 + 8 * N] -= r12[2] * f12[2];
    }
  }
}

void find_descriptor_for_lammps(
  NEP::ParaMB& paramb,
  NEP::ANN& annmb,
  int nlocal,
  int N,
  int* g_ilist,
  int* g_NN,
  int** g_NL,
  int* g_type,
  int* type_map,
  double** g_pos,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gn_radial,
  const double* g_gn_angular,
#endif
  double* g_Fp,
  double* g_sum_fxyz,
  double& g_total_potential,
  double* g_potential)
{
  for (int ii = 0; ii < N; ++ii) {
    int n1 = g_ilist[ii];
    int t1 = type_map[g_type[n1]]; // from LAMMPS to NEP convention
    double q[MAX_DIM] = {0.0};

    for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
      int n2 = g_NL[n1][i1];
      int t2 = type_map[g_type[n2]]; // from LAMMPS to NEP convention
      double rc = (paramb.rc_radial[t1] + paramb.rc_radial[t2]) * 0.5;
      double rcinv = 1.0 / rc;
      double r12[3] = {
        g_pos[n2][0] - g_pos[n1][0], g_pos[n2][1] - g_pos[n1][1], g_pos[n2][2] - g_pos[n1][2]};

      double d12sq = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
      if (d12sq >= rc * rc) {
        continue;
      }
      double d12 = sqrt(d12sq);

#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
      int index_left, index_right;
      double weight_left, weight_right;
      find_index_and_weight(
        d12 * paramb.rcinv_radial, index_left, index_right, weight_left, weight_right);
      int t12 = t1 * paramb.num_types + t2;
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        q[n] +=
          g_gn_radial[(index_left * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
            weight_left +
          g_gn_radial[(index_right * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
            weight_right;
      }
#else
      double fc12;
      find_fc(rc, rcinv, d12, fc12);
      double fn12[MAX_NUM_N];
      find_fn(paramb.basis_size_radial, rcinv, d12, fc12, fn12);
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        double gn12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_radial; ++k) {
          int c_index = (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2;
          gn12 += fn12[k] * annmb.c[c_index];
        }
        q[n] += gn12;
      }
#endif
    }

    for (int n = 0; n <= paramb.n_max_angular; ++n) {
      double s[NUM_OF_ABC] = {0.0};
      for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
        int n2 = g_NL[n1][i1];
        int t2 = type_map[g_type[n2]]; // from LAMMPS to NEP convention
        double rc = (paramb.rc_angular[t1] + paramb.rc_angular[t2]) * 0.5;
        double rcinv = 1.0 / rc;

        double r12[3] = {
          g_pos[n2][0] - g_pos[n1][0], g_pos[n2][1] - g_pos[n1][1], g_pos[n2][2] - g_pos[n1][2]};

        double d12sq = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
        if (d12sq >= rc * rc) {
          continue;
        }
        double d12 = sqrt(d12sq);


#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
        int index_left, index_right;
        double weight_left, weight_right;
        find_index_and_weight(
          d12 * paramb.rcinv_angular, index_left, index_right, weight_left, weight_right);
        int t12 = t1 * paramb.num_types + t2;
        double gn12 =
          g_gn_angular[(index_left * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n] *
            weight_left +
          g_gn_angular[(index_right * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n] *
            weight_right;
        accumulate_s(paramb.L_max, d12, r12[0], r12[1], r12[2], gn12, s);
#else
        double fc12;
        find_fc(rc, rcinv, d12, fc12);
        double fn12[MAX_NUM_N];
        find_fn(paramb.basis_size_angular, rcinv, d12, fc12, fn12);
        double gn12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_angular; ++k) {
          int c_index = (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2 + paramb.num_c_radial;
          gn12 += fn12[k] * annmb.c[c_index];
        }
        accumulate_s(paramb.L_max, d12, r12[0], r12[1], r12[2], gn12, s);
#endif
      }
      find_q(
        paramb.L_max, paramb.num_L, paramb.n_max_angular + 1, n, s, q + (paramb.n_max_radial + 1));
      for (int abc = 0; abc < NUM_OF_ABC; ++abc) {
        g_sum_fxyz[(n * NUM_OF_ABC + abc) * nlocal + n1] = s[abc];
      }
    }

    for (int d = 0; d < annmb.dim; ++d) {
      q[d] = q[d] * paramb.q_scaler[d];
    }

    double F = 0.0, Fp[MAX_DIM] = {0.0}, latent_space[MAX_NEURON] = {0.0};

    if (paramb.version == 5) {
      apply_ann_one_layer_nep5(
        annmb.dim, annmb.num_neurons1, annmb.w0[t1], annmb.b0[t1], annmb.w1[t1], annmb.b1, q, F, Fp,
        latent_space);
    } else {
      apply_ann_one_layer(
        annmb.dim, annmb.num_neurons1, annmb.w0[t1], annmb.b0[t1], annmb.w1[t1], annmb.b1, q, F, Fp,
        latent_space, false, nullptr);
    }

    g_total_potential += F; // always calculate this
    if (g_potential) {      // only calculate when required
      g_potential[n1] += F;
    }

    for (int d = 0; d < annmb.dim; ++d) {
      g_Fp[d * nlocal + n1] = Fp[d] * paramb.q_scaler[d];
    }
  }
}

void find_force_radial_for_lammps(
  NEP::ParaMB& paramb,
  NEP::ANN& annmb,
  int nlocal,
  int N,
  int* g_ilist,
  int* g_NN,
  int** g_NL,
  int* g_type,
  int* type_map,
  double** g_pos,
  double* g_Fp,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gnp_radial,
#endif
  double** g_force,
  double g_total_virial[6],
  double** g_virial)
{
  for (int ii = 0; ii < N; ++ii) {
    int n1 = g_ilist[ii];
    int t1 = type_map[g_type[n1]]; // from LAMMPS to NEP convention
    for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
      int n2 = g_NL[n1][i1];
      int t2 = type_map[g_type[n2]]; // from LAMMPS to NEP convention
      double rc = (paramb.rc_radial[t1] + paramb.rc_radial[t2]) * 0.5;
      double rcinv = 1.0 / rc;
      double r12[3] = {
        g_pos[n2][0] - g_pos[n1][0], g_pos[n2][1] - g_pos[n1][1], g_pos[n2][2] - g_pos[n1][2]};

      double d12sq = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
      if (d12sq >= rc * rc) {
        continue;
      }
      double d12 = sqrt(d12sq);
      double d12inv = 1.0 / d12;
      double f12[3] = {0.0};
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
      int index_left, index_right;
      double weight_left, weight_right;
      find_index_and_weight(
        d12 * paramb.rcinv_radial, index_left, index_right, weight_left, weight_right);
      int t12 = t1 * paramb.num_types + t2;
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        double gnp12 =
          g_gnp_radial[(index_left * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
            weight_left +
          g_gnp_radial[(index_right * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
            weight_right;
        double tmp12 = g_Fp[n1 + n * nlocal] * gnp12 * d12inv;
        for (int d = 0; d < 3; ++d) {
          f12[d] += tmp12 * r12[d];
        }
      }
#else
      double fc12, fcp12;
      find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);
      double fn12[MAX_NUM_N];
      double fnp12[MAX_NUM_N];
      find_fn_and_fnp(paramb.basis_size_radial, rcinv, d12, fc12, fcp12, fn12, fnp12);
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        double gnp12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_radial; ++k) {
          int c_index = (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2;
          gnp12 += fnp12[k] * annmb.c[c_index];
        }
        double tmp12 = g_Fp[n1 + n * nlocal] * gnp12 * d12inv;
        for (int d = 0; d < 3; ++d) {
          f12[d] += tmp12 * r12[d];
        }
      }
#endif

      g_force[n1][0] += f12[0];
      g_force[n1][1] += f12[1];
      g_force[n1][2] += f12[2];
      g_force[n2][0] -= f12[0];
      g_force[n2][1] -= f12[1];
      g_force[n2][2] -= f12[2];

      // always calculate the total virial:
      g_total_virial[0] -= r12[0] * f12[0]; // xx
      g_total_virial[1] -= r12[1] * f12[1]; // yy
      g_total_virial[2] -= r12[2] * f12[2]; // zz
      g_total_virial[3] -= r12[0] * f12[1]; // xy
      g_total_virial[4] -= r12[0] * f12[2]; // xz
      g_total_virial[5] -= r12[1] * f12[2]; // yz
      if (g_virial) {                       // only calculate the per-atom virial when required
        g_virial[n2][0] -= r12[0] * f12[0]; // xx
        g_virial[n2][1] -= r12[1] * f12[1]; // yy
        g_virial[n2][2] -= r12[2] * f12[2]; // zz
        g_virial[n2][3] -= r12[0] * f12[1]; // xy
        g_virial[n2][4] -= r12[0] * f12[2]; // xz
        g_virial[n2][5] -= r12[1] * f12[2]; // yz
        g_virial[n2][6] -= r12[1] * f12[0]; // yx
        g_virial[n2][7] -= r12[2] * f12[0]; // zx
        g_virial[n2][8] -= r12[2] * f12[1]; // zy
      }
    }
  }
}

void find_force_angular_for_lammps(
  NEP::ParaMB& paramb,
  NEP::ANN& annmb,
  int nlocal,
  int N,
  int* g_ilist,
  int* g_NN,
  int** g_NL,
  int* g_type,
  int* type_map,
  double** g_pos,
  double* g_Fp,
  double* g_sum_fxyz,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  const double* g_gn_angular,
  const double* g_gnp_angular,
#endif
  double** g_force,
  double g_total_virial[6],
  double** g_virial)
{
  for (int ii = 0; ii < N; ++ii) {
    int n1 = g_ilist[ii];
    double Fp[MAX_DIM_ANGULAR] = {0.0};
    double sum_fxyz[NUM_OF_ABC * MAX_NUM_N];
    for (int d = 0; d < paramb.dim_angular; ++d) {
      Fp[d] = g_Fp[(paramb.n_max_radial + 1 + d) * nlocal + n1];
    }
    for (int d = 0; d < (paramb.n_max_angular + 1) * NUM_OF_ABC; ++d) {
      sum_fxyz[d] = g_sum_fxyz[d * nlocal + n1];
    }

    int t1 = type_map[g_type[n1]]; // from LAMMPS to NEP convention

    for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
      int n2 = g_NL[n1][i1];
      int t2 = type_map[g_type[n2]]; // from LAMMPS to NEP convention
      double rc = (paramb.rc_angular[t1] + paramb.rc_angular[t2]) * 0.5;
      double rcinv = 1.0 / rc;
      double r12[3] = {
        g_pos[n2][0] - g_pos[n1][0], g_pos[n2][1] - g_pos[n1][1], g_pos[n2][2] - g_pos[n1][2]};

      double d12sq = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
      if (d12sq >= rc * rc) {
        continue;
      }
      double d12 = sqrt(d12sq);
      double f12[3] = {0.0};

#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
      int index_left, index_right;
      double weight_left, weight_right;
      find_index_and_weight(
        d12 * paramb.rcinv_angular, index_left, index_right, weight_left, weight_right);
      int t12 = t1 * paramb.num_types + t2;
      for (int n = 0; n <= paramb.n_max_angular; ++n) {
        int index_left_all =
          (index_left * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n;
        int index_right_all =
          (index_right * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n;
        double gn12 =
          g_gn_angular[index_left_all] * weight_left + g_gn_angular[index_right_all] * weight_right;
        double gnp12 = g_gnp_angular[index_left_all] * weight_left +
                       g_gnp_angular[index_right_all] * weight_right;
        accumulate_f12(
          paramb.L_max, paramb.num_L, n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp,
          sum_fxyz, f12);
      }
#else
      double fc12, fcp12;
      find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);
      double fn12[MAX_NUM_N];
      double fnp12[MAX_NUM_N];
      find_fn_and_fnp(paramb.basis_size_angular, rcinv, d12, fc12, fcp12, fn12, fnp12);
      for (int n = 0; n <= paramb.n_max_angular; ++n) {
        double gn12 = 0.0;
        double gnp12 = 0.0;
        for (int k = 0; k <= paramb.basis_size_angular; ++k) {
          int c_index = (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq;
          c_index += t1 * paramb.num_types + t2 + paramb.num_c_radial;
          gn12 += fn12[k] * annmb.c[c_index];
          gnp12 += fnp12[k] * annmb.c[c_index];
        }
        accumulate_f12(
          paramb.L_max, paramb.num_L, n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp,
          sum_fxyz, f12);
      }
#endif

      g_force[n1][0] += f12[0];
      g_force[n1][1] += f12[1];
      g_force[n1][2] += f12[2];
      g_force[n2][0] -= f12[0];
      g_force[n2][1] -= f12[1];
      g_force[n2][2] -= f12[2];
      // always calculate the total virial:
      g_total_virial[0] -= r12[0] * f12[0]; // xx
      g_total_virial[1] -= r12[1] * f12[1]; // yy
      g_total_virial[2] -= r12[2] * f12[2]; // zz
      g_total_virial[3] -= r12[0] * f12[1]; // xy
      g_total_virial[4] -= r12[0] * f12[2]; // xz
      g_total_virial[5] -= r12[1] * f12[2]; // yz
      if (g_virial) {                       // only calculate the per-atom virial when required
        g_virial[n2][0] -= r12[0] * f12[0]; // xx
        g_virial[n2][1] -= r12[1] * f12[1]; // yy
        g_virial[n2][2] -= r12[2] * f12[2]; // zz
        g_virial[n2][3] -= r12[0] * f12[1]; // xy
        g_virial[n2][4] -= r12[0] * f12[2]; // xz
        g_virial[n2][5] -= r12[1] * f12[2]; // yz
        g_virial[n2][6] -= r12[1] * f12[0]; // yx
        g_virial[n2][7] -= r12[2] * f12[0]; // zx
        g_virial[n2][8] -= r12[2] * f12[1]; // zy
      }
    }
  }
}

void find_force_ZBL_for_lammps(
  NEP::ParaMB& paramb,
  const NEP::ZBL& zbl,
  int N,
  int* g_ilist,
  int* g_NN,
  int** g_NL,
  int* g_type,
  int* type_map,
  double** g_pos,
  double** g_force,
  double g_total_virial[6],
  double** g_virial,
  double& g_total_potential,
  double* g_potential)
{
  for (int ii = 0; ii < N; ++ii) {
    int n1 = g_ilist[ii];
    int type1 = type_map[g_type[n1]]; // from LAMMPS to NEP convention
    int zi = paramb.atomic_numbers[type1] + 1;
    double pow_zi = pow(double(zi), 0.23);
    for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
      int n2 = g_NL[n1][i1];
      double r12[3] = {
        g_pos[n2][0] - g_pos[n1][0], g_pos[n2][1] - g_pos[n1][1], g_pos[n2][2] - g_pos[n1][2]};

      double d12sq = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
      double max_rc_outer = 2.5;
      if (d12sq >= max_rc_outer * max_rc_outer) {
        continue;
      }
      double d12 = sqrt(d12sq);

      double d12inv = 1.0 / d12;
      double f, fp;
      int type2 = type_map[g_type[n2]]; // from LAMMPS to NEP convention
      int zj = paramb.atomic_numbers[type2] + 1;
      double a_inv = (pow_zi + pow(double(zj), 0.23)) * 2.134563;
      double zizj = K_C_SP * zi * zj;
      if (zbl.flexibled) {
        int t1, t2;
        if (type1 < type2) {
          t1 = type1;
          t2 = type2;
        } else {
          t1 = type2;
          t2 = type1;
        }
        int zbl_index = t1 * zbl.num_types - (t1 * (t1 - 1)) / 2 + (t2 - t1);
        double ZBL_para[10];
        for (int i = 0; i < 10; ++i) {
          ZBL_para[i] = zbl.para[10 * zbl_index + i];
        }
        find_f_and_fp_zbl(ZBL_para, zizj, a_inv, d12, d12inv, f, fp);
      } else {
        double rc_inner = zbl.rc_inner;
        double rc_outer = zbl.rc_outer;
        if (paramb.use_typewise_cutoff_zbl) {
          // zi and zj start from 1, so need to minus 1 here
          rc_outer = std::min(
            (COVALENT_RADIUS[zi - 1] + COVALENT_RADIUS[zj - 1]) * paramb.typewise_cutoff_zbl_factor,
            rc_outer);
          rc_inner = 0.0;
        }
        find_f_and_fp_zbl(zizj, a_inv, rc_inner, rc_outer, d12, d12inv, f, fp);
      }
      double f2 = fp * d12inv * 0.5;
      double f12[3] = {r12[0] * f2, r12[1] * f2, r12[2] * f2};
      g_force[n1][0] += f12[0]; // accumulation here
      g_force[n1][1] += f12[1];
      g_force[n1][2] += f12[2];
      g_force[n2][0] -= f12[0];
      g_force[n2][1] -= f12[1];
      g_force[n2][2] -= f12[2];
      // always calculate the total virial:
      g_total_virial[0] -= r12[0] * f12[0]; // xx
      g_total_virial[1] -= r12[1] * f12[1]; // yy
      g_total_virial[2] -= r12[2] * f12[2]; // zz
      g_total_virial[3] -= r12[0] * f12[1]; // xy
      g_total_virial[4] -= r12[0] * f12[2]; // xz
      g_total_virial[5] -= r12[1] * f12[2]; // yz
      if (g_virial) {                       // only calculate the per-atom virial when required
        g_virial[n2][0] -= r12[0] * f12[0]; // xx
        g_virial[n2][1] -= r12[1] * f12[1]; // yy
        g_virial[n2][2] -= r12[2] * f12[2]; // zz
        g_virial[n2][3] -= r12[0] * f12[1]; // xy
        g_virial[n2][4] -= r12[0] * f12[2]; // xz
        g_virial[n2][5] -= r12[1] * f12[2]; // yz
        g_virial[n2][6] -= r12[1] * f12[0]; // yx
        g_virial[n2][7] -= r12[2] * f12[0]; // zx
        g_virial[n2][8] -= r12[2] * f12[1]; // zy
      }
      g_total_potential += f * 0.5; // always calculate this
      if (g_potential) {            // only calculate when required
        g_potential[n1] += f * 0.5;
      }
    }
  }
}

std::vector<std::string> get_tokens(std::ifstream& input)
{
  std::string line;
  std::getline(input, line);
  std::istringstream iss(line);
  std::vector<std::string> tokens{
    std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
  return tokens;
}

void print_tokens(const std::vector<std::string>& tokens)
{
  std::cout << "Line:";
  for (const auto& token : tokens) {
    std::cout << " " << token;
  }
  std::cout << std::endl;
}

int get_int_from_token(const std::string& token, const char* filename, const int line)
{
  int value = 0;
  try {
    value = std::stoi(token);
  } catch (const std::exception& e) {
    std::cout << "Standard exception:\n";
    std::cout << "    File:          " << filename << std::endl;
    std::cout << "    Line:          " << line << std::endl;
    std::cout << "    Error message: " << e.what() << std::endl;
    exit(1);
  }
  return value;
}

double get_double_from_token(const std::string& token, const char* filename, const int line)
{
  double value = 0;
  try {
    value = std::stod(token);
  } catch (const std::exception& e) {
    std::cout << "Standard exception:\n";
    std::cout << "    File:          " << filename << std::endl;
    std::cout << "    Line:          " << line << std::endl;
    std::cout << "    Error message: " << e.what() << std::endl;
    exit(1);
  }
  return value;
}

} // namespace

NEP::NEP() {}

NEP::NEP(const std::string& potential_filename) { init_from_file(potential_filename, true); }

void NEP::init_from_file(const std::string& potential_filename, const bool is_rank_0)
{
  std::ifstream input(potential_filename);
  if (!input.is_open()) {
    std::cout << "Failed to open " << potential_filename << std::endl;
    exit(1);
  }

  std::vector<std::string> tokens = get_tokens(input);
  if (tokens.size() < 3) {
    print_tokens(tokens);
    std::cout << "The first line of nep.txt should have at least 3 items." << std::endl;
    exit(1);
  }
  if (tokens[0] == "nep3") {
    paramb.model_type = 0;
    paramb.version = 3;
    zbl.enabled = false;
  } else if (tokens[0] == "nep3_zbl") {
    paramb.model_type = 0;
    paramb.version = 3;
    zbl.enabled = true;
  } else if (tokens[0] == "nep3_dipole") {
    paramb.model_type = 1;
    paramb.version = 3;
    zbl.enabled = false;
  } else if (tokens[0] == "nep3_polarizability") {
    paramb.model_type = 2;
    paramb.version = 3;
    zbl.enabled = false;
  } else if (tokens[0] == "nep4") {
    paramb.model_type = 0;
    paramb.version = 4;
    zbl.enabled = false;
  } else if (tokens[0] == "nep4_zbl") {
    paramb.model_type = 0;
    paramb.version = 4;
    zbl.enabled = true;
  } else if (tokens[0] == "nep4_dipole") {
    paramb.model_type = 1;
    paramb.version = 4;
    zbl.enabled = false;
  } else if (tokens[0] == "nep4_polarizability") {
    paramb.model_type = 2;
    paramb.version = 4;
    zbl.enabled = false;
  } else if (tokens[0] == "nep5") {
    paramb.model_type = 0;
    paramb.version = 5;
    zbl.enabled = false;
  } else if (tokens[0] == "nep5_zbl") {
    paramb.model_type = 0;
    paramb.version = 5;
    zbl.enabled = true;
  } else if (tokens[0] == "nep4_charge1") {
    paramb.model_type = 0;
    paramb.version = 4;
    zbl.enabled = false;
    paramb.charge_mode = 1;
  } else if (tokens[0] == "nep4_zbl_charge1") {
    paramb.model_type = 0;
    paramb.version = 4;
    zbl.enabled = true;
    paramb.charge_mode = 1;
  } else if (tokens[0] == "nep4_charge2") {
    paramb.model_type = 0;
    paramb.version = 4;
    zbl.enabled = false;
    paramb.charge_mode = 2;
  } else if (tokens[0] == "nep4_zbl_charge2") {
    paramb.model_type = 0;
    paramb.version = 4;
    zbl.enabled = true;
    paramb.charge_mode = 2;
  } else if (tokens[0] == "nep4_charge3") {
    paramb.model_type = 0;
    paramb.version = 4;
    zbl.enabled = false;
    paramb.charge_mode = 3;
  } else if (tokens[0] == "nep4_zbl_charge3") {
    paramb.model_type = 0;
    paramb.version = 4;
    zbl.enabled = true;
    paramb.charge_mode = 3;
  } else {
    std::cout << tokens[0] << " is an unsupported NEP model." << std::endl;
    exit(1);
  }

  paramb.num_types = get_int_from_token(tokens[1], __FILE__, __LINE__);
  if (tokens.size() != 2 + paramb.num_types) {
    print_tokens(tokens);
    std::cout << "The first line of nep.txt should have " << paramb.num_types << " atom symbols."
              << std::endl;
    exit(1);
  }

  element_list.resize(paramb.num_types);
  for (std::size_t n = 0; n < paramb.num_types; ++n) {
    int atomic_number = 0;
    element_list[n] = tokens[2 + n];
    for (int m = 0; m < NUM_ELEMENTS; ++m) {
      if (tokens[2 + n] == ELEMENTS[m]) {
        atomic_number = m;
        break;
      }
    }
    paramb.atomic_numbers[n] = atomic_number;
    dftd3.atomic_number[n] = atomic_number;
  }

  // zbl
  if (zbl.enabled) {
    tokens = get_tokens(input);
    if (tokens.size() != 3 && tokens.size() != 4) {
      print_tokens(tokens);
      std::cout << "This line should be zbl rc_inner rc_outer [zbl_factor]." << std::endl;
      exit(1);
    }
    zbl.rc_inner = get_double_from_token(tokens[1], __FILE__, __LINE__);
    zbl.rc_outer = get_double_from_token(tokens[2], __FILE__, __LINE__);
    if (zbl.rc_inner == 0 && zbl.rc_outer == 0) {
      zbl.flexibled = true;
    } else {
      if (tokens.size() == 4) {
        paramb.typewise_cutoff_zbl_factor = get_double_from_token(tokens[3], __FILE__, __LINE__);
        paramb.use_typewise_cutoff_zbl = true;
      }
    }
  }

  // cutoff
  tokens = get_tokens(input);
  if (tokens.size() != 5 && tokens.size() != paramb.num_types * 2 + 3) {
    print_tokens(tokens);
    std::cout << "cutoff should have 4 or num_types * 2 + 2 parameters.\n";
    exit(1);
  }
  if (tokens.size() == 5) {
    paramb.rc_radial[0] = get_double_from_token(tokens[1], __FILE__, __LINE__);
    paramb.rc_angular[0] = get_double_from_token(tokens[2], __FILE__, __LINE__);
    for (std::size_t n = 0; n < paramb.num_types; ++n) {
      paramb.rc_radial[n] = paramb.rc_radial[0];
      paramb.rc_angular[n] = paramb.rc_angular[0];
    }
  } else {
    for (std::size_t n = 0; n < paramb.num_types; ++n) {
      paramb.rc_radial[n] = get_double_from_token(tokens[1 + n * 2], __FILE__, __LINE__);
      paramb.rc_angular[n] = get_double_from_token(tokens[2 + n * 2], __FILE__, __LINE__);
    }
  }
  for (std::size_t n = 0; n < paramb.num_types; ++n) {
    if (paramb.rc_radial[n] > paramb.rc_radial_max) {
      paramb.rc_radial_max = paramb.rc_radial[n];
    }
    if (paramb.rc_angular[n] > paramb.rc_angular_max) {
      paramb.rc_angular_max = paramb.rc_angular[n];
    }
  }

  int MN_radial = get_int_from_token(tokens[tokens.size() - 2], __FILE__, __LINE__);
  int MN_angular = get_int_from_token(tokens[tokens.size() - 1], __FILE__, __LINE__);

  // n_max 10 8
  tokens = get_tokens(input);
  if (tokens.size() != 3) {
    print_tokens(tokens);
    std::cout << "This line should be n_max n_max_radial n_max_angular." << std::endl;
    exit(1);
  }
  paramb.n_max_radial = get_int_from_token(tokens[1], __FILE__, __LINE__);
  paramb.n_max_angular = get_int_from_token(tokens[2], __FILE__, __LINE__);

  // basis_size 10 8
  tokens = get_tokens(input);
  if (tokens.size() != 3) {
    print_tokens(tokens);
    std::cout << "This line should be basis_size basis_size_radial basis_size_angular."
              << std::endl;
    exit(1);
  }
  paramb.basis_size_radial = get_int_from_token(tokens[1], __FILE__, __LINE__);
  paramb.basis_size_angular = get_int_from_token(tokens[2], __FILE__, __LINE__);

  // l_max
  tokens = get_tokens(input);
  if (tokens.size() != 4) {
    print_tokens(tokens);
    std::cout << "This line should be l_max l_max_3body l_max_4body l_max_5body." << std::endl;
    exit(1);
  }

  paramb.L_max = get_int_from_token(tokens[1], __FILE__, __LINE__);
  paramb.num_L = paramb.L_max;

  int L_max_4body = get_int_from_token(tokens[2], __FILE__, __LINE__);
  int L_max_5body = get_int_from_token(tokens[3], __FILE__, __LINE__);
  if (L_max_4body == 2) {
    paramb.num_L += 1;
  }
  if (L_max_5body == 1) {
    paramb.num_L += 1;
  }

  paramb.dim_angular = (paramb.n_max_angular + 1) * paramb.num_L;

  // ANN
  tokens = get_tokens(input);
  if (tokens.size() != 3) {
    print_tokens(tokens);
    std::cout << "This line should be ANN num_neurons 0." << std::endl;
    exit(1);
  }
  annmb.num_neurons1 = get_int_from_token(tokens[1], __FILE__, __LINE__);
  annmb.dim = (paramb.n_max_radial + 1) + paramb.dim_angular;

  // calculated parameters:
  paramb.num_types_sq = paramb.num_types * paramb.num_types;
  if (paramb.version == 3) {
    annmb.num_para_ann = (annmb.dim + 2) * annmb.num_neurons1 + 1;
  } else if (paramb.version == 4) {
    annmb.num_para_ann = (annmb.dim + 2) * annmb.num_neurons1 * paramb.num_types + 1;
  } else {
    annmb.num_para_ann = ((annmb.dim + 2) * annmb.num_neurons1 + 1) * paramb.num_types + 1;
  }
  if (paramb.model_type == 2) {
    annmb.num_para_ann *= 2;
  }
  if (paramb.charge_mode > 0) {
    annmb.num_para_ann += annmb.num_neurons1 * paramb.num_types + 1;
  }
  int num_para_descriptor =
    paramb.num_types_sq * ((paramb.n_max_radial + 1) * (paramb.basis_size_radial + 1) +
                           (paramb.n_max_angular + 1) * (paramb.basis_size_angular + 1));
  annmb.num_para = annmb.num_para_ann + num_para_descriptor;

  paramb.num_c_radial =
    paramb.num_types_sq * (paramb.n_max_radial + 1) * (paramb.basis_size_radial + 1);

  // NN and descriptor parameters
  parameters.resize(annmb.num_para);
  for (int n = 0; n < annmb.num_para; ++n) {
    tokens = get_tokens(input);
    parameters[n] = get_double_from_token(tokens[0], __FILE__, __LINE__);
  }
  update_potential(parameters.data(), annmb);
  for (int d = 0; d < annmb.dim; ++d) {
    tokens = get_tokens(input);
    paramb.q_scaler[d] = get_double_from_token(tokens[0], __FILE__, __LINE__);
  }

  // flexible zbl potential parameters if (zbl.flexibled)
  if (zbl.flexibled) {
    int num_type_zbl = (paramb.num_types * (paramb.num_types + 1)) / 2;
    for (int d = 0; d < 10 * num_type_zbl; ++d) {
      tokens = get_tokens(input);
      zbl.para[d] = get_double_from_token(tokens[0], __FILE__, __LINE__);
    }
    zbl.num_types = paramb.num_types;
  }
  input.close();


  // charge related parameters and data
  if (paramb.charge_mode > 0) {
    charge_para.alpha = PI / paramb.rc_radial_max; // a good value
    ewald.initialize(charge_para.alpha);
    charge_para.two_alpha_over_sqrt_pi = 2.0 * charge_para.alpha / sqrt(PI);
    charge_para.A = erfc(PI) / (paramb.rc_radial_max * paramb.rc_radial_max);
    charge_para.A += charge_para.two_alpha_over_sqrt_pi * exp(-PI * PI) / paramb.rc_radial_max;
    charge_para.B = - erfc(PI) / paramb.rc_radial_max - charge_para.A * paramb.rc_radial_max;
  }

#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
  if (paramb.use_typewise_cutoff) {
    std::cout << "Cannot use tabulated radial functions with typewise cutoff." << std::endl;
    exit(1);
  }
  construct_table(parameters.data());
#endif

  // only report for rank_0
  if (is_rank_0) {

    if (paramb.charge_mode > 0) {
      if (paramb.num_types == 1) {
        std::cout << "Use the NEP4-Charge" << paramb.charge_mode << " potential with " << paramb.num_types
                  << " atom type.\n";
      } else {
        std::cout << "Use the NEP4-Charge" << paramb.charge_mode << " potential with " << paramb.num_types
                  << " atom types.\n";
      }
    } else {
      if (paramb.num_types == 1) {
        std::cout << "Use the NEP" << paramb.version << " potential with " << paramb.num_types
                  << " atom type.\n";
      } else {
        std::cout << "Use the NEP" << paramb.version << " potential with " << paramb.num_types
                  << " atom types.\n";
      }
    }

    for (std::size_t n = 0; n < paramb.num_types; ++n) {
      std::cout << "    type " << n << " (" << element_list[n]
                << " with Z = " << paramb.atomic_numbers[n] + 1 << ")"
                << " has cutoffs " << "(" << paramb.rc_radial[n] << " A, "
                << paramb.rc_angular[n] << " A).\n";
    }

    if (zbl.enabled) {
      if (zbl.flexibled) {
        std::cout << "    has flexible ZBL.\n";
      } else {
        std::cout << "    has universal ZBL with inner cutoff " << zbl.rc_inner
                  << " A and outer cutoff " << zbl.rc_outer << " A.\n";
        if (paramb.use_typewise_cutoff_zbl) {
          std::cout << "    ZBL typewise cutoff is enabled with factor "
                    << paramb.typewise_cutoff_zbl_factor << ".\n";
        }
      }
    }

    std::cout << "    n_max_radial = " << paramb.n_max_radial << ".\n";
    std::cout << "    n_max_angular = " << paramb.n_max_angular << ".\n";
    std::cout << "    basis_size_radial = " << paramb.basis_size_radial << ".\n";
    std::cout << "    basis_size_angular = " << paramb.basis_size_angular << ".\n";
    std::cout << "    l_max_3body = " << paramb.L_max << ".\n";
    std::cout << "    l_max_4body = " << (paramb.num_L >= 5 ? 2 : 0) << ".\n";
    std::cout << "    l_max_5body = " << (paramb.num_L >= 6 ? 1 : 0) << ".\n";
    std::cout << "    ANN = " << annmb.dim << "-" << annmb.num_neurons1 << "-1.\n";
    std::cout << "    number of neural network parameters = " << annmb.num_para_ann << ".\n";
    std::cout << "    number of descriptor parameters = " << num_para_descriptor << ".\n";
    std::cout << "    total number of parameters = " << annmb.num_para << ".\n";
  }
}

void NEP::update_type_map(const int ntype, int* type_map, char** elements)
{
  std::size_t n = 0;
  for (int itype = 0; itype < ntype + 1; ++itype) {
    // check if set NULL in lammps input file
    if (type_map[itype] == -1) {
      continue;
    }

    // find the same element name in potential file
    std::string element_name = elements[type_map[itype]];
    for (n = 0; n < paramb.num_types; ++n) {
      if (element_name == element_list[n]) {
        type_map[itype] = n;
        break;
      }
    }

    // check if no corresponding element
    if (n == paramb.num_types) {
      std::cout << "There is no element " << element_name << " in the potential file." << std::endl;
      exit(1);
    }
  }
}

void NEP::update_potential(double* parameters, ANN& ann)
{
  double* pointer = parameters;
  for (std::size_t t = 0; t < paramb.num_types; ++t) {
    if (t > 0 && paramb.version == 3) { // Use the same set of NN parameters for NEP3
      pointer -= (ann.dim + 2) * ann.num_neurons1;
    }
    ann.w0[t] = pointer;
    pointer += ann.num_neurons1 * ann.dim;
    ann.b0[t] = pointer;
    pointer += ann.num_neurons1;
    ann.w1[t] = pointer;
    if (paramb.charge_mode > 0) {
      pointer += ann.num_neurons1 * 2;
    } else {
      pointer += ann.num_neurons1;
    }
    
    if (paramb.version == 5) {
      pointer += 1; // one extra bias for NEP5 stored in ann.w1[t]
    }
  }

  if (paramb.charge_mode > 0) {
    ann.sqrt_epsilon_inf = pointer;
    pointer += 1;
  }

  ann.b1 = pointer;
  pointer += 1;

  if (paramb.model_type == 2) {
    for (std::size_t t = 0; t < paramb.num_types; ++t) {
      if (t > 0 && paramb.version == 3) { // Use the same set of NN parameters for NEP3
        pointer -= (ann.dim + 2) * ann.num_neurons1;
      }
      ann.w0_pol[t] = pointer;
      pointer += ann.num_neurons1 * ann.dim;
      ann.b0_pol[t] = pointer;
      pointer += ann.num_neurons1;
      ann.w1_pol[t] = pointer;
      pointer += ann.num_neurons1;
    }
    ann.b1_pol = pointer;
    pointer += 1;
  }

  ann.c = pointer;
}

#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
void NEP::construct_table(double* parameters)
{
  gn_radial.resize(table_length * paramb.num_types_sq * (paramb.n_max_radial + 1));
  gnp_radial.resize(table_length * paramb.num_types_sq * (paramb.n_max_radial + 1));
  gn_angular.resize(table_length * paramb.num_types_sq * (paramb.n_max_angular + 1));
  gnp_angular.resize(table_length * paramb.num_types_sq * (paramb.n_max_angular + 1));
  double* c_pointer = parameters + annmb.num_para_ann;
  construct_table_radial_or_angular(
    paramb.version, paramb.num_types, paramb.num_types_sq, paramb.n_max_radial,
    paramb.basis_size_radial, paramb.rc_radial, paramb.rcinv_radial, c_pointer, gn_radial.data(),
    gnp_radial.data());
  construct_table_radial_or_angular(
    paramb.version, paramb.num_types, paramb.num_types_sq, paramb.n_max_angular,
    paramb.basis_size_angular, paramb.rc_angular, paramb.rcinv_angular,
    c_pointer + paramb.num_c_radial, gn_angular.data(), gnp_angular.data());
}
#endif

void NEP::allocate_memory(const int N)
{
  if (num_atoms < N) {
    NN_radial.resize(N);
    NL_radial.resize(N * MN);
    NN_angular.resize(N);
    NL_angular.resize(N * MN);
    r12.resize(N * MN * 6);
    Fp.resize(N * annmb.dim);
    sum_fxyz.resize(N * (paramb.n_max_angular + 1) * NUM_OF_ABC);
    if (paramb.charge_mode > 0) {
      D_real.resize(N);
      charge_derivative.resize(N * annmb.dim);
    }
    dftd3.cn.resize(N);
    dftd3.dc6_sum.resize(N);
    dftd3.dc8_sum.resize(N);
    num_atoms = N;
  }
}

void NEP::compute(
  const std::vector<int>& type,
  const std::vector<double>& box,
  const std::vector<double>& position,
  std::vector<double>& potential,
  std::vector<double>& force,
  std::vector<double>& virial)
{
  if (paramb.model_type != 0) {
    std::cout << "Cannot compute potential using a non-potential NEP model.\n";
    exit(1);
  }

  if (paramb.charge_mode != 0) {
    std::cout << "Cannot use this compute for a qNEP model.\n";
    exit(1);
  }

  const std::size_t N = type.size();
  const int size_x12 = N * MN;

  if (N * 3 != position.size()) {
    std::cout << "Type and position sizes are inconsistent.\n";
    exit(1);
  }
  if (N != potential.size()) {
    std::cout << "Type and potential sizes are inconsistent.\n";
    exit(1);
  }
  if (N * 3 != force.size()) {
    std::cout << "Type and force sizes are inconsistent.\n";
    exit(1);
  }
  if (N * 9 != virial.size()) {
    std::cout << "Type and virial sizes are inconsistent.\n";
    exit(1);
  }

  allocate_memory(N);

  for (std::size_t n = 0; n < potential.size(); ++n) {
    potential[n] = 0.0;
  }
  for (std::size_t n = 0; n < force.size(); ++n) {
    force[n] = 0.0;
  }
  for (std::size_t n = 0; n < virial.size(); ++n) {
    virial[n] = 0.0;
  }

  find_neighbor_list_small_box(
    paramb.rc_radial_max, paramb.rc_angular_max, N, MN, box, position, num_cells, ebox, NN_radial, NL_radial,
    NN_angular, NL_angular, r12);

  find_descriptor_small_box(
    true, false, false, false, paramb, annmb, N, NN_radial.data(), NL_radial.data(),
    NN_angular.data(), NL_angular.data(), type.data(), r12.data(), r12.data() + size_x12,
    r12.data() + size_x12 * 2, r12.data() + size_x12 * 3, r12.data() + size_x12 * 4,
    r12.data() + size_x12 * 5,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gn_radial.data(), gn_angular.data(),
#endif
    Fp.data(), sum_fxyz.data(), potential.data(), nullptr, nullptr, nullptr, false, nullptr);

  find_force_radial_small_box(
    false, paramb, annmb, N, NN_radial.data(), NL_radial.data(), type.data(), r12.data(),
    r12.data() + size_x12, r12.data() + size_x12 * 2, Fp.data(),
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gnp_radial.data(),
#endif
    force.data(), force.data() + N, force.data() + N * 2, virial.data());

  find_force_angular_small_box(
    false, paramb, annmb, N, NN_angular.data(), NL_angular.data(), type.data(),
    r12.data() + size_x12 * 3, r12.data() + size_x12 * 4, r12.data() + size_x12 * 5, Fp.data(),
    sum_fxyz.data(),
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gn_angular.data(), gnp_angular.data(),
#endif
    force.data(), force.data() + N, force.data() + N * 2, virial.data());

  if (zbl.enabled) {
    find_force_ZBL_small_box(
      N, paramb, zbl, NN_angular.data(), NL_angular.data(), type.data(), r12.data() + size_x12 * 3,
      r12.data() + size_x12 * 4, r12.data() + size_x12 * 5, force.data(), force.data() + N,
      force.data() + N * 2, virial.data(), potential.data());
  }
}

void NEP::compute(
  const std::vector<int>& type,
  const std::vector<double>& box,
  const std::vector<double>& position,
  std::vector<double>& potential,
  std::vector<double>& force,
  std::vector<double>& virial,
  std::vector<double>& charge,
  std::vector<double>& bec)
{
  if (paramb.charge_mode == 0) {
    std::cout << "Can only use this compute for a qNEP model.\n";
    exit(1);
  }

  const std::size_t N = type.size();
  const std::size_t size_x12 = N * MN;

  if (N * 3 != position.size()) {
    std::cout << "Type and position sizes are inconsistent.\n";
    exit(1);
  }
  if (N != potential.size()) {
    std::cout << "Type and potential sizes are inconsistent.\n";
    exit(1);
  }
  if (N * 3 != force.size()) {
    std::cout << "Type and force sizes are inconsistent.\n";
    exit(1);
  }
  if (N * 9 != virial.size()) {
    std::cout << "Type and virial sizes are inconsistent.\n";
    exit(1);
  }
  if (N != charge.size()) {
    std::cout << "Type and charge sizes are inconsistent.\n";
    exit(1);
  }
  if (N * 9 != bec.size()) {
    std::cout << "Type and BEC sizes are inconsistent.\n";
    exit(1);
  }

  allocate_memory(N);

  for (std::size_t n = 0; n < potential.size(); ++n) {
    potential[n] = 0.0;
  }
  for (std::size_t n = 0; n < force.size(); ++n) {
    force[n] = 0.0;
  }
  for (std::size_t n = 0; n < virial.size(); ++n) {
    virial[n] = 0.0;
  }
  for (std::size_t n = 0; n < charge.size(); ++n) {
    charge[n] = 0.0;
  }
  for (std::size_t n = 0; n < bec.size(); ++n) {
    bec[n] = 0.0;
  }

  find_neighbor_list_small_box(
    paramb.rc_radial_max, paramb.rc_angular_max, N, MN, box, position, num_cells, ebox, NN_radial, NL_radial,
    NN_angular, NL_angular, r12);

  find_descriptor_small_box(
    true, false, paramb, annmb, N, NN_radial.data(), NL_radial.data(),
    NN_angular.data(), NL_angular.data(), type.data(), r12.data(), r12.data() + size_x12,
    r12.data() + size_x12 * 2, r12.data() + size_x12 * 3, r12.data() + size_x12 * 4,
    r12.data() + size_x12 * 5,
    Fp.data(), sum_fxyz.data(), charge.data(), charge_derivative.data(), potential.data(), nullptr);

  zero_total_charge(N, charge.data());

  find_bec_diagonal(N, charge.data(), bec.data());
  find_bec_radial_small_box(
    paramb,
    annmb,
    N,
    NN_radial.data(),
    NL_radial.data(),
    type.data(),
    r12.data(),
    r12.data() + size_x12,
    r12.data() + size_x12 * 2,
    charge_derivative.data(),
    bec.data());
  find_bec_angular_small_box(
    paramb,
    annmb,
    N,
    NN_angular.data(),
    NL_angular.data(),
    type.data(),
    r12.data() + size_x12 * 3,
    r12.data() + size_x12 * 4,
    r12.data() + size_x12 * 5,
    charge_derivative.data(),
    sum_fxyz.data(),
    bec.data());
  scale_bec(N, annmb.sqrt_epsilon_inf, bec.data());

  if (paramb.charge_mode == 1 || paramb.charge_mode == 2) {
    ewald.find_force(
      N,
      box.data(),
      charge,
      position,
      D_real,
      force,
      virial,
      potential);
  }

  if (paramb.charge_mode == 1) {
    find_force_charge_real_space_small_box(
      N,
      charge_para,
      NN_radial.data(),
      NL_radial.data(),
      charge.data(),
      r12.data(),
      r12.data() + size_x12,
      r12.data() + size_x12 * 2,
      force.data(),
      force.data() + N,
      force.data() + N * 2,
      virial.data(),
      potential.data(),
      D_real.data());
  }

  if (paramb.charge_mode == 3) {
    find_force_charge_real_space_only_small_box(
      N,
      charge_para,
      NN_radial.data(),
      NL_radial.data(),
      charge.data(),
      r12.data(),
      r12.data() + size_x12,
      r12.data() + size_x12 * 2,
      force.data(),
      force.data() + N,
      force.data() + N * 2,
      virial.data(),
      potential.data(),
      D_real.data());
  }

  find_force_radial_small_box(
    paramb, annmb, N, NN_radial.data(), NL_radial.data(), type.data(), r12.data(),
    r12.data() + size_x12, r12.data() + size_x12 * 2, Fp.data(),
    charge_derivative.data(), D_real.data(),
    force.data(), force.data() + N, force.data() + N * 2, virial.data());

  find_force_angular_small_box(
    paramb, annmb, N, NN_angular.data(), NL_angular.data(), type.data(),
    r12.data() + size_x12 * 3, r12.data() + size_x12 * 4, r12.data() + size_x12 * 5, Fp.data(),
    charge_derivative.data(), D_real.data(), sum_fxyz.data(),
    force.data(), force.data() + N, force.data() + N * 2, virial.data());

  if (zbl.enabled) {
    find_force_ZBL_small_box(
      N, paramb, zbl, NN_angular.data(), NL_angular.data(), type.data(), r12.data() + size_x12 * 3,
      r12.data() + size_x12 * 4, r12.data() + size_x12 * 5, force.data(), force.data() + N,
      force.data() + N * 2, virial.data(), potential.data());
  }
}

void NEP::compute_with_dftd3(
  const std::string& xc,
  const double rc_potential,
  const double rc_coordination_number,
  const std::vector<int>& type,
  const std::vector<double>& box,
  const std::vector<double>& position,
  std::vector<double>& potential,
  std::vector<double>& force,
  std::vector<double>& virial)
{
  compute(type, box, position, potential, force, virial);
  const std::size_t N = type.size();
  const std::size_t size_x12 = N * MN;
  set_dftd3_para_all(xc, rc_potential, rc_coordination_number);

  find_neighbor_list_small_box(
    dftd3.rc_radial, dftd3.rc_angular, N, MN, box, position, num_cells, ebox, NN_radial, NL_radial,
    NN_angular, NL_angular, r12);
  find_dftd3_coordination_number(
    dftd3, N, NN_angular.data(), NL_angular.data(), type.data(), r12.data() + size_x12 * 3,
    r12.data() + size_x12 * 4, r12.data() + size_x12 * 5);
  add_dftd3_force(
    dftd3, N, NN_radial.data(), NL_radial.data(), type.data(), r12.data() + size_x12 * 0,
    r12.data() + size_x12 * 1, r12.data() + size_x12 * 2, potential.data(), force.data(),
    virial.data());
  add_dftd3_force_extra(
    dftd3, N, NN_angular.data(), NL_angular.data(), type.data(), r12.data() + size_x12 * 3,
    r12.data() + size_x12 * 4, r12.data() + size_x12 * 5, force.data(), virial.data());
}

void NEP::compute_dftd3(
  const std::string& xc,
  const double rc_potential,
  const double rc_coordination_number,
  const std::vector<int>& type,
  const std::vector<double>& box,
  const std::vector<double>& position,
  std::vector<double>& potential,
  std::vector<double>& force,
  std::vector<double>& virial)
{
  if (paramb.model_type != 0) {
    std::cout << "Cannot compute potential using a non-potential NEP model.\n";
    exit(1);
  }

  const std::size_t N = type.size();
  const std::size_t size_x12 = N * MN;

  if (N * 3 != position.size()) {
    std::cout << "Type and position sizes are inconsistent.\n";
    exit(1);
  }
  if (N != potential.size()) {
    std::cout << "Type and potential sizes are inconsistent.\n";
    exit(1);
  }
  if (N * 3 != force.size()) {
    std::cout << "Type and force sizes are inconsistent.\n";
    exit(1);
  }
  if (N * 9 != virial.size()) {
    std::cout << "Type and virial sizes are inconsistent.\n";
    exit(1);
  }

  allocate_memory(N);

  for (std::size_t n = 0; n < potential.size(); ++n) {
    potential[n] = 0.0;
  }
  for (std::size_t n = 0; n < force.size(); ++n) {
    force[n] = 0.0;
  }
  for (std::size_t n = 0; n < virial.size(); ++n) {
    virial[n] = 0.0;
  }

  set_dftd3_para_all(xc, rc_potential, rc_coordination_number);

  find_neighbor_list_small_box(
    dftd3.rc_radial, dftd3.rc_angular, N, MN, box, position, num_cells, ebox, NN_radial, NL_radial,
    NN_angular, NL_angular, r12);
  find_dftd3_coordination_number(
    dftd3, N, NN_angular.data(), NL_angular.data(), type.data(), r12.data() + size_x12 * 3,
    r12.data() + size_x12 * 4, r12.data() + size_x12 * 5);
  add_dftd3_force(
    dftd3, N, NN_radial.data(), NL_radial.data(), type.data(), r12.data() + size_x12 * 0,
    r12.data() + size_x12 * 1, r12.data() + size_x12 * 2, potential.data(), force.data(),
    virial.data());
  add_dftd3_force_extra(
    dftd3, N, NN_angular.data(), NL_angular.data(), type.data(), r12.data() + size_x12 * 3,
    r12.data() + size_x12 * 4, r12.data() + size_x12 * 5, force.data(), virial.data());
}

void NEP::find_descriptor(
  const std::vector<int>& type,
  const std::vector<double>& box,
  const std::vector<double>& position,
  std::vector<double>& descriptor)
{
  const std::size_t N = type.size();
  const std::size_t size_x12 = N * MN;

  if (N * 3 != position.size()) {
    std::cout << "Type and position sizes are inconsistent.\n";
    exit(1);
  }
  if (N * annmb.dim != descriptor.size()) {
    std::cout << "Type and descriptor sizes are inconsistent.\n";
    exit(1);
  }

  allocate_memory(N);

  find_neighbor_list_small_box(
    paramb.rc_radial_max, paramb.rc_angular_max, N, MN, box, position, num_cells, ebox, NN_radial, NL_radial,
    NN_angular, NL_angular, r12);

  if (paramb.charge_mode > 0) {
    find_descriptor_small_box(
      false, true, paramb, annmb, N, NN_radial.data(), NL_radial.data(),
      NN_angular.data(), NL_angular.data(), type.data(), r12.data(), r12.data() + size_x12,
      r12.data() + size_x12 * 2, r12.data() + size_x12 * 3, r12.data() + size_x12 * 4,
      r12.data() + size_x12 * 5,
      Fp.data(), sum_fxyz.data(), nullptr, nullptr, nullptr, descriptor.data());
  } else {
    find_descriptor_small_box(
      false, true, false, false, paramb, annmb, N, NN_radial.data(), NL_radial.data(),
      NN_angular.data(), NL_angular.data(), type.data(), r12.data(), r12.data() + size_x12,
      r12.data() + size_x12 * 2, r12.data() + size_x12 * 3, r12.data() + size_x12 * 4,
      r12.data() + size_x12 * 5,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
      gn_radial.data(), gn_angular.data(),
#endif
      Fp.data(), sum_fxyz.data(), nullptr, descriptor.data(), nullptr, nullptr, false, nullptr);
  }
}

void NEP::find_latent_space(
  const std::vector<int>& type,
  const std::vector<double>& box,
  const std::vector<double>& position,
  std::vector<double>& latent_space)
{
  const std::size_t N = type.size();
  const std::size_t size_x12 = N * MN;

  if (N * 3 != position.size()) {
    std::cout << "Type and position sizes are inconsistent.\n";
    exit(1);
  }
  if (N * annmb.num_neurons1 != latent_space.size()) {
    std::cout << "Type and latent_space sizes are inconsistent.\n";
    exit(1);
  }

  allocate_memory(N);

  find_neighbor_list_small_box(
    paramb.rc_radial_max, paramb.rc_angular_max, N, MN, box, position, num_cells, ebox, NN_radial, NL_radial,
    NN_angular, NL_angular, r12);

  find_descriptor_small_box(
    false, false, true, false, paramb, annmb, N, NN_radial.data(), NL_radial.data(),
    NN_angular.data(), NL_angular.data(), type.data(), r12.data(), r12.data() + size_x12,
    r12.data() + size_x12 * 2, r12.data() + size_x12 * 3, r12.data() + size_x12 * 4,
    r12.data() + size_x12 * 5,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gn_radial.data(), gn_angular.data(),
#endif
    Fp.data(), sum_fxyz.data(), nullptr, nullptr, latent_space.data(), nullptr, false, nullptr);
}

void NEP::find_B_projection(
  const std::vector<int>& type,
  const std::vector<double>& box,
  const std::vector<double>& position,
  std::vector<double>& B_projection)
{
  const std::size_t N = type.size();
  const std::size_t size_x12 = N * MN;

  if (N * 3 != position.size()) {
    std::cout << "Type and position sizes are inconsistent.\n";
    exit(1);
  }
  if (N * annmb.num_neurons1 * (annmb.dim + 2) != B_projection.size()) {
    std::cout << "Type and B_projection sizes are inconsistent.\n";
    exit(1);
  }

  allocate_memory(N);
  find_neighbor_list_small_box(
    paramb.rc_radial_max, paramb.rc_angular_max, N, MN, box, position, num_cells, ebox, NN_radial, NL_radial,
    NN_angular, NL_angular, r12);

  find_descriptor_small_box(
    false, false, false, false, paramb, annmb, N, NN_radial.data(), NL_radial.data(),
    NN_angular.data(), NL_angular.data(), type.data(), r12.data(), r12.data() + size_x12,
    r12.data() + size_x12 * 2, r12.data() + size_x12 * 3, r12.data() + size_x12 * 4,
    r12.data() + size_x12 * 5,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gn_radial.data(), gn_angular.data(),
#endif
    Fp.data(), sum_fxyz.data(), nullptr, nullptr, nullptr, nullptr, true, B_projection.data());
}

void NEP::find_dipole(
  const std::vector<int>& type,
  const std::vector<double>& box,
  const std::vector<double>& position,
  std::vector<double>& dipole)
{
  if (paramb.model_type != 1) {
    std::cout << "Cannot compute dipole using a non-dipole NEP model.\n";
    exit(1);
  }

  const std::size_t N = type.size();
  const std::size_t size_x12 = N * MN;

  if (N * 3 != position.size()) {
    std::cout << "Type and position sizes are inconsistent.\n";
    exit(1);
  }

  allocate_memory(N);
  std::vector<double> potential(N);  // not used but needed for find_descriptor_small_box
  std::vector<double> virial(N * 3); // need the 3 diagonal components only

  for (std::size_t n = 0; n < potential.size(); ++n) {
    potential[n] = 0.0;
  }
  for (std::size_t n = 0; n < virial.size(); ++n) {
    virial[n] = 0.0;
  }

  find_neighbor_list_small_box(
    paramb.rc_radial_max, paramb.rc_angular_max, N, MN, box, position, num_cells, ebox, NN_radial, NL_radial,
    NN_angular, NL_angular, r12);

  find_descriptor_small_box(
    true, false, false, false, paramb, annmb, N, NN_radial.data(), NL_radial.data(),
    NN_angular.data(), NL_angular.data(), type.data(), r12.data(), r12.data() + size_x12,
    r12.data() + size_x12 * 2, r12.data() + size_x12 * 3, r12.data() + size_x12 * 4,
    r12.data() + size_x12 * 5,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gn_radial.data(), gn_angular.data(),
#endif
    Fp.data(), sum_fxyz.data(), potential.data(), nullptr, nullptr, nullptr, false, nullptr);

  find_force_radial_small_box(
    true, paramb, annmb, N, NN_radial.data(), NL_radial.data(), type.data(), r12.data(),
    r12.data() + size_x12, r12.data() + size_x12 * 2, Fp.data(),
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gnp_radial.data(),
#endif
    nullptr, nullptr, nullptr, virial.data());

  find_force_angular_small_box(
    true, paramb, annmb, N, NN_angular.data(), NL_angular.data(), type.data(),
    r12.data() + size_x12 * 3, r12.data() + size_x12 * 4, r12.data() + size_x12 * 5, Fp.data(),
    sum_fxyz.data(),
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gn_angular.data(), gnp_angular.data(),
#endif
    nullptr, nullptr, nullptr, virial.data());

  for (int d = 0; d < 3; ++d) {
    dipole[d] = 0.0;
    for (std::size_t n = 0; n < N; ++n) {
      dipole[d] += virial[d * N + n];
    }
  }
}

void NEP::find_polarizability(
  const std::vector<int>& type,
  const std::vector<double>& box,
  const std::vector<double>& position,
  std::vector<double>& polarizability)
{
  if (paramb.model_type != 2) {
    std::cout << "Cannot compute polarizability using a non-polarizability NEP model.\n";
    exit(1);
  }

  const std::size_t N = type.size();
  const std::size_t size_x12 = N * MN;

  if (N * 3 != position.size()) {
    std::cout << "Type and position sizes are inconsistent.\n";
    exit(1);
  }

  allocate_memory(N);
  std::vector<double> potential(N);  // not used but needed for find_descriptor_small_box
  std::vector<double> virial(N * 9); // per-atom polarizability

  for (std::size_t n = 0; n < potential.size(); ++n) {
    potential[n] = 0.0;
  }
  for (std::size_t n = 0; n < virial.size(); ++n) {
    virial[n] = 0.0;
  }

  find_neighbor_list_small_box(
    paramb.rc_radial_max, paramb.rc_angular_max, N, MN, box, position, num_cells, ebox, NN_radial, NL_radial,
    NN_angular, NL_angular, r12);

  find_descriptor_small_box(
    true, false, false, true, paramb, annmb, N, NN_radial.data(), NL_radial.data(),
    NN_angular.data(), NL_angular.data(), type.data(), r12.data(), r12.data() + size_x12,
    r12.data() + size_x12 * 2, r12.data() + size_x12 * 3, r12.data() + size_x12 * 4,
    r12.data() + size_x12 * 5,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gn_radial.data(), gn_angular.data(),
#endif
    Fp.data(), sum_fxyz.data(), potential.data(), nullptr, nullptr, virial.data(), false, nullptr);

  find_force_radial_small_box(
    false, paramb, annmb, N, NN_radial.data(), NL_radial.data(), type.data(), r12.data(),
    r12.data() + size_x12, r12.data() + size_x12 * 2, Fp.data(),
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gnp_radial.data(),
#endif
    nullptr, nullptr, nullptr, virial.data());

  find_force_angular_small_box(
    false, paramb, annmb, N, NN_angular.data(), NL_angular.data(), type.data(),
    r12.data() + size_x12 * 3, r12.data() + size_x12 * 4, r12.data() + size_x12 * 5, Fp.data(),
    sum_fxyz.data(),
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gn_angular.data(), gnp_angular.data(),
#endif
    nullptr, nullptr, nullptr, virial.data());

  for (int d = 0; d < 6; ++d) {
    polarizability[d] = 0.0;
  }
  for (std::size_t n = 0; n < N; ++n) {
    polarizability[0] += virial[0 * N + n]; // xx
    polarizability[1] += virial[4 * N + n]; // yy
    polarizability[2] += virial[8 * N + n]; // zz
    polarizability[3] += virial[1 * N + n]; // xy
    polarizability[4] += virial[5 * N + n]; // yz
    polarizability[5] += virial[6 * N + n]; // zx
  }
}

void NEP::compute_for_lammps(
  int nlocal,
  int N,
  int* ilist,
  int* NN,
  int** NL,
  int* type,
  int* type_map,
  double** pos,
  double& total_potential,
  double total_virial[6],
  double* potential,
  double** force,
  double** virial)
{
  if (num_atoms < nlocal) {
    Fp.resize(nlocal * annmb.dim);
    sum_fxyz.resize(nlocal * (paramb.n_max_angular + 1) * NUM_OF_ABC);
    num_atoms = nlocal;
  }
  find_descriptor_for_lammps(
    paramb, annmb, nlocal, N, ilist, NN, NL, type, type_map, pos,
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gn_radial.data(), gn_angular.data(),
#endif
    Fp.data(), sum_fxyz.data(), total_potential, potential);
  find_force_radial_for_lammps(
    paramb, annmb, nlocal, N, ilist, NN, NL, type, type_map, pos, Fp.data(),
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gnp_radial.data(),
#endif
    force, total_virial, virial);
  find_force_angular_for_lammps(
    paramb, annmb, nlocal, N, ilist, NN, NL, type, type_map, pos, Fp.data(), sum_fxyz.data(),
#ifdef USE_TABLE_FOR_RADIAL_FUNCTIONS
    gn_angular.data(), gnp_angular.data(),
#endif
    force, total_virial, virial);
  if (zbl.enabled) {
    find_force_ZBL_for_lammps(
      paramb, zbl, N, ilist, NN, NL, type, type_map, pos, force, total_virial, virial,
      total_potential, potential);
  }
}

bool NEP::set_dftd3_para_one(
  const std::string& functional_input,
  const std::string& functional_library,
  const double s6,
  const double a1,
  const double s8,
  const double a2)
{
  if (functional_input == functional_library) {
    dftd3.s6 = s6;
    dftd3.a1 = a1;
    dftd3.s8 = s8;
    dftd3.a2 = a2 * dftd3para::Bohr;
    return true;
  }
  return false;
}

void NEP::set_dftd3_para_all(
  const std::string& functional_input,
  const double rc_potential,
  const double rc_coordination_number)
{

  dftd3.rc_radial = rc_potential;
  dftd3.rc_angular = rc_coordination_number;

  std::string functional = functional_input;
  std::transform(functional.begin(), functional.end(), functional.begin(), [](unsigned char c) {
    return std::tolower(c);
  });

  bool valid = false;
  valid = valid || set_dftd3_para_one(functional, "b1b95", 1.000, 0.2092, 1.4507, 5.5545);
  valid = valid || set_dftd3_para_one(functional, "b2gpplyp", 0.560, 0.0000, 0.2597, 6.3332);
  valid = valid || set_dftd3_para_one(functional, "b2plyp", 0.640, 0.3065, 0.9147, 5.0570);
  valid = valid || set_dftd3_para_one(functional, "b3lyp", 1.000, 0.3981, 1.9889, 4.4211);
  valid = valid || set_dftd3_para_one(functional, "b3pw91", 1.000, 0.4312, 2.8524, 4.4693);
  valid = valid || set_dftd3_para_one(functional, "b97d", 1.000, 0.5545, 2.2609, 3.2297);
  valid = valid || set_dftd3_para_one(functional, "bhlyp", 1.000, 0.2793, 1.0354, 4.9615);
  valid = valid || set_dftd3_para_one(functional, "blyp", 1.000, 0.4298, 2.6996, 4.2359);
  valid = valid || set_dftd3_para_one(functional, "bmk", 1.000, 0.1940, 2.0860, 5.9197);
  valid = valid || set_dftd3_para_one(functional, "bop", 1.000, 0.4870, 3.295, 3.5043);
  valid = valid || set_dftd3_para_one(functional, "bp86", 1.000, 0.3946, 3.2822, 4.8516);
  valid = valid || set_dftd3_para_one(functional, "bpbe", 1.000, 0.4567, 4.0728, 4.3908);
  valid = valid || set_dftd3_para_one(functional, "camb3lyp", 1.000, 0.3708, 2.0674, 5.4743);
  valid = valid || set_dftd3_para_one(functional, "dsdblyp", 0.500, 0.0000, 0.2130, 6.0519);
  valid = valid || set_dftd3_para_one(functional, "hcth120", 1.000, 0.3563, 1.0821, 4.3359);
  valid = valid || set_dftd3_para_one(functional, "hf", 1.000, 0.3385, 0.9171, 2.883);
  valid = valid || set_dftd3_para_one(functional, "hse-hjs", 1.000, 0.3830, 2.3100, 5.685);
  valid = valid || set_dftd3_para_one(functional, "lc-wpbe08", 1.000, 0.3919, 1.8541, 5.0897);
  valid = valid || set_dftd3_para_one(functional, "lcwpbe", 1.000, 0.3919, 1.8541, 5.0897);
  valid = valid || set_dftd3_para_one(functional, "m11", 1.000, 0.0000, 2.8112, 10.1389);
  valid = valid || set_dftd3_para_one(functional, "mn12l", 1.000, 0.0000, 2.2674, 9.1494);
  valid = valid || set_dftd3_para_one(functional, "mn12sx", 1.000, 0.0983, 1.1674, 8.0259);
  valid = valid || set_dftd3_para_one(functional, "mpw1b95", 1.000, 0.1955, 1.0508, 6.4177);
  valid = valid || set_dftd3_para_one(functional, "mpwb1k", 1.000, 0.1474, 0.9499, 6.6223);
  valid = valid || set_dftd3_para_one(functional, "mpwlyp", 1.000, 0.4831, 2.0077, 4.5323);
  valid = valid || set_dftd3_para_one(functional, "n12sx", 1.000, 0.3283, 2.4900, 5.7898);
  valid = valid || set_dftd3_para_one(functional, "olyp", 1.000, 0.5299, 2.6205, 2.8065);
  valid = valid || set_dftd3_para_one(functional, "opbe", 1.000, 0.5512, 3.3816, 2.9444);
  valid = valid || set_dftd3_para_one(functional, "otpss", 1.000, 0.4634, 2.7495, 4.3153);
  valid = valid || set_dftd3_para_one(functional, "pbe", 1.000, 0.4289, 0.7875, 4.4407);
  valid = valid || set_dftd3_para_one(functional, "pbe0", 1.000, 0.4145, 1.2177, 4.8593);
  valid = valid || set_dftd3_para_one(functional, "pbe38", 1.000, 0.3995, 1.4623, 5.1405);
  valid = valid || set_dftd3_para_one(functional, "pbesol", 1.000, 0.4466, 2.9491, 6.1742);
  valid = valid || set_dftd3_para_one(functional, "ptpss", 0.750, 0.000, 0.2804, 6.5745);
  valid = valid || set_dftd3_para_one(functional, "pw6b95", 1.000, 0.2076, 0.7257, 6.375);
  valid = valid || set_dftd3_para_one(functional, "pwb6k", 1.000, 0.1805, 0.9383, 7.7627);
  valid = valid || set_dftd3_para_one(functional, "pwpb95", 0.820, 0.0000, 0.2904, 7.3141);
  valid = valid || set_dftd3_para_one(functional, "revpbe", 1.000, 0.5238, 2.3550, 3.5016);
  valid = valid || set_dftd3_para_one(functional, "revpbe0", 1.000, 0.4679, 1.7588, 3.7619);
  valid = valid || set_dftd3_para_one(functional, "revpbe38", 1.000, 0.4309, 1.4760, 3.9446);
  valid = valid || set_dftd3_para_one(functional, "revssb", 1.000, 0.4720, 0.4389, 4.0986);
  valid = valid || set_dftd3_para_one(functional, "rpbe", 1.000, 0.1820, 0.8318, 4.0094);
  valid = valid || set_dftd3_para_one(functional, "rpw86pbe", 1.000, 0.4613, 1.3845, 4.5062);
  valid = valid || set_dftd3_para_one(functional, "scan", 1.000, 0.5380, 0.0000, 5.42);
  valid = valid || set_dftd3_para_one(functional, "sogga11x", 1.000, 0.1330, 1.1426, 5.7381);
  valid = valid || set_dftd3_para_one(functional, "ssb", 1.000, -0.0952, -0.1744, 5.2170);
  valid = valid || set_dftd3_para_one(functional, "tpss", 1.000, 0.4535, 1.9435, 4.4752);
  valid = valid || set_dftd3_para_one(functional, "tpss0", 1.000, 0.3768, 1.2576, 4.5865);
  valid = valid || set_dftd3_para_one(functional, "tpssh", 1.000, 0.4529, 2.2382, 4.6550);
  valid = valid || set_dftd3_para_one(functional, "b2kplyp", 0.64, 0.0000, 0.1521, 7.1916);
  valid = valid || set_dftd3_para_one(functional, "dsd-pbep86", 0.418, 0.0000, 0.0000, 5.6500);
  valid = valid || set_dftd3_para_one(functional, "b97m", 1.0000, -0.0780, 0.1384, 5.5946);
  valid = valid || set_dftd3_para_one(functional, "wb97x", 1.0000, 0.0000, 0.2641, 5.4959);
  valid = valid || set_dftd3_para_one(functional, "wb97m", 1.0000, 0.5660, 0.3908, 3.1280);

  if (!valid) {
    std::cout << "The " << functional
              << " functional is not supported for DFT-D3 with BJ damping.\n"
              << std::endl;
    exit(1);
  }
};
