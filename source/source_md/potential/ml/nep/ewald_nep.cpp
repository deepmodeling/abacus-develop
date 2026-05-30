/*
    Copyright 2017 Zheyong Fan and GPUMD development team
    This file is part of GPUMD.
    GPUMD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    GPUMD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with GPUMD.  If not, see <http://www.gnu.org/licenses/>.
*/

/*----------------------------------------------------------------------------80
The k-space part of the Ewald summation
------------------------------------------------------------------------------*/

#include "ewald_nep.h"
#include <cmath>
#include <vector>

namespace {

const double K_C_SP_ = 14.399645; // 1/(4*PI*epsilon_0)

void cross_product(const double a[3], const double b[3], double c[3])
{
  c[0] =  a[1] * b [2] - a[2] * b [1];
  c[1] =  a[2] * b [0] - a[0] * b [2];
  c[2] =  a[0] * b [1] - a[1] * b [0];
}

double get_area(const double* a, const double* b)
{
  const double s1 = a[1] * b[2] - a[2] * b[1];
  const double s2 = a[2] * b[0] - a[0] * b[2];
  const double s3 = a[0] * b[1] - a[1] * b[0];
  return sqrt(s1 * s1 + s2 * s2 + s3 * s3);
}

void find_structure_factor(
  const int num_kpoints,
  const int N,
  const double* g_charge,
  const double* g_x,
  const double* g_y,
  const double* g_z,
  const double* g_kx,
  const double* g_ky,
  const double* g_kz,
  double* g_S_real,
  double* g_S_imag)
{
  for (int nk = 0; nk < num_kpoints; ++nk) {
    double S_real = 0.0;
    double S_imag = 0.0;
    for (int n = 0; n < N; ++n) {
      double kr = g_kx[nk] * double(g_x[n]) + g_ky[nk] * double(g_y[n]) + g_kz[nk] * double(g_z[n]);
      const double charge = g_charge[n];
      double sin_kr = sin(kr);
      double cos_kr = cos(kr);
      S_real += charge * cos_kr;
      S_imag -= charge * sin_kr;
    }
    g_S_real[nk] = S_real;
    g_S_imag[nk] = S_imag;
  }
}

void find_force_charge_reciprocal_space(
  const int N,
  const int num_kpoints,
  const double alpha_factor,
  const double* g_charge,
  const double* g_x,
  const double* g_y,
  const double* g_z,
  const double* g_kx,
  const double* g_ky,
  const double* g_kz,
  const double* g_G,
  const double* g_S_real,
  const double* g_S_imag,
  double* g_D_real,
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial,
  double* g_pe)
{
  for (int n = 0; n < N; ++n) {
    const double q = g_charge[n];
    double temp_energy_sum = 0.0;
    double temp_virial_sum[6] = {0.0};
    double temp_force_sum[3] = {0.0};
    double temp_D_real_sum = 0.0;
    for (int nk = 0; nk < num_kpoints; ++nk) {
      const double kx = g_kx[nk];
      const double ky = g_ky[nk];
      const double kz = g_kz[nk];
      const double kr = kx * g_x[n] + ky * g_y[n] + kz * g_z[n];
      const double G = g_G[nk];
      const double S_real = g_S_real[nk];
      const double S_imag = g_S_imag[nk];
      double sin_kr = sin(kr);
      double cos_kr = cos(kr);
      const double imag_term = G * (S_real * sin_kr + S_imag * cos_kr);
      const double GSE = G * (S_real * cos_kr - S_imag * sin_kr);
      const double qGSE = q * GSE;
      temp_energy_sum += qGSE;
      const double alpha_k_factor = 2.0 * alpha_factor + 2.0 / (kx * kx + ky * ky + kz * kz);
      temp_virial_sum[0] += qGSE * (1.0 - alpha_k_factor * kx * kx); // xx
      temp_virial_sum[1] += qGSE * (1.0 - alpha_k_factor * ky * ky); // yy
      temp_virial_sum[2] += qGSE * (1.0 - alpha_k_factor * kz * kz); // zz
      temp_virial_sum[3] -= qGSE * (alpha_k_factor * kx * ky); // xy
      temp_virial_sum[4] -= qGSE * (alpha_k_factor * ky * kz); // yz
      temp_virial_sum[5] -= qGSE * (alpha_k_factor * kz * kx); // zx
      temp_D_real_sum += GSE;
      temp_force_sum[0] += kx * imag_term;
      temp_force_sum[1] += ky * imag_term;
      temp_force_sum[2] += kz * imag_term;
    }
    g_pe[n] += K_C_SP_ * temp_energy_sum;
    g_virial[n + 0 * N] += K_C_SP_ * temp_virial_sum[0];
    g_virial[n + 1 * N] += K_C_SP_ * temp_virial_sum[3];
    g_virial[n + 2 * N] += K_C_SP_ * temp_virial_sum[5];
    g_virial[n + 3 * N] += K_C_SP_ * temp_virial_sum[3];
    g_virial[n + 4 * N] += K_C_SP_ * temp_virial_sum[1];
    g_virial[n + 5 * N] += K_C_SP_ * temp_virial_sum[4];
    g_virial[n + 6 * N] += K_C_SP_ * temp_virial_sum[5];
    g_virial[n + 7 * N] += K_C_SP_ * temp_virial_sum[4];
    g_virial[n + 8 * N] += K_C_SP_ * temp_virial_sum[2];
    g_D_real[n] = 2.0 * K_C_SP_ * temp_D_real_sum;
    const double charge_factor = K_C_SP_ * 2.0 * q;
    g_fx[n] += charge_factor * temp_force_sum[0];
    g_fy[n] += charge_factor * temp_force_sum[1];
    g_fz[n] += charge_factor * temp_force_sum[2];
  }
}

}

EwaldNep::EwaldNep()
{
  // nothing
}

EwaldNep::~EwaldNep()
{
  // nothing
}

void EwaldNep::initialize(const double alpha_input)
{
  alpha = alpha_input;
  alpha_factor = 0.25 / (alpha * alpha);
  kx.resize(num_kpoints_max);
  ky.resize(num_kpoints_max);
  kz.resize(num_kpoints_max);
  G.resize(num_kpoints_max);
  S_real.resize(num_kpoints_max);
  S_imag.resize(num_kpoints_max);
}

void EwaldNep::find_k_and_G(const double* box)
{
  double a1[3] = {0.0};
  double a2[3] = {0.0};
  double a3[3] = {0.0};
  double det = box[0] * (box[4] * box[8] - box[5] * box[7]) +
    box[1] * (box[5] * box[6] - box[3] * box[8]) +
    box[2] * (box[3] * box[7] - box[4] * box[6]);
  a1[0] = box[0];
  a1[1] = box[3];
  a1[2] = box[6];
  a2[0] = box[1];
  a2[1] = box[4];
  a2[2] = box[7];
  a3[0] = box[2];
  a3[1] = box[5];
  a3[2] = box[8];
  double b1[3] = {0.0};
  double b2[3] = {0.0};
  double b3[3] = {0.0};
  cross_product(a2, a3, b1);
  cross_product(a3, a1, b2);
  cross_product(a1, a2, b3);

  const double two_pi = 6.2831853;
  const double two_pi_over_det = two_pi / det;
  for (int d = 0; d < 3; ++d) {
    b1[d] *= two_pi_over_det;
    b2[d] *= two_pi_over_det;
    b3[d] *= two_pi_over_det;
  }

  const double volume_k = two_pi * two_pi * two_pi / std::abs(det);
  int n1_max = alpha * two_pi * get_area(b2, b3) / volume_k;
  int n2_max = alpha * two_pi * get_area(b3, b1) / volume_k;
  int n3_max = alpha * two_pi * get_area(b1, b2) / volume_k;
  double ksq_max = two_pi * two_pi * alpha * alpha;

  kx.clear();
  ky.clear();
  kz.clear();
  G.clear();

  for (int n1 = 0; n1 <= n1_max; ++n1) {
    for (int n2 = - n2_max; n2 <= n2_max; ++n2) {
      for (int n3 = - n3_max; n3 <= n3_max; ++n3) {
        const int nsq = n1 * n1 + n2 * n2 + n3 * n3;
        if (nsq == 0 || (n1 == 0 && n2 < 0) || (n1 == 0 && n2 == 0 && n3 < 0)) continue;
        const double kx_temp = n1 * b1[0] + n2 * b2[0] + n3 * b3[0];
        const double ky_temp = n1 * b1[1] + n2 * b2[1] + n3 * b3[1];
        const double kz_temp = n1 * b1[2] + n2 * b2[2] + n3 * b3[2];
        const double ksq = kx_temp * kx_temp + ky_temp * ky_temp + kz_temp * kz_temp;
        if (ksq < ksq_max) {
          kx.emplace_back(kx_temp);
          ky.emplace_back(ky_temp);
          kz.emplace_back(kz_temp);
          G.emplace_back(2.0 * fabs(two_pi_over_det) / ksq * exp(-ksq * alpha_factor));
        }
      }
    }
  }

  num_kpoints = int(G.size());

  if (num_kpoints > num_kpoints_max) {
    num_kpoints_max = num_kpoints;
    S_real.resize(num_kpoints_max);
    S_imag.resize(num_kpoints_max);
  }
}

void EwaldNep::find_force(
  const int N,
  const double* box,
  const std::vector<double>& charge,
  const std::vector<double>& position_per_atom,
  std::vector<double>& D_real,
  std::vector<double>& force_per_atom,
  std::vector<double>& virial_per_atom,
  std::vector<double>& potential_per_atom)
{
  find_k_and_G(box);
  find_structure_factor(
    num_kpoints,
    N,
    charge.data(),
    position_per_atom.data(),
    position_per_atom.data() + N,
    position_per_atom.data() + N * 2,
    kx.data(),
    ky.data(),
    kz.data(),
    S_real.data(),
    S_imag.data());

  find_force_charge_reciprocal_space(
    N,
    num_kpoints,
    alpha_factor,
    charge.data(),
    position_per_atom.data(),
    position_per_atom.data() + N,
    position_per_atom.data() + N * 2,
    kx.data(),
    ky.data(),
    kz.data(),
    G.data(),
    S_real.data(),
    S_imag.data(),
    D_real.data(),
    force_per_atom.data(),
    force_per_atom.data() + N,
    force_per_atom.data() + N * 2,
    virial_per_atom.data(),
    potential_per_atom.data());
}
