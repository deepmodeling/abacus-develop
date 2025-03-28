#ifndef DEEPKS_VDELTA_H
#define DEEPKS_VDELTA_H

#ifdef __DEEPKS
#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include "module_basis/module_ao/parallel_orbitals.h"

namespace DeePKS_domain
{
//------------------------
// deepks_vdelta.cpp
//------------------------

// This file contains 2 subroutine for calculating e_delta_bands
// 1. cal_e_delta_band : calculates e_delta_bands
// 2. collect_h_mat, which collect H(k) data from different processes

/// calculate tr(\rho V_delta)
template <typename TK>
void cal_e_delta_band(const std::vector<std::vector<TK>>& dm,
                      const std::vector<std::vector<TK>>& V_delta,
                      const int nks,
                      const Parallel_Orbitals* pv,
                      double& e_delta_band);

// Collect data in h_in to matrix h_out. Note that left lower trianger in h_out is filled
template <typename TK, typename TH>
void collect_h_mat(const Parallel_Orbitals& pv,
                    const std::vector<std::vector<TK>>& h_in,
                    std::vector<TH>& h_out,
                    const int nlocal,
                    const int nks);
} // namespace DeePKS_domain
#endif
#endif