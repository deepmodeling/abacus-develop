#ifndef DEEPKS_VDELTA_H
#define DEEPKS_VDELTA_H

#ifdef __DEEPKS
#include "module_basis/module_ao/parallel_orbitals.h"

namespace DeePKS_domain
{
//------------------------
// deepks_vdelta.cpp
//------------------------

// This file contains 1 subroutine for calculating e_delta_bands
// 1. cal_e_delta_band : calculates e_delta_bands

/// calculate tr(\rho V_delta)
template <typename TK>
void cal_e_delta_band(const std::vector<std::vector<TK>>& dm,
                      const std::vector<std::vector<TK>>& H_V_delta,
                      const int nks,
                      const Parallel_Orbitals* pv,
                      double& e_delta_band);
} // namespace DeePKS_domain
#endif
#endif