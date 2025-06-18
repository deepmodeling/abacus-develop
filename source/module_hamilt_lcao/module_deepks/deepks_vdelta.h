#ifndef DEEPKS_VDELTA_H
#define DEEPKS_VDELTA_H

#ifdef __MLALGO
#include "source_base/complexmatrix.h"
#include "source_base/matrix.h"
#include "module_basis/module_ao/parallel_orbitals.h"

// break the circular dependency of HamiltLCAO
namespace hamilt
{
template <typename TK, typename TR>
class HamiltLCAO;
}
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
                      const int nspin,
                      const Parallel_Orbitals* pv,
                      double& e_delta_band);

// Collect data in h_in to matrix h_out. Note that left lower trianger in h_out is filled
template <typename TK, typename TH>
void collect_h_mat(const Parallel_Orbitals& pv,
                   const std::vector<std::vector<TK>>& h_in,
                   std::vector<TH>& h_out,
                   const int nlocal,
                   const int nks);

// Get H(k) or S(k) matrix from p_hamilt and store it in h_tot
template <typename TK, typename TH, typename TR>
void get_h_tot(const Parallel_Orbitals& pv,
               hamilt::HamiltLCAO<TK, TR>* p_ham,
               std::vector<TH>& h_tot,
               const int nlocal,
               const int nks,
               const char matrix_type); // 'H' for H(k), 'S' for S(k)
} // namespace DeePKS_domain
#endif
#endif