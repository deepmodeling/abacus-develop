/**
 * @file lowf_new.h
 * @author kirk0830 (you@domain.com)
 * @brief This file belongs to the project of removal of class Local_Orbital_Wfc (LOWF)
 * @version 0.1
 * @date 2024-06-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef READ_WFC_LCAO_H
#define READ_WFC_LCAO_H

// serial
#include <string>
#include <vector>
#include <complex>
#include "module_base/vector3.h"
#ifdef __MPI
// parallelization
#include "module_basis/module_ao/parallel_2d.h"
#include "module_base/scalapack_connector.h"
#endif

/**
 * @brief This class has two functions: restart psi from the previous calculation, and write psi to the disk.
 * 
 */
namespace ModuleIO
{
// function 1: restart psi from the previous calculation (one single kpoint)
// in: flowf, 
// out: lowf(lowf_glb), ekb, occ, kvec_c, wk
// only rank 0 will do this, but for larger wavefunction, it is needed to seperately
// read the file
template<typename T>
void read_abacus_lowf(const std::string& flowf, 
                      int& ik,
                      ModuleBase::Vector3<double>& kvec_c, 
                      int& nbands, 
                      int& nbasis, 
                      std::vector<std::complex<T>>& lowf, 
                      std::vector<double>& ekb, 
                      std::vector<double>& occ,
                      double& wk);
template<typename T>
void read_abacus_lowf(const std::string& flowf, 
                      int& ik,
                      ModuleBase::Vector3<double>& kvec_c, 
                      int& nbands, 
                      int& nbasis, 
                      std::vector<T>& lowf, 
                      std::vector<double>& ekb, 
                      std::vector<double>& occ,
                      double& wk);
// the two functions above will return nbands, nbasis, lowf, ekb, occ and wk.
// the lowf is actually lowf_glb, which means the global matrix (ScaLAPACK convention), need to distribute
// to the local matrix (2D-block-cyclic parallel distribution) in the following function.

// only-MPI-visible function, because the use of comm_world
#ifdef __MPI
template <typename T>
void restart_from_file(const std::string& out_dir, // hard-code the file name to be LOWF_K_*.txt?
                       const Parallel_2D& p2d,
                       const int& nks,
                       int& nbands,
                       int& nbasis,
                       std::vector<T>& lowf_loc,
                       std::vector<double>& ekb,
                       std::vector<double>& occ,
                       std::vector<ModuleBase::Vector3<double>>& kvec_c,
                       std::vector<double>& wk);
#endif
// serial version, can always present
template <typename T>
void restart_from_file(const std::string& out_dir, // hard-code the file name to be LOWF_K_*.txt?
                       const int& nks,
                       int& nbands,
                       int& nbasis,
                       std::vector<T>& lowf,
                       std::vector<double>& ekb,
                       std::vector<double>& occ,
                       std::vector<ModuleBase::Vector3<double>>& kvec_c,
                       std::vector<double>& wk);
} // namespace ModuleIO
#endif