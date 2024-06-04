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
#ifndef LOWF_NEW
#define LOWF_NEW

// serial
#include <string>
#include <vector>
#include <complex>
#include "module_base/vector3.h"
// parallelization
#ifdef __MPI
#include "mpi.h"
#include "module_basis/module_ao/parallel_2d.h"
#include "module_base/blacs_connector.h"
#endif
/**
 * @brief This class has two functions: restart psi from the previous calculation, and write psi to the disk.
 * 
 */
namespace lowf
{
// function 1: restart psi from the previous calculation (one single kpoint)
// in: flowf, out: lowf(lowf_glb), ekb, occ, kvec_c, wk

void read_abacus_lowf(const std::string& flowf, 
                      int& nbands, 
                      int& nbasis, 
                      std::vector<std::complex<double>>& lowf, 
                      std::vector<double>& ekb, 
                      std::vector<double>& occ,
                      ModuleBase::Vector3<double>& kvec_c, 
                      double& wk);
void read_abacus_lowf(const std::string& flowf, 
                      int& nbands, 
                      int& nbasis, 
                      std::vector<double>& lowf, 
                      std::vector<double>& ekb, 
                      std::vector<double>& occ, 
                      ModuleBase::Vector3<double>& kvec_c, 
                      double& wk);
// the two functions above will return nbands, nbasis, lowf, ekb, occ and wk.
// the lowf is actually lowf_glb, which means the global matrix (ScaLAPACK convention), need to distribute
// to the local matrix (2D-block-cyclic parallel distribution) in the following function.

#ifdef __MPI
// the following function will be used to restart/assign the value of psi_2d
void scatter_lowf(const int& nbands, 
                  const int& nbasis, 
                  const std::vector<std::complex<double>>& lowf_glb,
                  const MPI_Comm& comm,  // MPI comm world for 2d bcd
                  const int& nb,         // 64, 32 or 1, performance tunning param
                  const int& blacs_ctxt, // the one of orb_con.ParaV
                  std::vector<std::complex<double>>& lowf_loc)
{
    Parallel_2D para2d_glb, para2d_loc;
    para2d_glb.init(nbands, nbasis, nbasis, comm);
    para2d_loc.set(nbands, nbasis, nb, comm, blacs_ctxt);
    lowf_loc.resize(nb*nb);
    const int ctxt = Csys2blacs_handle(comm);
}
#endif
// then will construct psi like:
// psi::Psi psi_loc(lowf_loc.data(), nks, nbands, nbasis, ...);

} // namespace lowf
#endif