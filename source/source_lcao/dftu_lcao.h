#ifndef DFTU_LCAO_H
#define DFTU_LCAO_H

#include <string>
#include <vector>
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"

namespace ModuleESolver
{

/**
 * @brief Initialize DFT+U for LCAO method in iter_init
 *
 * This function handles the DFT+U initialization during the SCF iteration.
 * It sets the density matrix and calculates Slater integrals if needed.
 *
 * @param istep Current ionic step
 * @param iter Current SCF iteration
 * @param dft_plus_u DFT+U mode (0=disabled, 1=old, 2=new)
 * @param dftu DFT+U object
 * @param dm Density matrix
 * @param ucell Unit cell
 * @param rho Charge density
 * @param nrxx Number of real space grid points
 */
template <typename TK>
void init_dftu_lcao(const int istep,
                     const int iter,
                     int dft_plus_u,
                     void* dftu,
                     void* dm,
                     const UnitCell& ucell,
                     double** rho,
                     const int nrxx);

/**
 * @brief Finish DFT+U calculation for LCAO method in iter_finish
 *
 * This function handles the DFT+U finalization during the SCF iteration.
 * It calculates the occupation matrix and energy correction if needed.
 *
 * @param iter Current SCF iteration
 * @param conv_esolver Whether ESolver has converged
 * @param dft_plus_u DFT+U mode (0=disabled, 1=old, 2=new)
 * @param out_chg Whether to output dm_onsite.txt
 * @param dftu DFT+U object
 * @param ucell Unit cell
 * @param dm_vec Density matrix vector
 * @param kv K-vectors
 * @param mixing_beta Mixing beta parameter
 * @param hamilt_lcao Hamiltonian LCAO object
 * @param global_out_dir Output directory for dm_onsite.txt
 * @param nspin Number of spin channels (1, 2, or 4)
 * @param npol Number of polarizations
 */
template <typename TK>
void finish_dftu_lcao(const int iter,
                       const bool conv_esolver,
                       int dft_plus_u,
                       bool out_chg,
                       void* dftu,
                       const UnitCell& ucell,
                       const std::vector<std::vector<TK>>& dm_vec,
                       const K_Vectors& kv,
                       const double mixing_beta,
                       void* hamilt_lcao,
                       const std::string& global_out_dir,
                       int nspin,
                       int npol);

} // namespace ModuleESolver

#endif // DFTU_LCAO_H
