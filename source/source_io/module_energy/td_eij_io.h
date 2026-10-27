#ifndef MODULE_IO_TD_EIJ_IO_H
#define MODULE_IO_TD_EIJ_IO_H

#include "source_cell/klist.h"
#include "source_hsolver/module_rt_pw/td_orbital_ops.h"

#include <complex>
#include <ostream>

namespace ModuleIO
{
/** @brief Eij output settings and one-based identifiers for one electronic state. */
struct EijOutputInfo
{
    double threshold;
    int electronic_step;
    int scf_iteration;
    int kpoint;
    int nspin;
    int spin;
};

/** @brief Print selected elements of a column-major Hamiltonian in the orbital basis, in Ry. */
void print_eij(const std::complex<double>* eij,
               const int nband,
               const int leading_dimension,
               const EijOutputInfo& output_info,
               std::ostream& ofs);
/** @brief Collect and print Eij from an already bound endpoint Hamiltonian.
 * All ranks participate; both k-point cursors are restored before returning.
 */
template <typename T, typename Device>
void print_eij_pw(hsolver::HSOperator<T, Device>& op,
                  hsolver::TDOrbOps<T, Device>& orbitals,
                  psi::Psi<T, Device>& current,
                  const K_Vectors& kv,
                  int nspin,
                  int step,
                  int iter,
                  double threshold,
                  std::ostream& log);
} // namespace ModuleIO
#endif
