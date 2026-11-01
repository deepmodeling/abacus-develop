#ifndef TD_EIJ_PW_H
#define TD_EIJ_PW_H

#include "source_cell/klist.h"
#include "source_pw/module_rt/td_orbital_ops.h"

#include <ostream>

namespace ModuleIO
{
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
