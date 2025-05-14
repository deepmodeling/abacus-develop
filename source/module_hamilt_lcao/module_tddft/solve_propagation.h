#ifndef TD_SOLVE_PROPAGATION_H
#define TD_SOLVE_PROPAGATION_H

#include "module_basis/module_ao/parallel_orbitals.h"
#include <complex>

namespace module_tddft
{
#ifdef __MPI
/**
*  @brief solve propagation equation A@c(t+dt) = B@c(t)
*
* @param[in] pv information of parallel
* @param[in] nband number of bands
* @param[in] nlocal number of orbitals
* @param[in] dt time interval
* @param[in] Stmp overlap matrix S(t+dt/2)
* @param[in] Htmp H(t+dt/2)
* @param[in] psi_k_laststep psi of last step
* @param[out] psi_k psi of this step
*/
void solve_propagation(const Parallel_Orbitals* pv,
                        const int nband,
                        const int nlocal,
                        const double dt,
                        const std::complex<double>* Stmp,
                        const std::complex<double>* Htmp,
                        const std::complex<double>* psi_k_laststep,
                        std::complex<double>* psi_k);

#endif
} // namespace module_tddft

#endif // TD_SOLVE_H