#ifndef TD_SOLVE_PROPAGATION_H
#define TD_SOLVE_PROPAGATION_H

#include <complex>

class BlockCyclicDistribution;

namespace module_rt
{
#ifdef __MPI
/**
*  @brief solve propagation equation A@c(t+dt) = B@c(t)
*
* @param[in] matrix_layout block-cyclic layout of the H/S matrices
* @param[in] wfc_layout block-cyclic layout of the wavefunction
* @param[in] nband number of bands
* @param[in] nlocal number of orbitals
* @param[in] dt time interval
* @param[in] Stmp overlap matrix S(t+dt/2)
* @param[in] Htmp H(t+dt/2)
* @param[in] psi_k_laststep psi of last step
* @param[out] psi_k psi of this step
*/
void solve_propagation(const BlockCyclicDistribution& matrix_layout,
                        const BlockCyclicDistribution& wfc_layout,
                        const int nband,
                        const int nlocal,
                        const double dt,
                        const std::complex<double>* Stmp,
                        const std::complex<double>* Htmp,
                        const std::complex<double>* psi_k_laststep,
                        std::complex<double>* psi_k);

/**
 * @brief solve propagation equation A@c(t+dt) = B@c(t) with moving spatial gauge P_k
 *
 * @param[in] P_k moving spatial gauge matrix
 */
void solve_propagation(const BlockCyclicDistribution& matrix_layout,
                       const BlockCyclicDistribution& wfc_layout,
                       const int nband,
                       const int nlocal,
                       const double dt,
                       const std::complex<double>* Stmp,
                       const std::complex<double>* Htmp,
                       const std::complex<double>* P_k,
                       const std::complex<double>* psi_k_laststep,
                       std::complex<double>* psi_k);
#endif
} // namespace module_rt

#endif // TD_SOLVE_H
