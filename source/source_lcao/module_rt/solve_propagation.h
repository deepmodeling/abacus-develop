#ifndef TD_SOLVE_PROPAGATION_H
#define TD_SOLVE_PROPAGATION_H

#include "source_base/matrix2d_view.h"

#include <complex>

namespace module_rt
{
#ifdef __MPI
/**
*  @brief solve propagation equation A@c(t+dt) = B@c(t)
*
* @param[in] dt time interval
* @param[in] overlap overlap matrix S(t+dt/2)
* @param[in] hamiltonian H(t+dt/2)
* @param[in] psi_last psi of last step
* @param[out] psi psi of this step
*/
void solve_propagation(
    const double dt,
    Matrix2dView<const std::complex<double>> overlap,
    Matrix2dView<const std::complex<double>> hamiltonian,
    Matrix2dView<const std::complex<double>> psi_last,
    Matrix2dView<std::complex<double>> psi);

/**
 * @brief solve propagation equation A@c(t+dt) = B@c(t) with moving spatial gauge P_k
 *
 * @param[in] gauge moving spatial gauge matrix
 */
void solve_propagation(
    const double dt,
    Matrix2dView<const std::complex<double>> overlap,
    Matrix2dView<const std::complex<double>> hamiltonian,
    Matrix2dView<const std::complex<double>> gauge,
    Matrix2dView<const std::complex<double>> psi_last,
    Matrix2dView<std::complex<double>> psi);
#endif
} // namespace module_rt

#endif // TD_SOLVE_PROPAGATION_H
