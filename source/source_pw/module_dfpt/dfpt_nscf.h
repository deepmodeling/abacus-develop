#ifndef DFPT_NSCF_H
#define DFPT_NSCF_H

#include "dfpt_stern.h"

#include <complex>
#include <vector>

namespace ModuleDFPT
{

/**
 * @brief Block-Davidson eigensolver for the full k+q Hamiltonian.
 *
 * QE aligns the k+q occupied states with an nscf diagonalization
 * (run_nscf -> non_scf_ph) rather than any copy from the ground-state
 * list; this class provides the ABACUS counterpart on top of the
 * DFPT_Stern::LinearOperator interface (DFPT_HamiltShift with shift 0
 * applies the full H(k+q)).
 *
 * The iteration keeps a subspace V (n x m), Rayleigh-Ritz projects with
 * zheev, and expands with the diagonal-preconditioned residuals; the
 * subspace is compressed back to the current Ritz vectors when it grows
 * past a multiple of nwant. The preconditioner is the kinetic diagonal
 * tpiba^2 |G+k+q|^2 (+ a small floor), the standard highest-kinetic
 * filter for plane-wave Davidson.
 */
class DFPT_NSCF
{
  public:
    /**
     * @brief Solve the lowest nwant eigenpairs of the operator.
     * @param hop  full-Hamiltonian operator (apply must be H, no shift)
     * @param precond diagonal preconditioner, size = hop.dimension()
     * @param nwant number of lowest eigenpairs requested (<= dimension)
     * @param max_iter outer Davidson iterations
     * @param eig_tol convergence threshold on the residual norm (Ry)
     * @param psi output eigenvectors [nwant][n]
     * @param eig output eigenvalues ascending [nwant]
     * @return 0 converged, >0 iterations used without full convergence
     */
    static int lowest_eigpairs(const DFPT_Stern::LinearOperator& hop,
                               const std::vector<double>& precond,
                               int nwant,
                               int max_iter,
                               double eig_tol,
                               std::vector<std::vector<std::complex<double>>>& psi,
                               std::vector<double>& eig);
};

} // namespace ModuleDFPT

#endif // DFPT_NSCF_H
