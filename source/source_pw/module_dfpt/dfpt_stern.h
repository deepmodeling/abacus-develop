// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_STERN_H
#define DFPT_STERN_H

#include <complex>
#include <vector>

namespace ModuleDFPT {

/**
 * @brief Projected conjugate-gradient solver of the Sternheimer equation (C2).
 *
 * Solves the insulating valence response
 *   (H(k+q) - eps_n) P_c |dpsi_n> = -P_c |dV psi_n>,
 *   P_c = 1 - sum_{m occ} |u_m(k+q)><u_m(k+q)|,
 * entirely inside the complement of the occupied states at k+q (the metallic
 * occupation / dmu branch is reserved for C4).
 *
 * The shifted Hamiltonian action is injected through LinearOperator so the
 * solver core stays decoupled from the ground-state operator chain: the
 * production adapter reuses hamilt::Hamilt::ops->hPsi at the k+q point
 * (wired in C7), while unit tests supply analytic operators.
 */
class DFPT_Stern {
public:
    DFPT_Stern();
    ~DFPT_Stern();

    /// Hermitian linear action y = (H(k+q) - eps) x on the k+q basis; the
    /// eigenvalue shift is carried inside the implementation.
    class LinearOperator {
    public:
        virtual ~LinearOperator() = default;
        virtual int dimension() const = 0;
        virtual void apply(const std::complex<double>* x, std::complex<double>* y) const = 0;
    };

    /**
     * @brief Solve one projected Sternheimer system with conjugate gradients.
     *
     * @param aop       shifted Hamiltonian (H(k+q) - eps_n), Hermitian
     * @param occ_kq    orthonormal occupied states at k+q (may be empty)
     * @param b         right-hand side -dV|psi_n>; projected internally
     * @param max_iter  linear-solver iteration cap (> 0)
     * @param conv_thr  relative residual threshold ||P_c r|| / ||P_c b||
     * @param dpsi      output P_c|dpsi_n> (zero inside the occ subspace)
     * @param residual  achieved relative residual
     * @return          number of iterations used
     */
    int solve(const LinearOperator& aop,
              const std::vector<std::vector<std::complex<double>>>& occ_kq,
              const std::vector<std::complex<double>>& b,
              int max_iter,
              double conv_thr,
              std::vector<std::complex<double>>& dpsi,
              double& residual) const;

private:
    /// P_c x by modified Gram-Schmidt against the occupied states; safe for
    /// px to alias x (projection coefficients are collected before subtracting)
    void apply_pv(const std::vector<std::vector<std::complex<double>>>& occ_kq,
                  const std::vector<std::complex<double>>& x,
                  std::vector<std::complex<double>>& px) const;
};

} // namespace ModuleDFPT

#endif // DFPT_STERN_H
