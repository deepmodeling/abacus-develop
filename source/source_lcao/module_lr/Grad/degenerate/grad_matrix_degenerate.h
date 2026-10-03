#pragma once
#include <cstddef>
#include <utility>
#include <vector>

/// @file
/// Algebra of the degenerate-subspace gradient matrix
/// $G^{(A\alpha)}_{kl}=\langle X_k|\partial A/\partial R_{A\alpha}|X_l\rangle$.
///
/// At a $d$-fold degeneracy no single state has a gradient vector: the branch slopes along a
/// displacement $u$ are the eigenvalues of $M(u)=\sum_{A\alpha}u_{A\alpha}G^{(A\alpha)}$, whose
/// eigenvectors depend on $u$. The full first-order object is $3N$ matrices of size $d\times d$,
/// and only $\operatorname{Tr}M$ is basis-invariant.
///
/// The per-state analytic gradient implemented in `esolver_lr_grad.cpp` is, restricted to one
/// multiplet, a *quadratic form* in the excitation vector: every ingredient ($D^X$, $T$, $R$, $Z$,
/// $D^Z$, $W^X$, $\Lambda$) is built from $X\otimes X$, the CPSCF operator $L$ depends only on the
/// ground state, and $\Omega$ is a constant inside the multiplet. Write it $\mathcal F[X]$. Because
/// every $X$ in the degenerate subspace is itself a legitimate eigenvector with the same $\Omega$,
/// $\mathcal F[X]=G[X,X]$ holds on the whole diagonal -- and a symmetric bilinear form is
/// determined by its diagonal. So the existing code already contains all of $G$; it has only ever
/// been evaluated on the basis the diagonalizer happened to return.
///
/// Recovering the off-diagonal elements is then the polarization identity. With
/// $X_\pm=(X_k\pm X_l)/\sqrt2$ and $\mathcal F[\alpha X]=\alpha^2\mathcal F[X]$,
///     (A1)  $G_{kl}=\tfrac12(\mathcal F[X_+]-\mathcal F[X_-])$,
///     (A2)  $G_{kl}=\mathcal F[X_+]-\tfrac12(G_{kk}+G_{ll})$,
/// and (A2) reuses the $d$ diagonal gradients the code already computes, so the whole matrix costs
/// $d(d+1)/2$ evaluations -- exactly its number of independent components.
///
/// The construction admits two checks that need no finite differences and are worth running on
/// any new case: $\mathcal F[2X]=4\mathcal F[X]$ (it is a quadratic form at all), and
/// covariance under a rotation of the subspace basis, $\mathcal F[X'_k]=(U^\top GU)_{kk}$ with
/// $X'_k=\sum_lU_{lk}X_l$. Both are exercised in `test/test_grad_matrix_degenerate.cpp`.
///
/// The alternative -- deriving an explicitly bilinear Z-vector right-hand side that takes two
/// different $X$ -- yields the same $G$, but it has to re-derive every factor and hand-polarize
/// the $g^{xc}$ potential, so it is the more expensive and more error-prone route.
///
/// This header holds only the basis-independent bookkeeping, so that it is unit-testable without a
/// ground state: grouping states into multiplets, enumerating the pairs, forming the normalized
/// combination, and assembling $G$. Everything needing the parallel layout or the Z-vector solver
/// stays in `esolver_lr_grad.cpp`.

namespace LR
{
    /// @brief Split states into multiplets of (near-)degenerate excitation energies.
    ///
    /// A state joins the group it is within `thr` of *the group's first member*, not of its
    /// predecessor: chaining on consecutive gaps would let a run of small steps span a spread far
    /// larger than `thr`, and "degenerate" has to mean a bounded total spread.
    ///
    /// The threshold alone does NOT decide whether route (A2) applies: an accidental near-degeneracy
    /// falls inside any loose threshold, yet there $\Omega$ differs between the states, so they are
    /// not one quadratic form and the combinations $X_\pm$ are not eigenvectors. The discriminator
    /// is whether both the analytic and the finite-difference side split by the same amount (see
    /// section 2(B) of the document above); this function only proposes the candidates.
    ///
    /// @param omega  excitation energies, any order (Ry)
    /// @param thr    maximum spread inside one multiplet (Ry); <= 0 puts every state alone
    /// @return  groups of indices into `omega`, each group ascending, groups ordered by their
    ///          lowest energy
    std::vector<std::vector<int>> group_degenerate_states(const std::vector<double>& omega,
        const double thr);

    /// @brief The $(k,l)$, $k<l$ pairs of a $d$-dimensional multiplet, lexicographically.
    ///
    /// This fixes the order `assemble_grad_matrix` expects its `plus` argument in, and hence the
    /// order the combinations have to be handed to the gradient evaluation in.
    std::vector<std::pair<int, int>> degenerate_pairs(const int d);

    /// @brief $X_+=(X_k+X_l)/\sqrt2$, normalized when $X_k$ and $X_l$ are orthonormal.
    ///
    /// Using the normalized combination rather than $X_k+X_l$ is what lets the gradient be
    /// evaluated by the untouched per-state path: $X_+$ is then a genuine normalized eigenvector,
    /// so every normalization, $\Omega$ and $W^X$ convention inside that path still holds.
    template <typename T>
    void combine_normalized(const T* const Xk, const T* const Xl, const size_t nloc, T* const Xplus)
    {
        // 1/sqrt(2) as a literal: `std::sqrt` is not constexpr under the C++11 baseline.
        const T inv_sqrt2 = static_cast<T>(0.70710678118654752440);
        for (size_t i = 0; i < nloc; ++i) { Xplus[i] = inv_sqrt2 * (Xk[i] + Xl[i]); }
    }

    /// @brief Assemble $G_{kl}$ of one multiplet from the diagonal gradients and the combinations,
    ///        i.e. route (A2) above.
    ///
    /// @param diag   $G_{kk}=\mathcal F[X_k]$, `d` entries
    /// @param plus   $\mathcal F[(X_k+X_l)/\sqrt2]$, one per entry of `pairs`, in that order
    /// @param pairs  as returned by `degenerate_pairs(d)`
    /// @return  G[k][l], symmetric by construction
    ///
    /// `TMat` needs `operator+`, `operator-` and `operator*(double)`; `ModuleBase::matrix` and
    /// plain `double` both qualify, which is what keeps this testable.
    template <typename TMat>
    std::vector<std::vector<TMat>> assemble_grad_matrix(const std::vector<TMat>& diag,
        const std::vector<TMat>& plus,
        const std::vector<std::pair<int, int>>& pairs)
    {
        const int d = static_cast<int>(diag.size());
        std::vector<std::vector<TMat>> g(d, std::vector<TMat>(d));
        for (int k = 0; k < d; ++k) { g[k][k] = diag[k]; }
        for (size_t ip = 0; ip < pairs.size(); ++ip)
        {
            const int k = pairs[ip].first;
            const int l = pairs[ip].second;
            // $G_{kl}=\mathcal F[X_+]-\tfrac12(G_{kk}+G_{ll})$
            const TMat off = plus[ip] - (diag[k] + diag[l]) * 0.5;
            g[k][l] = off;
            g[l][k] = off;
        }
        return g;
    }

    /// @brief The outcome of the Jahn-Teller search: the displacement direction that lowers one
    ///        branch of the multiplet fastest, and the electronic state that follows it.
    struct JTDirection
    {
        std::vector<double> displacement;  ///< the $3N$ direction, normalized, laid out (atom, xyz)
        std::vector<double> mixing;        ///< $v$, that branch's $d$ coefficients in the multiplet basis
        double slope = 0.0;                ///< $\|q(v)\|$, the branch's steepest descent rate
        int restarts_agreeing = 0;         ///< how many starting points reached this same optimum
        int iterations = 0;                ///< iterations used by the winning start
    };

    /// @brief Find the Jahn-Teller direction: the displacement that splits the multiplet and lowers
    ///        one branch as fast as possible.
    ///
    /// The Jahn-Teller theorem says a degenerate electronic state of a non-linear molecule is
    /// unstable against some symmetry-lowering displacement. Finding it is a JOINT optimization over
    /// the displacement and the mixing inside the subspace -- the two are determined together, which
    /// is why it cannot be done one Cartesian axis at a time:
    ///
    ///     $\min_{\|u\|=1}\lambda_{\min}\big(\sum_a u_a G^{(a)}\big)$.
    ///
    /// Written that way it looks like a non-convex problem on the unit sphere in $3N$ dimensions. It
    /// is not: since $\lambda_{\min}(M)=\min_{\|v\|=1}v^\top Mv$, the two minimizations can be
    /// swapped, and the inner one over $u$ has a closed form,
    ///
    ///     $\min_{\|u\|=1}\ u\cdot q(v)=-\|q(v)\|$, where $q(v)_a=v^\top G^{(a)}v$,
    ///
    /// leaving
    ///
    ///     $\min_{\|u\|=1}\lambda_{\min}\big(M(u)\big)=-\max_{\|v\|=1}\|q(v)\|$,
    ///     with the optimal direction $u^*=-q(v^*)/\|q(v^*)\|$.
    ///
    /// So the search runs over the $d$-dimensional subspace, not over $3N$ coordinates, and it has a
    /// plain physical reading: $q(v)$ is the gradient of the mixed state
    /// $|v\rangle=\sum_kv_k|X_k\rangle$, so the Jahn-Teller direction is the steepest-descent
    /// direction of whichever state in the multiplet has the largest gradient.
    ///
    /// Solved by alternating exact minimization -- $u\leftarrow-q(v)/\|q(v)\|$, then $v\leftarrow$
    /// the $\lambda_{\min}$ eigenvector of $M(u)$ -- which decreases the joint objective
    /// monotonically. The objective is homogeneous of degree one and only piecewise smooth, so
    /// several deterministic starting points are tried and the best kept; `restarts_agreeing`
    /// reports how many landed on it, which is the practical signal that the optimum is global.
    ///
    /// This is first order only. It gives the DIRECTION; the actual distortion amplitude needs the
    /// harmonic term as well ($Q\approx-g/k$), and a norm that means anything physically should be
    /// taken in mass-weighted coordinates rather than plain Cartesian ones.
    ///
    /// @param gflat  the gradient matrix, flattened as `gflat[a * d * d + k * d + l]`, symmetric in
    ///               (k, l). The sign convention is the caller's: feeding FORCES
    ///               ($-\partial\Omega/\partial R$) makes `displacement` point downhill directly,
    ///               while feeding gradients makes it point uphill.
    /// @param ncoord $3N$
    /// @param d      the multiplet's dimension, must be >= 2
    JTDirection find_jt_direction(const std::vector<double>& gflat, const int ncoord, const int d);

    /// @brief Split a branch gradient into the part shared by the whole multiplet and the part that
    ///        actually breaks the degeneracy.
    ///
    /// $q(v)=\bar q+\big(q(v)-\bar q\big)$ with $\bar q_a=\operatorname{Tr}G^{(a)}/d$. The first term
    /// is the multiplet average: totally symmetric, common to every branch, and it only relaxes the
    /// geometry without splitting anything. The second is the Jahn-Teller part. The distinction
    /// matters when reading the result -- at a stationary point of the average surface, which is
    /// where an `lr_relax_degen_mode = average` relaxation ends up, the first term vanishes and the
    /// whole direction is Jahn-Teller.
    ///
    /// @return the symmetric part $\bar q$; `jt_part` receives $q(v)-\bar q$
    std::vector<double> split_symmetric_part(const std::vector<double>& gflat,
        const int ncoord,
        const int d,
        const std::vector<double>& mixing,
        std::vector<double>& jt_part);
}
