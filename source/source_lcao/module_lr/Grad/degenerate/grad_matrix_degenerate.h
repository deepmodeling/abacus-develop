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
/// Derivation, the self-checks it admits, and why the alternative (an explicitly bilinear Z-vector
/// right-hand side) is the more expensive route:
/// `LR-Grad-formulas/2026-09-简并激发态梯度-实测和讨论.md` section 5.4.
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
}
