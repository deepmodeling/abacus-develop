#include "source_lcao/module_lr/lr_mpl.h"

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_lcao/module_lr/grad_degen.h"
#include "source_lcao/module_lr/gradient_output.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>

namespace LR
{

template <typename T>
MplLvc build_mpl(const ct::Tensor& xz, const std::vector<double>& omega, const std::vector<ModuleBase::matrix>& diag,
                 const int nloc, const int ispin, const std::vector<int>& states, const UnitCell& ucell,
                 const MplEval& eval, std::ofstream& log)
{
    const int d = static_cast<int>(states.size());
    assert(d >= 2);
    assert(diag.size() == states.size());
    assert(omega.size() == states.size());
    const std::vector<std::pair<int, int>> pairs = degenerate_pairs(d);

    // 1. the two preconditions of route (A2), reported rather than enforced: the combinations
    //    $X_\pm$ are normalized eigenvectors only if the members are orthonormal, and the gradient
    //    is a single quadratic form only if they share one $\Omega$. An accidental near-degeneracy
    //    passes the grouping threshold but fails the second, and its `omega_spread` says so.
    double max_ovlp_err = 0.0;
    for (int k = 0; k < d; ++k)
    {
        for (int l = k; l < d; ++l)
        {
            const T* const xk = xz.template data<T>() + static_cast<size_t>(k) * nloc;
            const T* const xl = xz.template data<T>() + static_cast<size_t>(l) * nloc;
            T loc = static_cast<T>(0);
            for (int i = 0; i < nloc; ++i)
            {
                loc += LR_Util::get_conj(xk[i]) * xl[i];
            }
            Parallel_Reduce::reduce_all(loc);
            const double ref = (k == l) ? 1.0 : 0.0;
            max_ovlp_err = std::max(max_ovlp_err, std::abs(loc - ref));
        }
    }
    MplLvc mpl;
    mpl.ispin = ispin;
    mpl.states = states;
    const double omega_max = *std::max_element(omega.begin(), omega.end());
    const double omega_min = *std::min_element(omega.begin(), omega.end());
    mpl.omega_spread = omega_max - omega_min;
    // one $\Omega$ for every combination: the members share it up to `omega_spread`, and the mean
    // is the neutral choice for a vector that belongs to no single member
    double omega_sum = 0.0;
    for (const double e : omega)
    {
        omega_sum += e;
    }
    mpl.omega0 = omega_sum / static_cast<double>(d);
    log << " Degenerate multiplet of " << d << " states, Omega = " << mpl.omega0 << " Ry, spread = "
        << mpl.omega_spread << " Ry, max |<X_k|X_l> - delta_kl| = " << max_ovlp_err << std::endl;
    // Dimensionless diagnostic tolerance for the orthonormality needed by X_plus.
    // This warning does not change the multiplet or the computed gradient.
    const double orthonormality_warn_tol = 1e-6;
    if (max_ovlp_err > orthonormality_warn_tol)
    {
        log << " WARNING: this multiplet's eigenvectors are not orthonormal to"
            << " tolerance " << orthonormality_warn_tol
            << " (dimensionless overlap error); (X_k+X_l)/sqrt(2) may not be normalized"
            << " and the assembled gradient matrix may consequently be inaccurate." << std::endl;
    }

    // 2. the $d(d-1)/2$ combinations $X_+=(X_k+X_l)/\sqrt2$, all in one Z-vector solve
    const int npair = static_cast<int>(pairs.size());
    ct::Tensor Xp = LR_Util::newTensor<T>({npair, nloc});
    Xp.zero();
    for (int ip = 0; ip < npair; ++ip)
    {
        LR::combine_normalized(xz.template data<T>() + static_cast<size_t>(pairs[ip].first) * nloc,
                               xz.template data<T>() + static_cast<size_t>(pairs[ip].second) * nloc,
                               static_cast<size_t>(nloc), Xp.template data<T>() + static_cast<size_t>(ip) * nloc);
    }
    const std::vector<double> omega_pair(npair, mpl.omega0);
    const std::vector<ModuleBase::matrix> plus = eval(Xp, omega_pair, states[0]);

    // 3. $G_{kl}=\mathcal F[X_+]-\tfrac12(G_{kk}+G_{ll})$
    mpl.g = assemble_grad_matrix(diag, plus, pairs);
    print_grad_matrix(mpl.g, states, ucell, log);
    return mpl;
}

JtForce eval_jt_force(const std::vector<std::vector<ModuleBase::matrix>>& g, std::ostream& log)
{
    const int d = static_cast<int>(g.size());
    const int nat = g[0][0].nr;
    const int ncoord = nat * 3;
    // flatten to the layout `find_jt_direction` takes: one d x d block per nuclear coordinate.
    // FORCES go in, not gradients, so the direction that comes back points downhill.
    std::vector<double> gflat(static_cast<size_t>(ncoord) * d * d, 0.0);
    for (int iat = 0; iat < nat; ++iat)
    {
        for (int ixyz = 0; ixyz < 3; ++ixyz)
        {
            const int a = iat * 3 + ixyz;
            for (int k = 0; k < d; ++k)
            {
                for (int l = 0; l < d; ++l)
                {
                    gflat[static_cast<size_t>(a) * d * d + k * d + l] = g[k][l](iat, ixyz);
                }
            }
        }
    }
    const JTDirection jt = find_jt_direction(gflat, ncoord, d);
    std::vector<double> jt_part;
    const std::vector<double> sym = split_symmetric_part(gflat, ncoord, d, jt.mixing, jt_part);

    const double fac = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    log << " lr_degen_mode=jt: descending the steepest branch of the multiplet." << std::endl
        << "   |F| of that branch  = " << jt.slope * fac << " eV/Angstrom" << std::endl
        << "   mixing v            =";
    for (int k = 0; k < d; ++k)
    {
        log << " " << jt.mixing[k];
    }
    log << std::endl
        << "   starts agreeing     = " << jt.restarts_agreeing << " of " << (d + (1 << (d - 1))) << " (iterations "
        << jt.iterations << ")" << std::endl;
    // The two halves matter separately: the symmetric part is common to the whole multiplet and
    // only relaxes the geometry, while the remainder is what actually breaks the degeneracy. At a
    // stationary point of the average surface the first is zero and the whole force is Jahn-Teller.
    double nsym = 0.0;
    double njt = 0.0;
    for (int a = 0; a < ncoord; ++a)
    {
        nsym += sym[a] * sym[a];
        njt += jt_part[a] * jt_part[a];
    }
    log << "   |symmetric part|    = " << std::sqrt(nsym) * fac
        << " eV/Angstrom (common to the"
           " multiplet; relaxes the geometry without splitting it)"
        << std::endl
        << "   |Jahn-Teller part|  = " << std::sqrt(njt) * fac
        << " eV/Angstrom (the symmetry-"
           "breaking remainder)"
        << std::endl;
    if (std::sqrt(njt) < 1e-8)
    {
        log << " WARNING: the symmetry-breaking part vanishes -- every branch of this multiplet has"
               " the same gradient, so there is no Jahn-Teller direction to descend here. This is the"
               " expected outcome for a linear molecule, where the effect is second order"
               " (Renner-Teller) and no first-order term exists."
            << std::endl;
    }
    log << " NOTE: this is the first-order DIRECTION only. The distortion amplitude also needs the"
           " harmonic term, and the step norm is Cartesian rather than mass-weighted."
        << std::endl;

    // The force handed back is that of the descending branch, q(v) with v the optimal mixing --
    // i.e. the branch's own gradient, which is what a relaxation must follow. Once a step has
    // split the multiplet, `resolve_group` finds no group and the ordinary single-state path
    // takes over, so this mode is self-limiting by construction.
    JtForce out;
    out.force.create(nat, 3);
    out.mixing = jt.mixing;
    for (int iat = 0; iat < nat; ++iat)
    {
        for (int ixyz = 0; ixyz < 3; ++ixyz)
        {
            const int a = iat * 3 + ixyz;
            out.force(iat, ixyz) = sym[a] + jt_part[a];
        }
    }
    return out;
}

template MplLvc build_mpl<double>(const ct::Tensor&, const std::vector<double>&, const std::vector<ModuleBase::matrix>&,
                                  int, int, const std::vector<int>&, const UnitCell&, const MplEval&, std::ofstream&);
template MplLvc build_mpl<std::complex<double>>(const ct::Tensor&, const std::vector<double>&,
                                                const std::vector<ModuleBase::matrix>&, int, int,
                                                const std::vector<int>&, const UnitCell&, const MplEval&,
                                                std::ofstream&);

} // namespace LR
