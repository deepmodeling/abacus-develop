#include "dfpt_nscf.h"

#include "source_base/complexmatrix.h"
#include "source_base/global_function.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <random>

namespace ModuleDFPT
{

namespace
{

// column-major hermitian product V^dagger W -> m x m
void hermitian_project(
    const std::vector<std::vector<std::complex<double>>> &v,
    const std::vector<std::vector<std::complex<double>>> &w,
    int n,
    int m,
    ModuleBase::ComplexMatrix &out)
{
    out.zero_out();
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < m; ++i)
        {
            std::complex<double> acc(0.0, 0.0);
            const std::vector<std::complex<double>> &vi = v[i];
            const std::vector<std::complex<double>> &wj = w[j];
            for (int l = 0; l < n; ++l)
            {
                acc += std::conj(vi[l]) * wj[l];
            }
            out(i, j) = acc;
        }
    }
}

// M <- V * S for column-combination matrix S (n x m out from n x m_in).
// The two LapackConnector::zheev transposes cancel: the round-trip leaves
// the C++ matrix holding the eigenvector entries un-conjugated, so S is
// consumed exactly as returned.
void combine(
    const std::vector<std::vector<std::complex<double>>> &v,
    const ModuleBase::ComplexMatrix &s,
    int n,
    int m_in,
    int m_out,
    std::vector<std::vector<std::complex<double>>> &out)
{
    out.assign(m_out, std::vector<std::complex<double>>(n, std::complex<double>(0.0, 0.0)));
    for (int j = 0; j < m_out; ++j)
    {
        std::vector<std::complex<double>> &oj = out[j];
        for (int k = 0; k < m_in; ++k)
        {
            const std::complex<double> sk = s(k, j);
            if (sk == std::complex<double>(0.0, 0.0))
            {
                continue;
            }
            const std::vector<std::complex<double>> &vk = v[k];
            for (int l = 0; l < n; ++l)
            {
                oj[l] += sk * vk[l];
            }
        }
    }
}

// in-place modified Gram-Schmidt orthonormalization of the head columns of
// vecs against the orthonormal basis; returns false when a column is
// linearly dependent (norm below the floor) after two passes
bool orthogonalize(
    std::vector<std::vector<std::complex<double>>> &vecs,
    int head,
    int against_from,
    int against_to,
    int n)
{
    const double dep_floor = 1.0e-8; ///< empirical parameter: Gram-Schmidt dependence floor
    for (int pass = 0; pass < 2; ++pass)
    {
        for (int k = against_from; k < against_to; ++k)
        {
            std::complex<double> dot(0.0, 0.0);
            const std::vector<std::complex<double>> &vk = vecs[k];
            for (int l = 0; l < n; ++l)
            {
                dot += std::conj(vk[l]) * vecs[head][l];
            }
            for (int l = 0; l < n; ++l)
            {
                vecs[head][l] -= dot * vk[l];
            }
        }
    }
    double nrm = 0.0;
    for (int l = 0; l < n; ++l)
    {
        nrm += std::norm(vecs[head][l]);
    }
    nrm = std::sqrt(nrm);
    if (nrm < dep_floor * std::sqrt(static_cast<double>(n)))
    {
        return false;
    }
    const double inv = 1.0 / nrm;
    for (int l = 0; l < n; ++l)
    {
        vecs[head][l] *= inv;
    }
    return true;
}

} // namespace

int DFPT_NSCF::lowest_eigpairs(
    const DFPT_Stern::LinearOperator &hop,
    const std::vector<double> &precond,
    int nwant,
    int max_iter,
    double eig_tol,
    std::vector<std::vector<std::complex<double>>> &psi,
    std::vector<double> &eig)
{
    ModuleBase::TITLE("DFPT_NSCF", "lowest_eigpairs");
    ModuleBase::timer::start("DFPT_NSCF", "lowest_eigpairs");
    const int n = hop.dimension();
    const bool dbg = (std::getenv("DFPT_NSCF_DEBUG") != nullptr);
    psi.clear();
    eig.clear();
    if (n <= 0 || nwant <= 0 || nwant > n || static_cast<int>(precond.size()) != n || max_iter <= 0)
    {
        ModuleBase::timer::end("DFPT_NSCF", "lowest_eigpairs");
        return -1;
    }

    // subspace cap: a generous multiple of nwant bounds the zheev cost;
    // a hard restart (compressing to the current Ritz block) discards the
    // accumulated search directions and stalls the convergence short of
    // tight thresholds, so the cap is large enough that most problems
    // converge inside one window
    const int subspace_cap = nwant * 8 + 32;
    const int m_cap = std::min(n, subspace_cap);

    // initial guess: deterministic pseudo-random vectors filtered by the
    // kinetic preconditioner (high-|G| damping), then orthonormalized;
    // the fixed seed keeps runs reproducible for the reference tests
    std::vector<std::vector<std::complex<double>>> basis;
    basis.reserve(m_cap);
    std::mt19937 gen(20261125u);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    for (int m = 0; m < nwant; ++m)
    {
        std::vector<std::complex<double>> v(n);
        for (int l = 0; l < n; ++l)
        {
            const double re = uni(gen);
            const double im = uni(gen);
            const double pc = precond[l] > 1.0 ? std::sqrt(precond[l]) : 1.0;
            v[l] = std::complex<double>(re / pc, im / pc);
        }
        basis.push_back(std::move(v));
        if (!orthogonalize(basis, m, 0, m, n))
        {
            // degenerate initial draw (not expected with random data);
            // redraw with a fresh batch of randoms instead
            for (int l = 0; l < n; ++l)
            {
                const double pc = precond[l] > 1.0 ? std::sqrt(precond[l]) : 1.0;
                basis[m][l] = std::complex<double>(uni(gen) / pc, uni(gen) / pc);
            }
            orthogonalize(basis, m, 0, m, n);
        }
    }

    std::vector<std::vector<std::complex<double>>> wspace;
    int used = 0;
    bool converged = false;
    // progressive tightening: the compressed-restart limit cycle can pin
    // the residuals just above a one-shot target, so the tolerance starts
    // loose and steps down once met; every intermediate level is reachable
    // inside the window cycle
    double cur_tol = 1.0e-3; ///< empirical parameter: initial residual target
    for (int outer = 0; outer < max_iter && !converged; ++outer)
    {
        const int m = static_cast<int>(basis.size());
        used = outer + 1;
        // W = H V
        wspace.assign(m, std::vector<std::complex<double>>());
        for (int j = 0; j < m; ++j)
        {
            wspace[j].assign(n, std::complex<double>(0.0, 0.0));
            hop.apply(basis[j].data(), wspace[j].data());
        }
        // small hermitian eigenproblem M = V^dagger W
        ModuleBase::ComplexMatrix mat(m, m, true);
        hermitian_project(basis, wspace, n, m, mat);
        std::vector<double> theta(m, 0.0);
        std::vector<std::complex<double>> work(2 * m + 10);
        std::vector<double> rwork(3 * m + 10);
        int info = 0;
        LapackConnector::zheev(
            'V', 'U', m, mat, m, theta.data(), work.data(), static_cast<int>(work.size()), rwork.data(), &info);
        if (info != 0)
        {
            ModuleBase::timer::end("DFPT_NSCF", "lowest_eigpairs");
            return -2;
        }
        if (dbg)
        {
            std::cout << "NSCF outer=" << outer << " m=" << m << " theta:";
            for (int j = 0; j < std::min(m, nwant + 2); ++j)
            {
                std::cout << " " << theta[j];
            }
            std::cout << std::endl;
            double worst = 0.0;
            for (int a = 0; a < m; ++a)
            {
                for (int b = a + 1; b < m; ++b)
                {
                    std::complex<double> dot(0.0, 0.0);
                    for (int l = 0; l < n; ++l)
                    {
                        dot += std::conj(basis[a][l]) * basis[b][l];
                    }
                    worst = std::max(worst, std::abs(dot));
                }
            }
            std::cout << "NSCF basis worst off-diag overlap = " << worst << std::endl;
        }
        // Ritz vectors of the wanted block
        std::vector<std::vector<std::complex<double>>> ritz;
        combine(basis, mat, n, m, nwant, ritz);
        // residuals res_j = W s_j - theta_j * psi_j (see the combine()
        // note: the round-trip needs no conjugation)
        std::vector<std::vector<std::complex<double>>> res(nwant);
        converged = true;
        for (int j = 0; j < nwant; ++j)
        {
            res[j].assign(n, std::complex<double>(0.0, 0.0));
            for (int k = 0; k < m; ++k)
            {
                const std::complex<double> sk = mat(k, j);
                const std::vector<std::complex<double>> &wk = wspace[k];
                for (int l = 0; l < n; ++l)
                {
                    res[j][l] += sk * wk[l];
                }
            }
            for (int l = 0; l < n; ++l)
            {
                res[j][l] -= theta[j] * ritz[j][l];
            }
            double rn = 0.0;
            for (int l = 0; l < n; ++l)
            {
                rn += std::norm(res[j][l]);
            }
            rn = std::sqrt(rn);
            const double scale = std::max(1.0, std::fabs(theta[j]));
            if (rn > cur_tol * scale)
            {
                converged = false;
            }
        }
        if (converged)
        {
            if (cur_tol <= eig_tol)
            {
                break;
            }
            // met the intermediate level: step toward the final target
            // and keep iterating
            const double tightened_tol = cur_tol / 10.0;
            cur_tol = std::max(eig_tol, tightened_tol);
            converged = false;
        }
        // expansion: preconditioned residuals, orthogonalized against the
        // basis and each other; compress when at cap -- keeping a few
        // suboptimal Ritz directions beyond the wanted block (their
        // conjugate search history covers the near-threshold bands) and
        // re-orthonormalizing the whole compressed basis: zheev in a
        // numerically degenerate subspace returns Ritz vectors that are
        // only approximately orthogonal, and the contamination compounds
        // through the project matrices of later rounds
        if (m + nwant > m_cap)
        {
            const int keep_count = nwant + 4;
            const int keep = std::min(m, keep_count);
            std::vector<std::vector<std::complex<double>>> keep_basis;
            combine(basis, mat, n, m, keep, keep_basis);
            for (int k = 1; k < keep; ++k)
            {
                orthogonalize(keep_basis, k, 0, k, n);
            }
            // the compressed basis may have lost the faint components of
            // symmetry blocks that the initial guess only touched through
            // their random overlap; once those are purged, every future
            // residual lives in the retained block and the frozen spectrum
            // can never recover. Re-seeding keeps the search global.
            std::vector<std::complex<double>> t(n);
            for (int l = 0; l < n; ++l)
            {
                t[l] = std::complex<double>(uni(gen), uni(gen));
            }
            keep_basis.push_back(std::move(t));
            if (keep < n && orthogonalize(keep_basis, keep, 0, keep, n))
            {
                basis = std::move(keep_basis);
            }
            else
            {
                // the reseed has no independent component (keep == n, the
                // compressed basis already spans the full space) or it is
                // linearly dependent: keeping either case would leave a
                // rank-deficient basis with a spurious zero Ritz pair
                keep_basis.pop_back();
                basis = std::move(keep_basis);
            }
        }
        int kept_count = 0;
        for (int j = 0; j < nwant; ++j)
        {
            std::vector<std::complex<double>> t(n);
            for (int l = 0; l < n; ++l)
            {
                const double denom = precond[l] - theta[j];
                const double floor_den = 0.1; ///< empirical parameter: denominator clamp
                const double use = std::fabs(denom) < floor_den ? (denom >= 0.0 ? floor_den : -floor_den) : denom;
                t[l] = res[j][l] / use;
            }
            if (static_cast<int>(basis.size()) >= n)
            {
                break; // the subspace already spans the whole space
            }
            const int base = static_cast<int>(basis.size());
            basis.push_back(std::move(t));
            const bool kept = orthogonalize(basis, base, 0, base, n);
            if (dbg)
            {
                std::cout << "NSCF   kept[" << j << "]=" << kept << std::endl;
            }
            if (!kept)
            {
                // dependent direction: drop it (the remaining residuals
                // still expand the subspace)
                basis.pop_back();
            }
            else
            {
                ++kept_count;
            }
        }
        if (kept_count == 0 && static_cast<int>(basis.size()) < n)
        {
            // every expansion direction lived inside the span (possible
            // when the preconditioner matches the operator structure);
            // inject a random orthogonal direction so the iteration can
            // still make progress instead of stagnating
            std::vector<std::complex<double>> t(n);
            for (int l = 0; l < n; ++l)
            {
                t[l] = std::complex<double>(uni(gen), uni(gen));
            }
            const int base = static_cast<int>(basis.size());
            basis.push_back(std::move(t));
            if (!orthogonalize(basis, base, 0, base, n))
            {
                basis.pop_back();
            }
        }
    }

    // final Rayleigh-Ritz on the surviving subspace
    {
        const int m = static_cast<int>(basis.size());
        wspace.assign(m, std::vector<std::complex<double>>());
        for (int j = 0; j < m; ++j)
        {
            wspace[j].assign(n, std::complex<double>(0.0, 0.0));
            hop.apply(basis[j].data(), wspace[j].data());
        }
        ModuleBase::ComplexMatrix mat(m, m, true);
        hermitian_project(basis, wspace, n, m, mat);
        std::vector<double> theta(m, 0.0);
        std::vector<std::complex<double>> work(2 * m + 10);
        std::vector<double> rwork(3 * m + 10);
        int info = 0;
        LapackConnector::zheev(
            'V', 'U', m, mat, m, theta.data(), work.data(), static_cast<int>(work.size()), rwork.data(), &info);
        if (info != 0)
        {
            ModuleBase::timer::end("DFPT_NSCF", "lowest_eigpairs");
            return -2;
        }
        psi.assign(nwant, std::vector<std::complex<double>>());
        combine(basis, mat, n, m, nwant, psi);
        eig.assign(theta.begin(), theta.begin() + nwant);
        if (dbg)
        {
            for (int j = 0; j < nwant; ++j)
            {
                double nrm = 0.0;
                std::vector<std::complex<double>> hp(n, std::complex<double>(0.0, 0.0));
                hop.apply(psi[j].data(), hp.data());
                std::complex<double> hv(0.0, 0.0);
                for (int l = 0; l < n; ++l)
                {
                    nrm += std::norm(psi[j][l]);
                    hv += std::conj(psi[j][l]) * hp[l];
                }
                std::cout << "NSCF final psi[" << j << "] norm=" << std::sqrt(nrm) << " <psi|H|psi>=" << hv.real()
                          << " theta=" << theta[j] << std::endl;
            }
        }
    }
    ModuleBase::timer::end("DFPT_NSCF", "lowest_eigpairs");
    return converged ? 0 : used;
}

} // namespace ModuleDFPT
