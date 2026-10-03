#include "grad_matrix_degenerate.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

namespace LR
{
    std::vector<std::vector<int>> group_degenerate_states(const std::vector<double>& omega,
        const double thr)
    {
        const int nst = static_cast<int>(omega.size());
        std::vector<std::vector<int>> groups;
        if (nst == 0) { return groups; }

        // Casida output happens to come out ascending, but nothing here should depend on that.
        std::vector<int> order(nst);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(),
            [&omega](const int a, const int b) { return omega[a] < omega[b]; });

        if (thr <= 0.0)
        {
            for (int i = 0; i < nst; ++i) { groups.push_back(std::vector<int>(1, order[i])); }
            return groups;
        }

        double anchor = omega[order[0]];
        groups.push_back(std::vector<int>(1, order[0]));
        for (int i = 1; i < nst; ++i)
        {
            const int ist = order[i];
            // measured against the group's first (lowest) member, so the spread inside a group is
            // bounded by `thr` however many members it collects
            if (omega[ist] - anchor < thr)
            {
                groups.back().push_back(ist);
            }
            else
            {
                groups.push_back(std::vector<int>(1, ist));
                anchor = omega[ist];
            }
        }
        for (std::vector<int>& g : groups) { std::sort(g.begin(), g.end()); }
        return groups;
    }

    std::vector<std::pair<int, int>> degenerate_pairs(const int d)
    {
        std::vector<std::pair<int, int>> pairs;
        if (d < 2) { return pairs; }
        pairs.reserve(static_cast<size_t>(d) * (d - 1) / 2);
        for (int k = 0; k < d; ++k)
        {
            for (int l = k + 1; l < d; ++l) { pairs.push_back(std::make_pair(k, l)); }
        }
        return pairs;
    }

    namespace
    {
        /// $q(v)_a=v^\top G^{(a)}v$: the gradient of the mixed state $\sum_kv_k|X_k\rangle$.
        std::vector<double> branch_gradient(const std::vector<double>& gflat,
            const int ncoord,
            const int d,
            const std::vector<double>& v)
        {
            std::vector<double> q(ncoord, 0.0);
            for (int a = 0; a < ncoord; ++a)
            {
                const double* const block = gflat.data() + static_cast<size_t>(a) * d * d;
                double s = 0.0;
                for (int k = 0; k < d; ++k)
                {
                    for (int l = 0; l < d; ++l) { s += v[k] * block[k * d + l] * v[l]; }
                }
                q[a] = s;
            }
            return q;
        }

        double norm2(const std::vector<double>& x)
        {
            double s = 0.0;
            for (size_t i = 0; i < x.size(); ++i) { s += x[i] * x[i]; }
            return std::sqrt(s);
        }

        /// $M(u)=\sum_a u_aG^{(a)}$, as a dense $d\times d$ row-major block.
        std::vector<double> contract_direction(const std::vector<double>& gflat,
            const int ncoord,
            const int d,
            const std::vector<double>& u)
        {
            std::vector<double> m(static_cast<size_t>(d) * d, 0.0);
            for (int a = 0; a < ncoord; ++a)
            {
                const double* const block = gflat.data() + static_cast<size_t>(a) * d * d;
                for (int i = 0; i < d * d; ++i) { m[i] += u[a] * block[i]; }
            }
            return m;
        }

        /// The eigenvector of the smallest eigenvalue of a small symmetric matrix, by Jacobi
        /// rotations. Written out rather than taken from LAPACK so that this file stays free of
        /// the parallel/linear-algebra layer and can be unit-tested on its own; $d$ is the
        /// dimension of an electronic multiplet, so 2 or 3 in practice and never large.
        std::vector<double> min_eigenvector(std::vector<double> m, const int d)
        {
            std::vector<double> ev(static_cast<size_t>(d) * d, 0.0);
            for (int i = 0; i < d; ++i) { ev[i * d + i] = 1.0; }
            for (int sweep = 0; sweep < 100; ++sweep)
            {
                double off = 0.0;
                for (int i = 0; i < d; ++i)
                {
                    for (int j = i + 1; j < d; ++j) { off += m[i * d + j] * m[i * d + j]; }
                }
                if (off < 1e-30) { break; }
                for (int i = 0; i < d; ++i)
                {
                    for (int j = i + 1; j < d; ++j)
                    {
                        const double aij = m[i * d + j];
                        if (std::abs(aij) < 1e-300) { continue; }
                        const double theta = 0.5 * (m[j * d + j] - m[i * d + i]) / aij;
                        const double t = (theta >= 0.0 ? 1.0 : -1.0)
                            / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                        const double c = 1.0 / std::sqrt(t * t + 1.0);
                        const double s = t * c;
                        for (int k = 0; k < d; ++k)
                        {
                            const double mik = m[i * d + k];
                            const double mjk = m[j * d + k];
                            m[i * d + k] = c * mik - s * mjk;
                            m[j * d + k] = s * mik + c * mjk;
                        }
                        for (int k = 0; k < d; ++k)
                        {
                            const double mki = m[k * d + i];
                            const double mkj = m[k * d + j];
                            m[k * d + i] = c * mki - s * mkj;
                            m[k * d + j] = s * mki + c * mkj;
                            const double eki = ev[k * d + i];
                            const double ekj = ev[k * d + j];
                            ev[k * d + i] = c * eki - s * ekj;
                            ev[k * d + j] = s * eki + c * ekj;
                        }
                    }
                }
            }
            int best = 0;
            for (int i = 1; i < d; ++i)
            {
                if (m[i * d + i] < m[best * d + best]) { best = i; }
            }
            std::vector<double> v(d, 0.0);
            for (int k = 0; k < d; ++k) { v[k] = ev[k * d + best]; }
            const double n = norm2(v);
            for (int k = 0; k < d; ++k) { v[k] /= n; }
            return v;
        }

        /// Deterministic starting points: the unit vectors, then the normalized all-ones-with-signs
        /// patterns. Deterministic because a relaxation must give the same answer twice.
        std::vector<std::vector<double>> jt_start_points(const int d)
        {
            std::vector<std::vector<double>> starts;
            for (int k = 0; k < d; ++k)
            {
                std::vector<double> v(d, 0.0);
                v[k] = 1.0;
                starts.push_back(v);
            }
            const int nsign = 1 << (d - 1);   // fix the first sign: v and -v give the same q(v)
            for (int mask = 0; mask < nsign; ++mask)
            {
                std::vector<double> v(d, 1.0 / std::sqrt(static_cast<double>(d)));
                for (int k = 1; k < d; ++k)
                {
                    if ((mask >> (k - 1)) & 1) { v[k] = -v[k]; }
                }
                starts.push_back(v);
            }
            return starts;
        }
    }

    JTDirection find_jt_direction(const std::vector<double>& gflat, const int ncoord, const int d)
    {
        assert(d >= 2);
        assert(gflat.size() == static_cast<size_t>(ncoord) * d * d);
        JTDirection best;
        const std::vector<std::vector<double>> starts = jt_start_points(d);
        for (size_t is = 0; is < starts.size(); ++is)
        {
            std::vector<double> v = starts[is];
            double obj = -1.0;
            int it = 0;
            std::vector<double> q;
            for (; it < 200; ++it)
            {
                q = branch_gradient(gflat, ncoord, d, v);
                const double nq = norm2(q);
                // A vanishing gradient means this branch is already stationary: there is no
                // direction to report from this start, so leave it to the others.
                if (nq < 1e-14) { break; }
                if (nq - obj < 1e-12 * std::max(1.0, nq)) { obj = nq; break; }
                obj = nq;
                // u minimizes u.q(v) at fixed v; v then minimizes v^T M(u) v at fixed u. Both are
                // exact, so the joint objective -\|q\| decreases monotonically.
                std::vector<double> u(ncoord);
                for (int a = 0; a < ncoord; ++a) { u[a] = -q[a] / nq; }
                v = min_eigenvector(contract_direction(gflat, ncoord, d, u), d);
            }
            if (obj <= 0.0) { continue; }
            if (obj > best.slope * (1.0 + 1e-9))
            {
                best.slope = obj;
                best.mixing = v;
                best.displacement.assign(ncoord, 0.0);
                for (int a = 0; a < ncoord; ++a) { best.displacement[a] = q[a] / obj; }
                best.iterations = it;
                best.restarts_agreeing = 1;
            }
            else if (obj > best.slope * (1.0 - 1e-6))
            {
                ++best.restarts_agreeing;
            }
        }
        return best;
    }

    std::vector<double> split_symmetric_part(const std::vector<double>& gflat,
        const int ncoord,
        const int d,
        const std::vector<double>& mixing,
        std::vector<double>& jt_part)
    {
        const std::vector<double> q = branch_gradient(gflat, ncoord, d, mixing);
        std::vector<double> sym(ncoord, 0.0);
        jt_part.assign(ncoord, 0.0);
        for (int a = 0; a < ncoord; ++a)
        {
            const double* const block = gflat.data() + static_cast<size_t>(a) * d * d;
            double tr = 0.0;
            for (int k = 0; k < d; ++k) { tr += block[k * d + k]; }
            sym[a] = tr / static_cast<double>(d);
            jt_part[a] = q[a] - sym[a];
        }
        return sym;
    }
}
