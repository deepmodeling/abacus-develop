// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_stern.h"

#include <cmath>

namespace ModuleDFPT {

DFPT_Stern::DFPT_Stern() {}

DFPT_Stern::~DFPT_Stern() {}

namespace {

double real_vdot(const std::vector<std::complex<double>>& a,
                 const std::vector<std::complex<double>>& b)
{
    // Re <a|b> = Re sum_i conj(a_i) b_i (the CG scalar products of a
    // Hermitian operator are real up to roundoff)
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        s += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    }
    return s;
}

} // namespace

void DFPT_Stern::apply_pv(const std::vector<std::vector<std::complex<double>>>& occ_kq,
                          const std::vector<std::complex<double>>& x,
                          std::vector<std::complex<double>>& px) const
{
    px = x;
    // two modified Gram-Schmidt sweeps keep the complement exact enough for
    // long CG chains even when the occupied set is only machine-orthonormal;
    // each projection collects its coefficient before subtracting, so px may
    // alias x
    for (int sweep = 0; sweep < 2; ++sweep)
    {
        for (size_t m = 0; m < occ_kq.size(); ++m)
        {
            const std::vector<std::complex<double>>& u = occ_kq[m];
            std::complex<double> c(0.0, 0.0);
            for (size_t i = 0; i < u.size(); ++i)
            {
                c += std::conj(u[i]) * px[i];
            }
            for (size_t i = 0; i < u.size(); ++i)
            {
                px[i] -= c * u[i];
            }
        }
    }
}

int DFPT_Stern::solve(const LinearOperator& aop,
                      const std::vector<std::vector<std::complex<double>>>& occ_kq,
                      const std::vector<std::complex<double>>& b,
                      int max_iter,
                      double conv_thr,
                      std::vector<std::complex<double>>& dpsi,
                      double& residual) const
{
    const int n = aop.dimension();
    dpsi.assign(n, std::complex<double>(0.0, 0.0));
    if (n == 0 || static_cast<int>(b.size()) != n || max_iter <= 0)
    {
        residual = 0.0;
        return 0;
    }
    for (size_t m = 0; m < occ_kq.size(); ++m)
    {
        if (static_cast<int>(occ_kq[m].size()) != n)
        {
            residual = 0.0;
            return 0;
        }
    }

    std::vector<std::complex<double>> pb(n);
    apply_pv(occ_kq, b, pb);
    const double bnorm = std::sqrt(real_vdot(pb, pb));
    if (bnorm < 1.0e-300)
    {
        // the right-hand side lies inside the occupied subspace: the
        // projected system is homogeneous and dpsi = 0 solves it exactly
        residual = 0.0;
        return 0;
    }

    std::vector<std::complex<double>> r = pb;
    std::vector<std::complex<double>> p = pb;
    std::vector<std::complex<double>> ap(n);
    std::vector<std::complex<double>> pap(n);
    std::vector<std::complex<double>> tmp(n);
    double rnorm2 = real_vdot(r, r);
    int used = 0;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        used = iter + 1;
        aop.apply(p.data(), ap.data());
        apply_pv(occ_kq, ap, pap);
        double pAp = real_vdot(p, pap);
        if (pAp <= 0.0)
        {
            // loss of positive definiteness along p (roundoff drift out of
            // the complement): restart the direction from the residual
            p = r;
            aop.apply(p.data(), ap.data());
            apply_pv(occ_kq, ap, pap);
            pAp = real_vdot(p, pap);
            if (pAp <= 0.0)
            {
                rnorm2 = real_vdot(r, r);
                break;
            }
        }
        const double alpha = rnorm2 / pAp;
        for (int i = 0; i < n; ++i)
        {
            dpsi[i] += alpha * p[i];
            r[i] -= alpha * pap[i];
        }
        const double rnew2 = real_vdot(r, r);
        if (std::sqrt(rnew2) / bnorm < conv_thr)
        {
            rnorm2 = rnew2;
            break;
        }
        const double beta = rnew2 / rnorm2;
        for (int i = 0; i < n; ++i)
        {
            p[i] = r[i] + beta * p[i];
        }
        apply_pv(occ_kq, p, p); // in-place re-projection of the search direction
        rnorm2 = rnew2;
    }
    // final hygiene: remove any occupied-subspace leakage of the solution
    apply_pv(occ_kq, dpsi, tmp);
    dpsi.swap(tmp);
    residual = std::sqrt(rnorm2) / bnorm;
    return used;
}

} // namespace ModuleDFPT
