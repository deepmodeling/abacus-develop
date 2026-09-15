// KB-projector construction of DFPT_Pert, split out of dfpt_pert.cpp:
// the radial vq integral, the real spherical harmonics (l <= 2) with
// their gradients, and the vkb / dvkb builders on the (k+q) basis. All
// formulas are moved verbatim from the original body; the per-l
// spherical-harmonic channels are factored into file-local helpers.

#include "dfpt_pert.h"

#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/math_sphbes.h"
#include "source_base/tool_quit.h"
#include "source_cell/atom_pseudo.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace ModuleDFPT
{

namespace
{

/// Y_{1,m} channel of the orthonormal real spherical harmonics (m in [-1, 1])
double ylm_l1(int m, double nx, double ny, double nz)
{
    switch (m)
    {
    case -1:
        return -0.5 * std::sqrt(3.0 / ModuleBase::PI) * ny;
    case 0:
        return 0.5 * std::sqrt(3.0 / ModuleBase::PI) * nz;
    case 1:
        return -0.5 * std::sqrt(3.0 / ModuleBase::PI) * nx;
    }
    return 0.0;
}

/// Y_{2,m} channel of the orthonormal real spherical harmonics (m in [-2, 2])
double ylm_l2(int m, double nx, double ny, double nz)
{
    switch (m)
    {
    case -2:
        return 0.5 * std::sqrt(15.0 / ModuleBase::PI) * nx * ny;
    case -1:
        return -0.5 * std::sqrt(15.0 / ModuleBase::PI) * nz * ny;
    case 0:
        return 0.25 * std::sqrt(5.0 / ModuleBase::PI) * (3.0 * nz * nz - 1.0);
    case 1:
        return -0.5 * std::sqrt(15.0 / ModuleBase::PI) * nz * nx;
    case 2:
        return 0.25 * std::sqrt(15.0 / ModuleBase::PI) * (nx * nx - ny * ny);
    }
    return 0.0;
}

/// gradient of the Y_{1,m} channel (m in [-1, 1])
void grad_l1(int m, double c1, double* grad)
{
    switch (m)
    {
    case -1:
        grad[1] = -c1;
        return;
    case 0:
        grad[2] = c1;
        return;
    case 1:
        grad[0] = -c1;
        return;
    }
}

/// gradient of the Y_{2,m} channel (m in [-2, 2])
void grad_l2(int m, double c2, double c20, double x, double y, double z, double* grad)
{
    switch (m)
    {
    case -2:
        grad[0] = c2 * y;
        grad[1] = c2 * x;
        return;
    case -1:
        grad[1] = -c2 * z;
        grad[2] = -c2 * y;
        return;
    case 0:
        grad[2] = 6.0 * c20 * z;
        return;
    case 1:
        grad[0] = -c2 * z;
        grad[2] = -c2 * x;
        return;
    case 2:
        grad[0] = 2.0 * c20 * x;
        grad[1] = -2.0 * c20 * y;
        return;
    }
}

} // namespace

double DFPT_Pert::radial_vq(int it, int ib, double g) const
{
    ModuleBase::TITLE("DFPT_Pert", "radial_vq");
    ModuleBase::timer::start("DFPT_Pert", "radial_vq");
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    const int l = ncpp.lll[ib];
    int kkbeta = ncpp.kkbeta;
    if (kkbeta > 0 && (kkbeta % 2 == 0))
    {
        --kkbeta;
    }
    std::vector<double> jl(kkbeta);
    std::vector<double> aux(kkbeta);
    ModuleBase::Sphbes::Spherical_Bessel(kkbeta, ncpp.r.data(), g, l, jl.data());
    for (int ir = 0; ir < kkbeta; ++ir)
    {
        aux[ir] = ncpp.betar(ib, ir) * jl[ir] * ncpp.r[ir];
    }
    double v = 0.0;
    ModuleBase::Integral::Simpson_Integral(kkbeta, aux.data(), ncpp.rab.data(), v);
    // tab convention from vnl_pw.cpp: (4pi/sqrt(Omega)) * integral
    ModuleBase::timer::end("DFPT_Pert", "radial_vq");
    return v * ModuleBase::FOUR_PI / std::sqrt(ucell_->omega);
}

double DFPT_Pert::real_ylm(int l, int m, const ModuleBase::Vector3<double>& ghat) const
{
    ModuleBase::TITLE("DFPT_Pert", "real_ylm");
    ModuleBase::timer::start("DFPT_Pert", "real_ylm");
    // orthonormal real spherical harmonics Y_{l,m} for l <= 2 with the
    // standard convention, m in [-l, l]:
    //   Y_{l,0}      = sqrt((2l+1)/4pi) P_l^0(cos0)
    //   Y_{l,m>0}    = sqrt(2 (2l+1)/4pi (l-m)!/(l+m)!) P_l^m(cos0) cos(m phi)
    //   Y_{l,m<0}    = sqrt(2 (2l+1)/4pi (l-|m|)!/(l+|m|)!) P_l^{|m|}(cos0) sin(|m| phi)
    // with the associated Legendre convention P_1^1 = -sin0, P_2^1 = -3 sin0 cos0,
    // P_2^2 = 3 sin^2 0. The ABACUS GS vkb applies an additional (-1)^|m| phase
    // for the m>0 channels; exact GS parity is reconciled in the diamond
    // end-to-end test (C7), while the C1 identity test is convention-independent.
    const double ghat_zero_tol = 1.0e-12; ///< empirical parameter: |ghat| floor for the zero-direction limit
    const double x = ghat.x;
    const double y = ghat.y;
    const double z = ghat.z;
    const double r = std::sqrt(x * x + y * y + z * z);
    if (r < ghat_zero_tol)
    {
        ModuleBase::timer::end("DFPT_Pert", "real_ylm");
        return (l == 0) ? 0.5 * std::sqrt(1.0 / ModuleBase::PI) : 0.0;
    }
    const double nx = x / r;
    const double ny = y / r;
    const double nz = z / r;
    double ylm = 0.0;
    switch (l)
    {
    case 0:
        ylm = 0.5 * std::sqrt(1.0 / ModuleBase::PI);
        break;
    case 1:
        ylm = ylm_l1(m, nx, ny, nz);
        break;
    case 2:
        ylm = ylm_l2(m, nx, ny, nz);
        break;
    default:
        ModuleBase::WARNING_QUIT("DFPT_Pert::real_ylm", "real_ylm implemented for l<=2 only (DFPT NC path).");
    }
    ModuleBase::timer::end("DFPT_Pert", "real_ylm");
    return ylm;
}

void DFPT_Pert::build_vkb(int it,
                          int ia,
                          const std::vector<ModuleBase::Vector3<double>>& gk,
                          std::vector<std::vector<std::complex<double>>>& vkb) const
{
    ModuleBase::TITLE("DFPT_Pert", "build_vkb");
    ModuleBase::timer::start("DFPT_Pert", "build_vkb");
    // per-type projector bookkeeping mirrors the ground-state vnl_pw.cpp layout:
    // every radial beta (nbeta) with angular momentum l spins out (2l+1)
    // projectors with combined index lm = l^2 + m, m in 0..2l (i.e. the real
    // harmonic m channels -l..l walked as m' = (-1)^(m+1) ceil... ABACUS ylm
    // block: m=0, +1, -1, +2, -2, ...). We use the signed m' directly.
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    const int nh = ncpp.nh;
    const int ngk = static_cast<int>(gk.size());
    const ModuleBase::Vector3<double>& tau = ucell_->atoms[it].tau[ia];
    vkb.assign(nh, std::vector<std::complex<double>>(ngk, std::complex<double>(0.0, 0.0)));
    if (nh == 0)
    {
        ModuleBase::timer::end("DFPT_Pert", "build_vkb");
        return;
    }
    int mu = 0;
    for (int ib = 0; ib < ncpp.nbeta; ++ib)
    {
        const int l = ncpp.lll[ib];
        if (l > 2)
        {
            ModuleBase::WARNING_QUIT("DFPT_Pert::build_vkb", "DFPT NC projector path implemented for l<=2 only.");
        }
        const std::complex<double> pref = std::pow(std::complex<double>(0.0, -1.0), l); // (-i)^l
        for (int m = 0; m < 2 * l + 1; ++m)
        {
            // ABACUS real-harmonic walk over the m channels of this radial beta:
            // m=0 -> m'=0; m=1 -> m'=+1; m=2 -> m'=-1; m=3 -> m'=+2; m=4 -> m'=-2
            const int mr = (m == 0) ? 0 : ((m % 2 == 1) ? (m + 1) / 2 : -(m / 2));
            for (int ig = 0; ig < ngk; ++ig)
            {
                const ModuleBase::Vector3<double>& G = gk[ig];         // k(+q)+G, 2*pi/lat0
                const double gnorm = std::sqrt(G * G) * ucell_->tpiba; // bohr^-1
                // real_ylm handles the |G|=0 point itself (Y_00 is
                // direction-independent; l>0 channels vanish there together
                // with vq), so the raw vector is passed directly.
                const double ylm = real_ylm(l, mr, G);
                const double vq = radial_vq(it, ib, gnorm);
                // GS structure-factor convention (stru_fac.cpp get_sk /
                // eigts, ci_tpi = -2pi i): exp(-i 2pi (gk.tau))
                const double arg = -ModuleBase::TWO_PI * (G * tau);
                const std::complex<double> phase(std::cos(arg), std::sin(arg));
                vkb[mu][ig] = pref * ylm * vq * phase;
            }
            ++mu;
        }
    }
    ModuleBase::timer::end("DFPT_Pert", "build_vkb");
}

void DFPT_Pert::grad_real_ylm(int l, int m, const ModuleBase::Vector3<double>& ghat, double grad[3]) const
{
    ModuleBase::TITLE("DFPT_Pert", "grad_real_ylm");
    ModuleBase::timer::start("DFPT_Pert", "grad_real_ylm");
    // analytic gradients of the real_ylm polynomials (l <= 2), consistent
    // with the conventions documented above real_ylm
    const double x = ghat.x;
    const double y = ghat.y;
    const double z = ghat.z;
    const double c1 = 0.5 * std::sqrt(3.0 / ModuleBase::PI);
    const double c2 = 0.5 * std::sqrt(15.0 / ModuleBase::PI);
    const double c20 = 0.25 * std::sqrt(5.0 / ModuleBase::PI);
    grad[0] = grad[1] = grad[2] = 0.0;
    switch (l)
    {
    case 0:
        break;
    case 1:
        grad_l1(m, c1, grad);
        break;
    case 2:
        grad_l2(m, c2, c20, x, y, z, grad);
        break;
    default:
        ModuleBase::WARNING_QUIT("DFPT_Pert::grad_real_ylm", "grad_real_ylm implemented for l<=2 only (DFPT NC path).");
    }
    ModuleBase::timer::end("DFPT_Pert", "grad_real_ylm");
}

void DFPT_Pert::build_vkb_dk(int it,
                             int ia,
                             int dir,
                             const std::vector<ModuleBase::Vector3<double>>& gk,
                             std::vector<std::vector<std::complex<double>>>& vkb,
                             std::vector<std::vector<std::complex<double>>>& dvkb) const
{
    ModuleBase::TITLE("DFPT_Pert", "build_vkb_dk");
    ModuleBase::timer::start("DFPT_Pert", "build_vkb_dk");
    const pseudo& ncpp = ucell_->atoms[it].ncpp;
    const int nh = ncpp.nh;
    const int ngk = static_cast<int>(gk.size());
    const ModuleBase::Vector3<double>& tau = ucell_->atoms[it].tau[ia];
    if (static_cast<int>(vkb.size()) != nh || static_cast<int>(vkb[0].size()) != ngk)
    {
        ModuleBase::WARNING_QUIT("DFPT_Pert::build_vkb_dk", "vkb must be built on the same gk list first.");
    }
    dvkb.assign(nh, std::vector<std::complex<double>>(ngk, std::complex<double>(0.0, 0.0)));
    if (nh == 0)
    {
        ModuleBase::timer::end("DFPT_Pert", "build_vkb_dk");
        return;
    }
    const double dg = 1.0e-4; // bohr^-1, radial central-difference step
    const double gmag_zero_tol = 1.0e-10; ///< empirical parameter: |G| floor (2*pi/lat0) for the angular terms
    int mu = 0;
    for (int ib = 0; ib < ncpp.nbeta; ++ib)
    {
        const int l = ncpp.lll[ib];
        const std::complex<double> pref = std::pow(std::complex<double>(0.0, -1.0), l); // (-i)^l
        for (int m = 0; m < 2 * l + 1; ++m)
        {
            const int mr = (m == 0) ? 0 : ((m % 2 == 1) ? (m + 1) / 2 : -(m / 2));
            for (int ig = 0; ig < ngk; ++ig)
            {
                const ModuleBase::Vector3<double>& G = gk[ig];
                const double gmag = std::sqrt(G * G);      // 2*pi/lat0 units
                const double gnorm = gmag * ucell_->tpiba; // bohr^-1
                const double vq0 = radial_vq(it, ib, gnorm);
                const double dvq = (radial_vq(it, ib, gnorm + dg) - radial_vq(it, ib, std::max(0.0, gnorm - dg)))
                                   / (dg * (gnorm > dg ? 2.0 : 1.0));
                const double arg = -ModuleBase::TWO_PI * (G * tau);
                const std::complex<double> phase(std::cos(arg), std::sin(arg));
                const std::complex<double> dphase = std::complex<double>(0.0, -ModuleBase::TWO_PI * tau[dir]) * phase;
                double dy[3] = {0.0, 0.0, 0.0};
                double ylm = 0.0;
                if (gmag > gmag_zero_tol)
                {
                    const ModuleBase::Vector3<double> ghat = G * (1.0 / gmag);
                    ylm = real_ylm(l, mr, ghat);
                    grad_real_ylm(l, mr, ghat, dy);
                    const double gdir[3] = {ghat.x, ghat.y, ghat.z};
                    // chain rule dghat/dk_dir = (e_dir - ghat*ghat_dir)/|G|
                    double dylm_dir = 0.0;
                    for (int c = 0; c < 3; ++c)
                    {
                        dylm_dir += dy[c] * ((c == dir ? 1.0 : 0.0) - gdir[c] * gdir[dir]);
                    }
                    dylm_dir /= gmag;
                    // radial chain: dg/dk_dir = tpiba * ghat_dir
                    const double dradial = dvq * ucell_->tpiba * gdir[dir];
                    dvkb[mu][ig] = pref * phase * (dylm_dir * vq0 + ylm * dradial) + pref * ylm * vq0 * dphase;
                }
                else
                {
                    // degenerate |G| = 0: only the l = 0 channel survives
                    // (real_ylm convention); keep only the phase term
                    ylm = (l == 0) ? 0.5 * std::sqrt(1.0 / ModuleBase::PI) : 0.0;
                    dvkb[mu][ig] = pref * ylm * vq0 * dphase;
                }
            }
            ++mu;
        }
    }
    ModuleBase::timer::end("DFPT_Pert", "build_vkb_dk");
}

} // namespace ModuleDFPT
