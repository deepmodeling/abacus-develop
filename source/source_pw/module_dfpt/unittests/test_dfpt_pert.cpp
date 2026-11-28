#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <functional>
#include <map>
#include <vector>

// serial unit test of the first-order perturbation potential (C1).
// Everything here runs without __MPI: the plane-wave bases are built through
// the real serial initgrids/initparameters/setuptransform path on a shared
// FFT grid, exactly like the production setup_pwrho/setup_pwwfc sequence.

#include "dfpt_test_tool.h"
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/qlist.h"
#include "source_cell/stru_fac.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#include "source_pw/module_dfpt/dfpt_pw_data.h"
#include "source_pw/module_pwdft/dftu_base.h"

// ctor/dtor stubs for the cell/spepot/stru_fac link closures live in the
// shared test/dfpt_test_mocks.cpp compiled into every DFPT test binary.

/************************************************
 *  serial unit test of DFPT_Pert (C1)
 ***********************************************/

/**
 * - Tested Functions:
 *   - rho_gvec: reconstruction of rho-grid G vectors must coincide with the
 *     distributed gcar array filled by the real serial collect_local_pw().
 *   - dVloc_dtau: the analytic first-order local potential is validated
 *     against a finite difference of the full displaced potential
 *     vloc(|Delta+q|) * exp(i 2pi (Delta+q).tau) for q != 0.
 *   - build_dv + apply_dv: the stored reciprocal coefficients and the FFT
 *     convolution dpsi(G'') = sum_G' psi(G') c(G''-G') are checked against
 *     the analytic matrix elements on the k+q basis (no extra q-phase).

 *   - build_vkb: l=0 projector against an independent hand-rolled Simpson
 *     transform; pure phase response under a tau shift.
 *   - dVnl_dtau: the two-term NC identity is validated against a finite
 *     difference of the displaced separable operator; USPP is rejected.
 *   - build_dv under with_u()/u_active()==false (pure-PW DFT+U safety).
 */

// per-file self-contained fixture (gtest style): owns its objects and
// calls the dfpt_test free-function utilities
class DFPTPertSerialTest : public ::testing::Test
{
  protected:
    ModuleBase::Matrix3 latvec_;
    UnitCell ucell_;
    ModulePW::PW_Basis pw_rho_;
    ModulePW::PW_Basis_K pw_wfc_;
    ModuleCell::QList qlist_;
    ModuleDFPT::DFPT_PW_Data data_;
    ModuleBase::Vector3<double> q_cart_;
    Structure_Factor sf_;
    ModuleDFPT::DFPT_Pert pert_;

    // default (k, q): q generic, k = -q so k+q = 0; tau of atom 0
    const ModuleBase::Vector3<double> k_d_{-0.13, 0.0, -0.07};
    const ModuleBase::Vector3<double> q_d_{0.13, 0.0, 0.07};
    const ModuleBase::Vector3<double> tau_{1.1, 2.3, 0.7};

    void SetUp() override
    {
        dfpt_test::setup_cell(ucell_, latvec_);
        dfpt_test::setup_bases(pw_rho_, pw_wfc_, qlist_, data_, ucell_, k_d_, q_d_, 2, q_cart_);
        pert_.init(ucell_, &pw_rho_, &pw_wfc_, sf_);
    }

    void TearDown() override
    {
        delete[] ucell_.atoms;
        ucell_.atoms = nullptr;
    }
};

TEST_F(DFPTPertSerialTest, RhoGvecMatchesDistributedGcar)
{
    ASSERT_GT(pw_rho_.npw, 0);
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        ModuleBase::Vector3<double> g;
        pert_.rho_gvec(ig, g);
        EXPECT_DOUBLE_EQ(g.x, pw_rho_.gcar[ig].x);
        EXPECT_DOUBLE_EQ(g.y, pw_rho_.gcar[ig].y);
        EXPECT_DOUBLE_EQ(g.z, pw_rho_.gcar[ig].z);
    }
}

TEST_F(DFPTPertSerialTest, DVlocDtauMatchesFiniteDifference)
{
    const ModuleBase::Vector3<double> qs[2] = {q_cart_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};
    const double eps = 1.0e-6; // lat0 units
    for (int iq = 0; iq < 2; ++iq)
    {
        for (int dir = 0; dir < 3; dir += 2)
        {
            std::vector<std::complex<double>> dv;
            pert_.dVloc_dtau(0, dir, qs[iq], dv);
            ASSERT_EQ(dv.size(), static_cast<size_t>(pw_rho_.npw));
            ModuleBase::Vector3<double> d(0.0, 0.0, 0.0);
            d[dir] = 1.0;
            for (int ig = 0; ig < pw_rho_.npw; ++ig)
            {
                const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + qs[iq];
                if (w * w < 1.0e-12)
                {
                    EXPECT_EQ(dv[ig], std::complex<double>(0.0, 0.0));
                    continue;
                }
                // finite difference of vloc(|Delta+q|) e^{-i 2pi (Delta+q).tau}
                // per bohr of displacement (GS stru_fac phase convention)
                const double ap = -ModuleBase::TWO_PI * (w * (tau_ + eps * d));
                const double am = -ModuleBase::TWO_PI * (w * (tau_ - eps * d));
                const std::complex<double> fd =
                    dfpt_test::vloc_coulomb(ucell_.atoms[0].ncpp.zv, ucell_.omega, (w * w) * ucell_.tpiba2) *
                    (std::polar(1.0, ap) - std::polar(1.0, am)) / (2.0 * eps * dfpt_test::lat0);
                EXPECT_NEAR(dv[ig].real(), fd.real(), 1.0e-9);
                EXPECT_NEAR(dv[ig].imag(), fd.imag(), 1.0e-9);
            }
        }
    }
}

TEST_F(DFPTPertSerialTest, ApplyDvConvolutionMatchesAnalyticMatrixElement)
{
    pert_.build_dv(0, 0, 0, data_); // atom 0 displaced along x, q = q_d_
    EXPECT_EQ(data_.get_pert_atom(), 0);
    EXPECT_EQ(data_.get_pert_dir(), 0);

    // stored reciprocal coefficients equal the analytic ones
    const std::vector<std::complex<double>> dv_recip = data_.get_dv_recip_c(0, 0);
    ASSERT_EQ(dv_recip.size(), static_cast<size_t>(pw_rho_.npw));
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        const std::complex<double> expect = dfpt_test::analytic_dvloc(ucell_, tau_, 0, pw_rho_.gcar[ig] + q_cart_);
        EXPECT_NEAR(dv_recip[ig].real(), expect.real(), 1.0e-9);
        EXPECT_NEAR(dv_recip[ig].imag(), expect.imag(), 1.0e-9);
    }

    // wavefunctions: band 0 = single plane wave at k, band 1 = two components
    const int npwk = pw_wfc_.npwk[0];
    const ModuleBase::Vector3<double> kc = pw_wfc_.kvec_c[0];
    int ig_zero = -1, ig_x = -1, ig_y = -1;
    for (int ig = 0; ig < npwk; ++ig)
    {
        const long long key = dfpt_test::gkey(pw_wfc_.getgpluskcar(0, ig) - kc, dfpt_test::cubic_edge);
        if (key == dfpt_test::fkey(0, 0, 0))
        {
            ig_zero = ig;
        }
        else if (key == dfpt_test::fkey(1, 0, 0))
        {
            ig_x = ig;
        }
        else if (key == dfpt_test::fkey(0, 1, 0))
        {
            ig_y = ig;
        }
    }
    ASSERT_GE(ig_zero, 0);
    ASSERT_GE(ig_x, 0);
    ASSERT_GE(ig_y, 0);

    psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
    psi.zero_out();
    psi(0, 0, ig_zero) = std::complex<double>(1.0, 0.0);
    psi(0, 1, ig_x) = std::complex<double>(0.7, 0.0);
    psi(0, 1, ig_y) = std::complex<double>(0.3, 0.2);

    pert_.apply_dv(0, 0, psi, data_);

    // expected: dpsi(G'') = sum_G' psi(G') c(G''-G'), c = analytic dVloc
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, &pw_rho_, q_cart_, 0);
    const ModuleBase::Vector3<double> g1(0.1, 0.0, 0.0), g2(0.0, 0.1, 0.0);
    const std::vector<std::complex<double>> d0 = data_.get_dpsi(0, 0, 0);
    const std::vector<std::complex<double>> d1 = data_.get_dpsi(0, 0, 1);
    ASSERT_EQ(d0.size(), static_cast<size_t>(kq.get_npwk()));
    ASSERT_EQ(d1.size(), static_cast<size_t>(kq.get_npwk()));
    for (int igl = 0; igl < kq.get_npwk(); ++igl)
    {
        const ModuleBase::Vector3<double> gpp = kq.get_gcar(igl);
        const std::complex<double> e0 = dfpt_test::analytic_dvloc(ucell_, tau_, 0, gpp + q_cart_);
        const std::complex<double> e1 =
            0.7 * dfpt_test::analytic_dvloc(ucell_, tau_, 0, gpp - g1 + q_cart_) +
            std::complex<double>(0.3, 0.2) * dfpt_test::analytic_dvloc(ucell_, tau_, 0, gpp - g2 + q_cart_);
        EXPECT_NEAR(d0[igl].real(), e0.real(), 1.0e-8);
        EXPECT_NEAR(d0[igl].imag(), e0.imag(), 1.0e-8);
        EXPECT_NEAR(d1[igl].real(), e1.real(), 1.0e-8);
        EXPECT_NEAR(d1[igl].imag(), e1.imag(), 1.0e-8);
    }
}

TEST_F(DFPTPertSerialTest, BuildVkbL0MatchesIndependentSimpson)
{
    dfpt_test::make_nc_atom(ucell_);
    const int npwk = pw_wfc_.npwk[0];
    std::vector<ModuleBase::Vector3<double>> gk(npwk);
    for (int ig = 0; ig < npwk; ++ig)
    {
        gk[ig] = pw_wfc_.getgpluskcar(0, ig);
    }
    std::vector<std::vector<std::complex<double>>> vkb;
    pert_.build_vkb(0, 0, gk, vkb);
    ASSERT_EQ(vkb.size(), 4u); // l=0 gives one row, l=1 gives three rows

    const pseudo &p = ucell_.atoms[0].ncpp;
    const double dx = p.rab[0];
    const double pref = ModuleBase::FOUR_PI / std::sqrt(ucell_.omega);
    auto simpson = [&](const std::function<double(int)> &f, int n)
    {
        double s = f(0) + f(n - 1);
        for (int i = 1; i < n - 1; ++i)
        {
            s += f(i) * ((i % 2 == 1) ? 4.0 : 2.0);
        }
        return s * dx / 3.0;
    };

    for (int ig = 0; ig < npwk; ++ig)
    {
        const double g = std::sqrt(gk[ig] * gk[ig]) * ucell_.tpiba; // bohr^-1
        // independent j0 and Simpson transform (no ModuleBase Sphbes/Integral)
        auto f0 = [&](int i)
        {
            const double gr = g * p.r[i];
            const double j0 = (gr < 1.0e-12) ? 1.0 : std::sin(gr) / gr;
            return p.betar(0, i) * j0 * p.r[i];
        };
        const double vq = pref * simpson(f0, p.msh);
        const double arg = -ModuleBase::TWO_PI * (gk[ig] * tau_);
        const std::complex<double> expect =
            0.5 * std::sqrt(1.0 / ModuleBase::PI) * vq * std::complex<double>(std::cos(arg), std::sin(arg));
        EXPECT_NEAR(vkb[0][ig].real(), expect.real(), 1.0e-9 * std::max(1.0, std::abs(expect)));
        EXPECT_NEAR(vkb[0][ig].imag(), expect.imag(), 1.0e-9 * std::max(1.0, std::abs(expect)));
    }

    // a tau shift changes every projector by the pure phase e^{-i 2pi gk.dtau}
    const ModuleBase::Vector3<double> dtau(0.07, -0.11, 0.05);
    ucell_.atoms[0].tau[0] = tau_ + dtau;
    std::vector<std::vector<std::complex<double>>> vkb2;
    pert_.build_vkb(0, 0, gk, vkb2);
    ucell_.atoms[0].tau[0] = tau_;
    for (int mu = 0; mu < 4; ++mu)
    {
        for (int ig = 0; ig < npwk; ++ig)
        {
            if (std::abs(vkb[mu][ig]) < 1.0e-14)
            {
                continue;
            }
            const double arg = -ModuleBase::TWO_PI * (gk[ig] * dtau);
            const std::complex<double> expect(std::cos(arg), std::sin(arg));
            const std::complex<double> ratio = vkb2[mu][ig] / vkb[mu][ig];
            EXPECT_NEAR(ratio.real(), expect.real(), 1.0e-9);
            EXPECT_NEAR(ratio.imag(), expect.imag(), 1.0e-9);
        }
    }
}

namespace
{

// Reconfigure the fixture atom with a third l = 2 projector so the
// d-channel (m = -2..2) rows of the projector tables are live: the
// shared fixture pseudo stops at l = 1, which is exactly why the
// m = +2 gradient factor bug had no test coverage.
void make_nc_atom_with_d_channel(UnitCell &ucell)
{
    dfpt_test::make_nc_atom(ucell);
    Atom &at = ucell.atoms[0];
    pseudo &p = at.ncpp;
    p.nbeta = 3;
    p.lll = {0, 1, 2};
    p.nh = 9;
    p.betar.create(3, 121);
    for (int i = 0; i < 121; ++i)
    {
        const double r = p.r[i];
        p.betar(0, i) = std::exp(-std::pow(r - 1.0, 2) / (2.0 * 0.3 * 0.3));
        p.betar(1, i) = std::exp(-std::pow(r - 1.2, 2) / (2.0 * 0.35 * 0.35));
        p.betar(2, i) = std::exp(-std::pow(r - 0.9, 2) / (2.0 * 0.4 * 0.4));
    }
    p.dion.create(3, 3);
    p.dion(0, 0) = 0.8;
    p.dion(0, 1) = 0.15;
    p.dion(1, 0) = -0.25;
    p.dion(1, 1) = 1.1;
    p.dion(2, 2) = 0.6;
}

} // namespace

TEST_F(DFPTPertSerialTest, BuildVkbDkMatchesCentralDifference)
{
    // QE-style value/gradient consistency lock: the analytic k-derivative
    // of every projector row (angular gradient table + radial difference +
    // structure phase) must equal the central difference of build_vkb
    // itself on a shifted gk list. A hardcoded gradient entry that
    // disagrees with its own value table (the m = +2 factor bug this test
    // regression-locks) fails here at the full amplitude of the channel.
    make_nc_atom_with_d_channel(ucell_);
    const int npwk = pw_wfc_.npwk[0];
    std::vector<ModuleBase::Vector3<double>> gk(npwk);
    for (int ig = 0; ig < npwk; ++ig)
    {
        gk[ig] = pw_wfc_.getgpluskcar(0, ig);
    }
    std::vector<std::vector<std::complex<double>>> vkb;
    pert_.build_vkb(0, 0, gk, vkb);
    ASSERT_EQ(vkb.size(), 9u); // l=0: 1, l=1: 3, l=2: 5 rows

    const double delta = 1.0e-5;
    for (int dir = 0; dir < 3; ++dir)
    {
        std::vector<ModuleBase::Vector3<double>> gk_plus(npwk);
        std::vector<ModuleBase::Vector3<double>> gk_minus(npwk);
        for (int ig = 0; ig < npwk; ++ig)
        {
            gk_plus[ig] = gk[ig];
            gk_minus[ig] = gk[ig];
            gk_plus[ig][dir] += delta;
            gk_minus[ig][dir] -= delta;
        }
        std::vector<std::vector<std::complex<double>>> vkb_plus;
        std::vector<std::vector<std::complex<double>>> vkb_minus;
        pert_.build_vkb(0, 0, gk_plus, vkb_plus);
        pert_.build_vkb(0, 0, gk_minus, vkb_minus);
        std::vector<std::vector<std::complex<double>>> dvkb;
        pert_.build_vkb_dk(0, 0, dir, gk, vkb, dvkb);
        ASSERT_EQ(dvkb.size(), 9u);
        for (int mu = 0; mu < 9; ++mu)
        {
            for (int ig = 0; ig < npwk; ++ig)
            {
                const std::complex<double> fd = (vkb_plus[mu][ig] - vkb_minus[mu][ig]) / (2.0 * delta);
                const double scale = std::max(1.0, std::abs(fd));
                EXPECT_NEAR(dvkb[mu][ig].real(), fd.real(), 1.0e-6 * scale);
                EXPECT_NEAR(dvkb[mu][ig].imag(), fd.imag(), 1.0e-6 * scale);
            }
        }
    }
}

TEST_F(DFPTPertSerialTest, BuildVkbL2MatchesIndependentHarmonics)
{
    // Independent value lock for the d block: the m = +2 row against the
    // textbook real harmonic Y_{2,2} and an independent spherical Bessel
    // j2, mirroring the l = 0 Simpson lock for the s block.
    make_nc_atom_with_d_channel(ucell_);
    const int npwk = pw_wfc_.npwk[0];
    std::vector<ModuleBase::Vector3<double>> gk(npwk);
    for (int ig = 0; ig < npwk; ++ig)
    {
        gk[ig] = pw_wfc_.getgpluskcar(0, ig);
    }
    std::vector<std::vector<std::complex<double>>> vkb;
    pert_.build_vkb(0, 0, gk, vkb);
    ASSERT_EQ(vkb.size(), 9u);

    const pseudo &p = ucell_.atoms[0].ncpp;
    const double dx = p.rab[0];
    const double pref = ModuleBase::FOUR_PI / std::sqrt(ucell_.omega);
    auto simpson = [&](const std::function<double(int)> &f, int n)
    {
        double s = f(0) + f(n - 1);
        for (int i = 1; i < n - 1; ++i)
        {
            s += f(i) * ((i % 2 == 1) ? 4.0 : 2.0);
        }
        return s * dx / 3.0;
    };
    auto j2 = [](double x)
    {
        if (x < 1.0e-12)
        {
            return 0.0;
        }
        return (3.0 / (x * x * x) - 1.0 / x) * std::sin(x) - 3.0 * std::cos(x) / (x * x);
    };

    for (int ig = 0; ig < npwk; ++ig)
    {
        const double gnorm = std::sqrt(gk[ig] * gk[ig]) * ucell_.tpiba;
        auto f2 = [&](int i) { return p.betar(2, i) * j2(gnorm * p.r[i]) * p.r[i]; };
        const double vq = pref * simpson(f2, p.msh);
        const double gmag = std::sqrt(gk[ig] * gk[ig]);
        if (gmag < 1.0e-10)
        {
            continue; // degenerate direction guard, exercised elsewhere
        }
        const double nx = gk[ig].x / gmag;
        const double ny = gk[ig].y / gmag;
        // the l = 2 block occupies mu = 4..8; the m-walk m=3 -> mr=+2 lands
        // on mu = 7, the (x^2 - y^2) channel of the real harmonic
        const double y22 = 0.25 * std::sqrt(15.0 / ModuleBase::PI) * (nx * nx - ny * ny);
        const double arg = -ModuleBase::TWO_PI * (gk[ig] * tau_);
        const std::complex<double> phase(std::cos(arg), std::sin(arg));
        const std::complex<double> expect = -1.0 * y22 * vq * phase; // (-i)^2 = -1
        EXPECT_NEAR(vkb[7][ig].real(), expect.real(), 1.0e-9 * std::max(1.0, std::abs(expect)));
        EXPECT_NEAR(vkb[7][ig].imag(), expect.imag(), 1.0e-9 * std::max(1.0, std::abs(expect)));
    }
}

TEST_F(DFPTPertSerialTest, DVnlDtauMatchesOperatorFiniteDifference)
{
    dfpt_test::make_nc_atom(ucell_);
    const int npwk = pw_wfc_.npwk[0];
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, &pw_rho_, q_cart_, 0);
    const int npwkq = kq.get_npwk();
    std::vector<ModuleBase::Vector3<double>> gk_in(npwk), gk_out(npwkq);
    for (int ig = 0; ig < npwk; ++ig)
    {
        gk_in[ig] = pw_wfc_.getgpluskcar(0, ig);
    }
    for (int igl = 0; igl < npwkq; ++igl)
    {
        gk_out[igl] = kq.get_gpluskq(igl);
    }

    // deterministic pseudo-random wavefunctions, normalized per band
    psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
    unsigned seed = 20260814u;
    auto rnd = [&]()
    {
        seed = seed * 1664525u + 1013904223u;
        return ((seed >> 8) & 0xffffff) / 16777216.0 * 2.0 - 1.0;
    };
    for (int b = 0; b < 2; ++b)
    {
        double nrm = 0.0;
        for (int ig = 0; ig < npwk; ++ig)
        {
            psi(0, b, ig) = std::complex<double>(rnd(), rnd());
            nrm += std::norm(psi(0, b, ig));
        }
        nrm = std::sqrt(nrm);
        for (int ig = 0; ig < npwk; ++ig)
        {
            psi(0, b, ig) /= nrm;
        }
    }

    const int dir = 1; // y displacement
    std::vector<std::vector<std::complex<double>>> dv_psi;
    pert_.dVnl_dtau(0, dir, q_cart_, psi, 0, dv_psi);
    ASSERT_EQ(dv_psi.size(), 2u);
    ASSERT_EQ(dv_psi[0].size(), static_cast<size_t>(npwkq));

    // reference: finite difference of the displaced separable operator
    // <k+q+G''| Vnl(tau) |psi_k> = sum_mu vkb_out(mu,G'') (D becp)_mu
    const pseudo &p = ucell_.atoms[0].ncpp;
    const int nh = p.nh;
    std::vector<int> row_ib, row_m;
    for (int ib = 0; ib < p.nbeta; ++ib)
    {
        for (int m = 0; m < 2 * p.lll[ib] + 1; ++m)
        {
            row_ib.push_back(ib);
            row_m.push_back(m);
        }
    }
    ASSERT_EQ(static_cast<int>(row_ib.size()), nh);

    const double eps = 1.0e-5; // lat0 units
    ModuleBase::Vector3<double> d(0.0, 0.0, 0.0);
    d[dir] = 1.0;
    std::vector<std::vector<std::complex<double>>> vkb_in_p, vkb_in_m, vkb_out_p, vkb_out_m;
    ucell_.atoms[0].tau[0] = tau_ + eps * d;
    pert_.build_vkb(0, 0, gk_in, vkb_in_p);
    pert_.build_vkb(0, 0, gk_out, vkb_out_p);
    ucell_.atoms[0].tau[0] = tau_ - eps * d;
    pert_.build_vkb(0, 0, gk_in, vkb_in_m);
    pert_.build_vkb(0, 0, gk_out, vkb_out_m);
    ucell_.atoms[0].tau[0] = tau_;

    for (int b = 0; b < 2; ++b)
    {
        std::vector<std::complex<double>> fd(npwkq, std::complex<double>(0.0, 0.0));
        for (int side = 0; side < 2; ++side)
        {
            const std::vector<std::vector<std::complex<double>>> &vin = side == 0 ? vkb_in_p : vkb_in_m;
            const std::vector<std::vector<std::complex<double>>> &vout = side == 0 ? vkb_out_p : vkb_out_m;
            std::vector<std::complex<double>> becp(nh, std::complex<double>(0.0, 0.0));
            std::vector<std::complex<double>> dc(nh, std::complex<double>(0.0, 0.0));
            for (int nu = 0; nu < nh; ++nu)
            {
                for (int ig = 0; ig < npwk; ++ig)
                {
                    becp[nu] += std::conj(vin[nu][ig]) * psi(0, b, ig);
                }
            }
            for (int mu = 0; mu < nh; ++mu)
            {
                for (int nu = 0; nu < nh; ++nu)
                {
                    if (row_m[mu] == row_m[nu])
                    {
                        dc[mu] += p.dion(row_ib[mu], row_ib[nu]) * becp[nu];
                    }
                }
            }
            const double s = (side == 0) ? 1.0 : -1.0;
            for (int igl = 0; igl < npwkq; ++igl)
            {
                std::complex<double> m(0.0, 0.0);
                for (int mu = 0; mu < nh; ++mu)
                {
                    m += vout[mu][igl] * dc[mu];
                }
                fd[igl] += s * m;
            }
        }
        for (int igl = 0; igl < npwkq; ++igl)
        {
            fd[igl] /= (2.0 * eps * dfpt_test::lat0);
            EXPECT_NEAR(dv_psi[b][igl].real(), fd[igl].real(), 1.0e-7);
            EXPECT_NEAR(dv_psi[b][igl].imag(), fd[igl].imag(), 1.0e-7);
        }
    }
}

TEST_F(DFPTPertSerialTest, NonlocalPathRejectsUltrasoft)
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    dfpt_test::make_nc_atom(ucell_);
    ucell_.atoms[0].ncpp.tvanp = true;
    const int npwk = pw_wfc_.npwk[0];
    psi::Psi<std::complex<double>> psi(1, 1, npwk, npwk, true);
    psi.zero_out();
    std::vector<std::vector<std::complex<double>>> dv;
    EXPECT_EXIT(pert_.dVnl_dtau(0, 0, q_cart_, psi, 0, dv), ::testing::ExitedWithCode(1), ".*");
}

TEST_F(DFPTPertSerialTest, BuildDvWithInactiveDftuIsPurePW)
{
    // a wired but unusable provider (occupation matrices not initialized)
    // must not change the assembled first-order potential (U0 reservation);
    // DFPT_PW::init additionally rejects a wired provider outright
    ModuleDFPT::DFPT_PW_Data data_plain;
    data_plain.init(&qlist_, 1, 2, pw_wfc_.npwk_max, pw_rho_.nrxx, 1, 1, nullptr);
    pert_.build_dv(0, 0, 1, data_plain);

    Plus_U_Base dftu;
    ModuleDFPT::DFPT_PW_Data data_u;
    data_u.init(&qlist_, 1, 2, pw_wfc_.npwk_max, pw_rho_.nrxx, 1, 1, &dftu);
    EXPECT_TRUE(data_u.with_u());
    EXPECT_FALSE(data_u.u_active());
    pert_.build_dv(0, 0, 1, data_u);

    const std::vector<std::complex<double>> dv_a = data_plain.get_dv_rc(0, 0);
    const std::vector<std::complex<double>> dv_b = data_u.get_dv_rc(0, 0);
    ASSERT_EQ(dv_a.size(), static_cast<size_t>(pw_rho_.nrxx));
    ASSERT_EQ(dv_a.size(), dv_b.size());
    for (size_t i = 0; i < dv_a.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(dv_a[i].real(), dv_b[i].real());
        EXPECT_DOUBLE_EQ(dv_a[i].imag(), dv_b[i].imag());
    }
}
