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

#define private public
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/pseudo.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_cell/magnetism.h"
#include "source_pw/module_pwdft/stru_fac.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#undef private

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_psi/psi.h"

// test-support ctor/dtor stubs (see test/dfpt_pw_run_test.cpp); the DFPT
// serial test only needs default-constructible cell/sf objects.
pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}
Atom::Atom()
{
}
Atom::~Atom()
{
}
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
Structure_Factor::Structure_Factor()
{
}
Structure_Factor::~Structure_Factor()
{
}

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
 *   - build_efield: -E.r ramp layout and its closed-form discrete FT.
 *   - build_vkb: l=0 projector against an independent hand-rolled Simpson
 *     transform; pure phase response under a tau shift.
 *   - dVnl_dtau: the two-term NC identity is validated against a finite
 *     difference of the displaced separable operator; USPP is rejected.
 *   - build_dv under with_u()/u_active()==false (pure-PW DFT+U safety).
 */

class DFPTPertSerialTest : public testing::Test
{
  protected:
    const double lat0_ = 1.8897261254578281;
    const double ecutwfc_ = 2.5; // Ry
    // rho cutoff inflated to 9x ecutwfc so every Delta = G''-G' of the
    // convolution lies inside the rho ball and nothing aliases
    const double rho_mult_ = 9.0;

    ModuleBase::Matrix3 latvec_;
    UnitCell ucell_;
    ModulePW::PW_Basis pw_rho_;
    ModulePW::PW_Basis_K pw_wfc_;
    Structure_Factor sf_;
    ModuleDFPT::DFPT_Pert pert_;
    ModuleCell::QList qlist_;
    ModuleDFPT::DFPT_PW_Data data_;

    // q is generic; k = -q so k+q = 0: the k+q ball then stays inside the
    // ground-state G list (single-k limitation documented in DFPT_KQ_Basis)
    const ModuleBase::Vector3<double> q_d_{0.13, 0.0, 0.07};
    const ModuleBase::Vector3<double> k_d_{-0.13, 0.0, -0.07};
    ModuleBase::Vector3<double> q_cart_;
    const ModuleBase::Vector3<double> tau_{1.1, 2.3, 0.7}; // lat0 units

    void SetUp() override
    {
        latvec_ = ModuleBase::Matrix3(10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0);
        ucell_.ntype = 1;
        ucell_.nat = 1;
        ucell_.atoms = new Atom[1];
        ucell_.atoms[0].na = 1;
        ucell_.atoms[0].tau.resize(1);
        ucell_.atoms[0].tau[0] = tau_;
        ucell_.latvec = latvec_;
        ucell_.GT = latvec_.Inverse();
        ucell_.G = ucell_.GT.Transpose();
        ucell_.lat0 = lat0_;
        ucell_.tpiba = ModuleBase::TWO_PI / lat0_;
        ucell_.tpiba2 = ucell_.tpiba * ucell_.tpiba;
        ucell_.omega = 1000.0 * lat0_ * lat0_ * lat0_;
        MakeCoulombAtom();

        // shared-grid basis setup, mirroring setup_pwrho / setup_pwwfc
        pw_rho_.initgrids(lat0_, latvec_, rho_mult_ * ecutwfc_);
        pw_rho_.initparameters(false, rho_mult_ * ecutwfc_);
        pw_rho_.fft_bundle.initfftmode(0);
        pw_rho_.setuptransform();
        pw_rho_.collect_local_pw();

        const ModuleBase::Vector3<double> klist[1] = {k_d_};
        pw_wfc_.initgrids(lat0_, latvec_, pw_rho_.nx, pw_rho_.ny, pw_rho_.nz);
        pw_wfc_.initparameters(false, ecutwfc_, 1, klist);
        pw_wfc_.fft_bundle.initfftmode(0);
        pw_wfc_.setuptransform();
        pw_wfc_.collect_local_pw();

        qlist_.nkstot = 1;
        qlist_.kvec_d.push_back(q_d_);
        q_cart_ = q_d_ * ucell_.G;

        data_.init(&qlist_, 1, 2, pw_wfc_.npwk_max, pw_rho_.nrxx, 1, 1, nullptr);
        pert_.init(ucell_, &pw_rho_, &pw_wfc_, sf_);
    }

    void TearDown() override
    {
        delete[] ucell_.atoms;
        ucell_.atoms = nullptr;
    }

    void MakeCoulombAtom()
    {
        Atom& at = ucell_.atoms[0];
        at.label = "C";
        at.coulomb_potential = true;
        at.ncpp.zv = 4.0;
        at.ncpp.tvanp = false;
        at.ncpp.has_so = false;
        at.ncpp.nbeta = 0;
        at.ncpp.nh = 0;
        at.ncpp.msh = 0;
        at.ncpp.kkbeta = 0;
    }

    void MakeNCAtom()
    {
        Atom& at = ucell_.atoms[0];
        at.label = "Si";
        at.coulomb_potential = false;
        pseudo& p = at.ncpp;
        p.zv = 4.0;
        p.tvanp = false;
        p.has_so = false;
        p.nbeta = 2;
        p.lll = {0, 1};
        p.nh = 4;
        p.msh = 121;
        p.kkbeta = 121;
        p.r.resize(121);
        p.rab.resize(121);
        p.vloc_at.assign(121, 0.0);
        const double dx = 0.025;
        for (int i = 0; i < 121; ++i)
        {
            p.r[i] = i * dx;
            p.rab[i] = dx;
        }
        p.betar.create(2, 121);
        for (int i = 0; i < 121; ++i)
        {
            const double r = p.r[i];
            p.betar(0, i) = std::exp(-std::pow(r - 1.0, 2) / (2.0 * 0.3 * 0.3));
            p.betar(1, i) = std::exp(-std::pow(r - 1.2, 2) / (2.0 * 0.35 * 0.35));
        }
        p.dion.create(2, 2);
        p.dion(0, 0) = 0.8;
        p.dion(0, 1) = 0.15;
        p.dion(1, 0) = -0.25;
        p.dion(1, 1) = 1.1;
    }

    // key of an integer FFT triple (gcar * a is integral on the cubic cell)
    long long FKey(int ix, int iy, int iz) const
    {
        return (static_cast<long long>(ix + 64) * 128 + (iy + 64)) * 128 + (iz + 64);
    }
    long long GKey(const ModuleBase::Vector3<double>& g) const
    {
        const double a = 10.0;
        return FKey(static_cast<int>(std::llround(g.x * a)),
                    static_cast<int>(std::llround(g.y * a)),
                    static_cast<int>(std::llround(g.z * a)));
    }

    // analytic Coulomb local potential (Ry) at |g|^2 in bohr^-2, mirroring
    // vl_pw.cpp::vloc_coulomb independently of DFPT_Pert::vloc_at_g
    double VlocCoulomb(double g2_bohr) const
    {
        return -ucell_.atoms[0].ncpp.zv * ModuleBase::e2 * ModuleBase::FOUR_PI / ucell_.omega / g2_bohr;
    }

    // analytic dVloc/dtau_alpha coefficient at displacement vector w (1/lat0)
    std::complex<double> AnalyticDVloc(int dir, const ModuleBase::Vector3<double>& w) const
    {
        const double w2 = w * w;
        if (w2 < 1.0e-12)
        {
            return std::complex<double>(0.0, 0.0);
        }
        const double arg = ModuleBase::TWO_PI * (w * tau_);
        return std::complex<double>(0.0, 1.0) * (ucell_.tpiba * w[dir]) * VlocCoulomb(w2 * ucell_.tpiba2)
               * std::complex<double>(std::cos(arg), std::sin(arg));
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
                // finite difference of vloc(|Delta+q|) e^{i 2pi (Delta+q).tau}
                // per bohr of displacement
                const double ap = ModuleBase::TWO_PI * (w * (tau_ + eps * d));
                const double am = ModuleBase::TWO_PI * (w * (tau_ - eps * d));
                const std::complex<double> fd = VlocCoulomb((w * w) * ucell_.tpiba2)
                                                * (std::polar(1.0, ap) - std::polar(1.0, am))
                                                / (2.0 * eps * lat0_);
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
        const std::complex<double> expect = AnalyticDVloc(0, pw_rho_.gcar[ig] + q_cart_);
        EXPECT_NEAR(dv_recip[ig].real(), expect.real(), 1.0e-9);
        EXPECT_NEAR(dv_recip[ig].imag(), expect.imag(), 1.0e-9);
    }

    // wavefunctions: band 0 = single plane wave at k, band 1 = two components
    const int npwk = pw_wfc_.npwk[0];
    const ModuleBase::Vector3<double> kc = pw_wfc_.kvec_c[0];
    int ig_zero = -1, ig_x = -1, ig_y = -1;
    for (int ig = 0; ig < npwk; ++ig)
    {
        const long long key = GKey(pw_wfc_.getgpluskcar(0, ig) - kc);
        if (key == FKey(0, 0, 0))
        {
            ig_zero = ig;
        }
        else if (key == FKey(1, 0, 0))
        {
            ig_x = ig;
        }
        else if (key == FKey(0, 1, 0))
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
    kq.init(&pw_wfc_, q_cart_, 0);
    const ModuleBase::Vector3<double> g1(0.1, 0.0, 0.0), g2(0.0, 0.1, 0.0);
    const std::vector<std::complex<double>> d0 = data_.get_dpsi(0, 0, 0);
    const std::vector<std::complex<double>> d1 = data_.get_dpsi(0, 0, 1);
    ASSERT_EQ(d0.size(), static_cast<size_t>(kq.get_npwk()));
    ASSERT_EQ(d1.size(), static_cast<size_t>(kq.get_npwk()));
    for (int igl = 0; igl < kq.get_npwk(); ++igl)
    {
        const ModuleBase::Vector3<double> gpp = kq.get_gcar(igl);
        const std::complex<double> e0 = AnalyticDVloc(0, gpp + q_cart_);
        const std::complex<double> e1 = 0.7 * AnalyticDVloc(0, gpp - g1 + q_cart_)
                                      + std::complex<double>(0.3, 0.2) * AnalyticDVloc(0, gpp - g2 + q_cart_);
        EXPECT_NEAR(d0[igl].real(), e0.real(), 1.0e-8);
        EXPECT_NEAR(d0[igl].imag(), e0.imag(), 1.0e-8);
        EXPECT_NEAR(d1[igl].real(), e1.real(), 1.0e-8);
        EXPECT_NEAR(d1[igl].imag(), e1.imag(), 1.0e-8);
    }
}

TEST_F(DFPTPertSerialTest, BuildEfieldRampMatchesClosedForm)
{
    const double E0 = 0.5;
    const double Lx = 10.0 * lat0_;
    pert_.build_efield(ModuleBase::Vector3<double>(E0, 0.0, 0.0), data_);
    const std::vector<std::complex<double>> dv = data_.get_dv_rc(0, 0);
    ASSERT_EQ(dv.size(), static_cast<size_t>(pw_rho_.nrxx));
    const int nx = pw_rho_.nx;

    // real-space values form a pure x ramp -E0 (ix/nx) Lx
    for (int ir = 0; ir < pw_rho_.nrxx; ++ir)
    {
        const double j = std::llround(-dv[ir].real() / (E0 * Lx) * nx);
        EXPECT_GE(j, 0.0);
        EXPECT_LE(j, nx - 1.0);
        EXPECT_NEAR(dv[ir].real(), -E0 * (j / nx) * Lx, 1.0e-10);
        EXPECT_EQ(dv[ir].imag(), 0.0);
    }

    // closed-form discrete FT of the sawtooth discriminates the (ix,iy,iz)
    // layout: c(k) = E0 Lx / (nx (1 - w)), w = exp(-2 pi i k / nx), k != 0
    std::vector<std::complex<double>> recip(pw_rho_.npw);
    pw_rho_.real2recip(dv.data(), recip.data());
    int ig_p100 = -1, ig_m100 = -1, ig_p110 = -1, ig_p010 = -1;
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        const long long key = GKey(pw_rho_.gcar[ig]);
        if (key == FKey(1, 0, 0))
        {
            ig_p100 = ig;
        }
        else if (key == FKey(-1, 0, 0))
        {
            ig_m100 = ig;
        }
        else if (key == FKey(1, 1, 0))
        {
            ig_p110 = ig;
        }
        else if (key == FKey(0, 1, 0))
        {
            ig_p010 = ig;
        }
    }
    ASSERT_GE(ig_p100, 0);
    ASSERT_GE(ig_m100, 0);
    ASSERT_GE(ig_p110, 0);
    ASSERT_GE(ig_p010, 0);
    const std::complex<double> w1 = std::polar(1.0, -ModuleBase::TWO_PI / nx);
    const std::complex<double> nx_c(nx);
    const std::complex<double> c1 = E0 * Lx / (nx_c * (1.0 - w1));
    const std::complex<double> cm1 = E0 * Lx / (nx_c * (1.0 - std::conj(w1)));
    EXPECT_NEAR(recip[ig_p100].real(), c1.real(), 1.0e-8);
    EXPECT_NEAR(recip[ig_p100].imag(), c1.imag(), 1.0e-8);
    EXPECT_NEAR(recip[ig_m100].real(), cm1.real(), 1.0e-8);
    EXPECT_NEAR(recip[ig_m100].imag(), cm1.imag(), 1.0e-8);
    EXPECT_NEAR(std::abs(recip[ig_p110]), 0.0, 1.0e-10);
    EXPECT_NEAR(std::abs(recip[ig_p010]), 0.0, 1.0e-10);
}

TEST_F(DFPTPertSerialTest, BuildVkbL0MatchesIndependentSimpson)
{
    MakeNCAtom();
    const int npwk = pw_wfc_.npwk[0];
    std::vector<ModuleBase::Vector3<double>> gk(npwk);
    for (int ig = 0; ig < npwk; ++ig)
    {
        gk[ig] = pw_wfc_.getgpluskcar(0, ig);
    }
    std::vector<std::vector<std::complex<double>>> vkb;
    pert_.build_vkb(0, 0, gk, vkb);
    ASSERT_EQ(vkb.size(), 4u); // l=0 gives one row, l=1 gives three rows

    const pseudo& p = ucell_.atoms[0].ncpp;
    const double dx = p.rab[0];
    const double pref = ModuleBase::FOUR_PI / std::sqrt(ucell_.omega);
    auto simpson = [&](const std::function<double(int)>& f, int n)
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
        const double arg = ModuleBase::TWO_PI * (gk[ig] * tau_);
        const std::complex<double> expect = 0.5 * std::sqrt(1.0 / ModuleBase::PI) * vq
                                            * std::complex<double>(std::cos(arg), std::sin(arg));
        EXPECT_NEAR(vkb[0][ig].real(), expect.real(), 1.0e-9 * std::max(1.0, std::abs(expect)));
        EXPECT_NEAR(vkb[0][ig].imag(), expect.imag(), 1.0e-9 * std::max(1.0, std::abs(expect)));
    }

    // a tau shift changes every projector by the pure phase e^{i 2pi gk.dtau}
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
            const double arg = ModuleBase::TWO_PI * (gk[ig] * dtau);
            const std::complex<double> expect(std::cos(arg), std::sin(arg));
            const std::complex<double> ratio = vkb2[mu][ig] / vkb[mu][ig];
            EXPECT_NEAR(ratio.real(), expect.real(), 1.0e-9);
            EXPECT_NEAR(ratio.imag(), expect.imag(), 1.0e-9);
        }
    }
}

TEST_F(DFPTPertSerialTest, DVnlDtauMatchesOperatorFiniteDifference)
{
    MakeNCAtom();
    const int npwk = pw_wfc_.npwk[0];
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, q_cart_, 0);
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
    const pseudo& p = ucell_.atoms[0].ncpp;
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
            const auto& vin = side == 0 ? vkb_in_p : vkb_in_m;
            const auto& vout = side == 0 ? vkb_out_p : vkb_out_m;
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
            fd[igl] /= (2.0 * eps * lat0_);
            EXPECT_NEAR(dv_psi[b][igl].real(), fd[igl].real(), 1.0e-7);
            EXPECT_NEAR(dv_psi[b][igl].imag(), fd[igl].imag(), 1.0e-7);
        }
    }
}

TEST_F(DFPTPertSerialTest, NonlocalPathRejectsUltrasoft)
{
    MakeNCAtom();
    ucell_.atoms[0].ncpp.tvanp = true;
    const int npwk = pw_wfc_.npwk[0];
    psi::Psi<std::complex<double>> psi(1, 1, npwk, npwk, true);
    psi.zero_out();
    std::vector<std::vector<std::complex<double>>> dv;
    EXPECT_EXIT(pert_.dVnl_dtau(0, 0, q_cart_, psi, 0, dv), ::testing::ExitedWithCode(1), ".*");
}

TEST_F(DFPTPertSerialTest, BuildDvWithInactiveDftuIsPurePW)
{
    // a wired but unusable Plus_U (locale uninitialized) must not change the
    // assembled first-order potential (U0 reservation, pure-PW degradation)
    ModuleDFPT::DFPT_PW_Data data_plain;
    data_plain.init(&qlist_, 1, 2, pw_wfc_.npwk_max, pw_rho_.nrxx, 1, 1, nullptr);
    pert_.build_dv(0, 0, 1, data_plain);

    Plus_U dftu;
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
