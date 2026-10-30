#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <functional>
#include <vector>

// serial unit test of the first-order density response (C3).
// Everything here runs without __MPI: the plane-wave bases are built through
// the real serial initgrids/initparameters/setuptransform path on a shared
// FFT grid, exactly like the production setup_pwrho/setup_pwwfc sequence.

#include "dfpt_serial_fixture.h"
#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_psi/psi.h"
#include "source_pw/module_dfpt/dfpt_kq_basis.h"
#include "source_pw/module_dfpt/dfpt_pw_data.h"
#include "source_pw/module_dfpt/dfpt_rho.h"

/************************************************
 *  serial unit test of DFPT_Rho (C3)
 ***********************************************/

/**
 * - Tested Functions:
 *   - DFPT_Rho::compute_drho - the q-shifted response density coefficients
 *     A_Delta = sum_{kn occ} wg sum_G c*_G d_{G+Delta} against a brute-force
 *     G-space double sum built from independent G-keyed coefficient maps,
 *     and the real-space manifest 2 Re[e^{iqr} u* du] against direct sums.
 *   - charge conservation: the Delta = -q harmonic is dropped at q = Gamma
 *     (drho_g[ig(|g|=0)] == 0, grid sum of drho_r ~ 0).
 *   - DFPT_Rho::mix_drho - plain mixing drho_in + beta (out - in) with a
 *     zero first input, and the second-step combination; residual formula.
 *   - occupation gate: bands with wg < 1e-8 do not contribute.
 */

namespace
{

unsigned g_seed = 20260815u;
double test_rand()
{
    g_seed = g_seed * 1664525u + 1013904223u;
    return ((g_seed >> 8) & 0xffffff) / 16777216.0 * 2.0 - 1.0;
}

std::complex<double> crand()
{
    return std::complex<double>(test_rand(), test_rand());
}

} // namespace

class DFPTRhoSerialTest : public DFPTSerialBase
{
  protected:
    ModuleDFPT::DFPT_Rho rho_;

    static const int nbands_ = 2;

    void SetUp() override
    {
        DFPTSerialBase::SetUp();
        rho_.init({1, pw_rho_.nrxx, &pw_rho_, &pw_wfc_, G_, std::string("rho_plain"), 0.4, 0.0});
    }

    void FillRandomStates(psi::Psi<std::complex<double>>& psi,
                          std::vector<std::vector<std::vector<std::complex<double>>>>& dpsi)
    {
        psi::Psi<std::complex<double>> p(1, nbands_, pw_wfc_.npwk_max, pw_wfc_.npwk[0], true);
        ModuleDFPT::DFPT_KQ_Basis kq;
        kq.init(&pw_wfc_, &pw_rho_, q_cart_, 0);
        dpsi.assign(1, std::vector<std::vector<std::complex<double>>>(nbands_));
        for (int ib = 0; ib < nbands_; ++ib)
        {
            for (int igl = 0; igl < pw_wfc_.npwk[0]; ++igl)
            {
                p(0, ib, igl) = crand();
            }
            dpsi[0][ib].assign(kq.get_npwk(), std::complex<double>(0.0, 0.0));
            for (int igl = 0; igl < kq.get_npwk(); ++igl)
            {
                dpsi[0][ib][igl] = crand();
            }
            data_.set_dpsi(0, 0, ib, dpsi[0][ib]);
        }
        psi = p;
    }
};

TEST_F(DFPTRhoSerialTest, ComputeDrhoMatchesBruteForceGSpace)
{
    psi::Psi<std::complex<double>> psi;
    std::vector<std::vector<std::vector<std::complex<double>>>> dpsi;
    FillRandomStates(psi, dpsi);
    ModuleBase::matrix wg(1, nbands_);
    wg(0, 0) = 1.0;
    wg(0, 1) = 0.0; // unoccupied: must not contribute

    rho_.compute_drho(psi, wg, 0, data_);
    const std::vector<std::complex<double>> drho_g = data_.get_drho_g(0, 0);
    ASSERT_EQ(drho_g.size(), static_cast<size_t>(pw_rho_.npw));

    // enumerate the occupied-band coefficient list G -> c (PW_Basis_K::gcar
    // is per-k: entry igl pairs with psi(0,0,igl)) and the k+q list
    std::vector<ModuleBase::Vector3<double>> glist;
    std::vector<std::complex<double>> clist;
    for (int igl = 0; igl < pw_wfc_.npwk[0]; ++igl)
    {
        glist.push_back(pw_wfc_.gcar[igl]);
        clist.push_back(psi(0, 0, igl));
    }
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, &pw_rho_, q_cart_, 0);
    const std::vector<std::complex<double>> dvec = data_.get_dpsi(0, 0, 0);

    double err2 = 0.0;
    double ref2 = 0.0;
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        // rho-grid ig -> cartesian Delta through the stick/cell position
        const int isz = pw_rho_.ig2isz[ig];
        const int iz = isz % pw_rho_.nz;
        const int is = isz / pw_rho_.nz;
        const int ixy = pw_rho_.is2fftixy[is];
        const int ix = ixy / pw_rho_.fftny;
        const int iy = ixy % pw_rho_.fftny;
        const int mx = (ix <= pw_rho_.nx / 2) ? ix : ix - pw_rho_.nx;
        const int my = (iy <= pw_rho_.ny / 2) ? iy : iy - pw_rho_.ny;
        const int mz = (iz <= pw_rho_.nz / 2) ? iz : iz - pw_rho_.nz;
        const ModuleBase::Vector3<double> delta = ModuleBase::Vector3<double>(mx, my, mz) * G_;
        // A_Delta = (2 w / omega) * sum_G c*_G d_{G+Delta}: the spin factor
        // 2 sits in the band weight w1 = 2 w / omega (the QE incdrhoscf
        // convention, a915352cd), brute-forced over the lists
        std::complex<double> aref(0.0, 0.0);
        for (int jgl = 0; jgl < kq.get_npwk(); ++jgl)
        {
            const ModuleBase::Vector3<double> gq = kq.get_gcar(jgl);
            for (size_t j = 0; j < glist.size(); ++j)
            {
                if (std::abs(glist[j].x - (gq.x - delta.x)) < 1.0e-6 && std::abs(glist[j].y - (gq.y - delta.y)) < 1.0e-6
                    && std::abs(glist[j].z - (gq.z - delta.z)) < 1.0e-6)
                {
                    aref += std::conj(clist[j]) * dvec[jgl];
                    break;
                }
            }
        }
        aref *= 2.0 * wg(0, 0) / pw_rho_.omega;
        err2 += std::norm(drho_g[ig] - aref);
        ref2 += std::norm(aref);
    }
    EXPECT_LT(std::sqrt(err2 / ref2), 1.0e-10);
}

TEST_F(DFPTRhoSerialTest, ComputeDrhoRealSpaceMatchesDirectSum)
{
    psi::Psi<std::complex<double>> psi;
    std::vector<std::vector<std::vector<std::complex<double>>>> dpsi;
    FillRandomStates(psi, dpsi);
    ModuleBase::matrix wg(1, nbands_);
    wg(0, 0) = 1.0;
    wg(0, 1) = 0.0;

    rho_.compute_drho(psi, wg, 0, data_);
    const std::vector<double> drho_r = data_.get_drho_r(0, 0);
    ASSERT_EQ(drho_r.size(), static_cast<size_t>(pw_rho_.nrxx));

    // direct real-space sums u(r) = sum c_G e^{i 2pi G.r} at sample points;
    // PW_Basis_K::gcar is a per-k array indexed by ik*npwk_max + igl, so for
    // k=0 the entry igl pairs directly with the coefficient psi(0,0,igl)
    std::vector<ModuleBase::Vector3<double>> glist;
    std::vector<std::complex<double>> clist;
    for (int igl = 0; igl < pw_wfc_.npwk[0]; ++igl)
    {
        glist.push_back(pw_wfc_.gcar[igl]);
        clist.push_back(psi(0, 0, igl));
    }
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_wfc_, &pw_rho_, q_cart_, 0);
    std::vector<ModuleBase::Vector3<double>> dgl;
    std::vector<std::complex<double>> dcl;
    for (int jgl = 0; jgl < kq.get_npwk(); ++jgl)
    {
        dgl.push_back(kq.get_gcar(jgl));
        dcl.push_back(data_.get_dpsi(0, 0, 0)[jgl]);
    }

    const int samples[5][3] = {{0, 0, 0}, {1, 0, 0}, {0, 2, 0}, {0, 0, 3}, {2, 1, 1}};
    // cartesian grid point: r = (fx,fy,fz) . latvec (row-vector convention)
    for (int s = 0; s < 5; ++s)
    {
        const int ix = samples[s][0] % pw_rho_.nx;
        const int iy = samples[s][1] % pw_rho_.ny;
        const int iz = samples[s][2] % pw_rho_.nz;
        const double fx = static_cast<double>(ix) / pw_rho_.nx;
        const double fy = static_cast<double>(iy) / pw_rho_.ny;
        const double fz = static_cast<double>(iz) / pw_rho_.nz;
        const ModuleBase::Vector3<double> r_cart = ModuleBase::Vector3<double>(fx, fy, fz) * latvec_;
        std::complex<double> u(0.0, 0.0);
        for (size_t j = 0; j < glist.size(); ++j)
        {
            const double ph = ModuleBase::TWO_PI * (glist[j] * r_cart);
            u += clist[j] * std::complex<double>(std::cos(ph), std::sin(ph));
        }
        std::complex<double> du(0.0, 0.0);
        for (size_t j = 0; j < dgl.size(); ++j)
        {
            const double ph = ModuleBase::TWO_PI * (dgl[j] * r_cart);
            du += dcl[j] * std::complex<double>(std::cos(ph), std::sin(ph));
        }
        const double phq = ModuleBase::TWO_PI * (q_d_.x * fx + q_d_.y * fy + q_d_.z * fz);
        const std::complex<double> eq(std::cos(phq), std::sin(phq));
        // manifest density 2 Re[e^{iqr} A(r)] with A carrying the band
        // weight w1 = 2 w / omega (spin factor 2, a915352cd): the outer 2
        // Re and the inner 2 w / omega combine to 4 w / omega
        const double ref = 4.0 * (wg(0, 0) / pw_rho_.omega) * (std::conj(u) * du * eq).real();
        const int ir = (ix * pw_rho_.ny + iy) * pw_rho_.nz + iz;
        EXPECT_NEAR(drho_r[ir], ref, 1.0e-9);
    }
}

TEST_F(DFPTRhoSerialTest, ChargeConservationAtGamma)
{
    // rebuild the fixture bases with q = k = 0
    ModuleBase::Vector3<double> klist0[1] = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};
    ModulePW::PW_Basis_K pw_wfc0;
    pw_wfc0.initgrids(lat0_, latvec_, pw_rho_.nx, pw_rho_.ny, pw_rho_.nz);
    pw_wfc0.initparameters(false, ecutwfc_, 1, klist0);
    pw_wfc0.fft_bundle.initfftmode(0);
    pw_wfc0.setuptransform();
    pw_wfc0.collect_local_pw();

    ModuleCell::QList qlist0;
    qlist0.nkstot = 1;
    qlist0.kvec_d.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));

    ModuleDFPT::DFPT_PW_Data data0;
    data0.init(&qlist0, 1, nbands_, pw_wfc0.npwk_max, pw_rho_.nrxx, 1, 1, nullptr);
    ModuleDFPT::DFPT_Rho rho0;
    rho0.init({1, pw_rho_.nrxx, &pw_rho_, &pw_wfc0, G_, std::string("rho_plain"), 0.4, 0.0});

    psi::Psi<std::complex<double>> psi(1, nbands_, pw_wfc0.npwk_max, pw_wfc0.npwk[0], true);
    ModuleDFPT::DFPT_KQ_Basis kq0;
    kq0.init(&pw_wfc0, &pw_rho_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0);
    for (int ib = 0; ib < nbands_; ++ib)
    {
        for (int igl = 0; igl < pw_wfc0.npwk[0]; ++igl)
        {
            psi(0, ib, igl) = crand();
        }
        std::vector<std::complex<double>> dv(kq0.get_npwk());
        for (int igl = 0; igl < kq0.get_npwk(); ++igl)
        {
            dv[igl] = crand();
        }
        data0.set_dpsi(0, 0, ib, dv);
    }
    ModuleBase::matrix wg(1, nbands_);
    wg(0, 0) = 2.0;
    wg(0, 1) = 0.0;

    rho0.compute_drho(psi, wg, 0, data0);
    const std::vector<std::complex<double>> drho_g = data0.get_drho_g(0, 0);
    ASSERT_EQ(drho_g.size(), static_cast<size_t>(pw_rho_.npw));
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        if (pw_rho_.gcar[ig].norm() < 1.0e-10)
        {
            EXPECT_EQ(drho_g[ig], std::complex<double>(0.0, 0.0));
        }
    }
    // the manifest density integrates (grid-sums) to zero
    const std::vector<double> drho_r = data0.get_drho_r(0, 0);
    double sum = 0.0;
    double absmax = 0.0;
    for (int ir = 0; ir < pw_rho_.nrxx; ++ir)
    {
        sum += drho_r[ir];
        absmax = std::max(absmax, std::abs(drho_r[ir]));
    }
    EXPECT_LT(std::abs(sum) / (absmax * pw_rho_.nrxx), 1.0e-12);
}

TEST_F(DFPTRhoSerialTest, MixDrhoFirstStepIsScaledOutput)
{
    psi::Psi<std::complex<double>> psi;
    std::vector<std::vector<std::vector<std::complex<double>>>> dpsi;
    FillRandomStates(psi, dpsi);
    ModuleBase::matrix wg(1, nbands_);
    wg(0, 0) = 1.0;
    wg(0, 1) = 0.0;

    rho_.compute_drho(psi, wg, 0, data_);
    const std::vector<std::complex<double>> out = data_.get_drho_g(0, 0);
    rho_.mix_drho(0, data_);
    const std::vector<std::complex<double>> mixed = data_.get_drho_g(0, 0);
    ASSERT_EQ(mixed.size(), out.size());
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        EXPECT_NEAR(mixed[ig].real(), 0.4 * out[ig].real(), 1.0e-12);
        EXPECT_NEAR(mixed[ig].imag(), 0.4 * out[ig].imag(), 1.0e-12);
    }
    // zero input: ||out - 0|| / ||out|| == 1 exactly
    EXPECT_NEAR(rho_.get_residual(0, data_), 1.0, 1.0e-12);
}

TEST_F(DFPTRhoSerialTest, MixDrhoSecondStepCombinesCorrectly)
{
    psi::Psi<std::complex<double>> psi;
    std::vector<std::vector<std::vector<std::complex<double>>>> dpsi;
    FillRandomStates(psi, dpsi);
    ModuleBase::matrix wg(1, nbands_);
    wg(0, 0) = 1.0;
    wg(0, 1) = 0.0;

    rho_.compute_drho(psi, wg, 0, data_);
    rho_.mix_drho(0, data_);
    const std::vector<std::complex<double>> in1 = data_.get_drho_g(0, 0);

    // new response from fresh dpsi
    g_seed = 777u;
    std::vector<std::vector<std::vector<std::complex<double>>>> dpsi2;
    FillRandomStates(psi, dpsi2);
    rho_.compute_drho(psi, wg, 0, data_);
    const std::vector<std::complex<double>> out2 = data_.get_drho_g(0, 0);
    rho_.mix_drho(0, data_);
    const std::vector<std::complex<double>> mixed2 = data_.get_drho_g(0, 0);

    double dn2 = 0.0;
    double o2 = 0.0;
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        const std::complex<double> ref = in1[ig] + 0.4 * (out2[ig] - in1[ig]);
        EXPECT_NEAR(mixed2[ig].real(), ref.real(), 1.0e-12);
        EXPECT_NEAR(mixed2[ig].imag(), ref.imag(), 1.0e-12);
        dn2 += std::norm(out2[ig] - in1[ig]);
        o2 += std::norm(out2[ig]);
    }
    EXPECT_NEAR(rho_.get_residual(0, data_), std::sqrt(dn2 / o2), 1.0e-12);
}

TEST_F(DFPTRhoSerialTest, VHartreeQClosedFormAndZeroMode)
{
    // single-G amplitude: dv_ha_g[ig] = e2 4 pi / (tpiba2 |G+q|^2) drho_g[ig]
    const int ig_star = [this]() {
        for (int ig = 0; ig < pw_rho_.npw; ++ig)
        {
            if ((pw_rho_.gcar[ig] + q_cart_) * (pw_rho_.gcar[ig] + q_cart_) > 1.0e-4)
            {
                return ig;
            }
        }
        return -1;
    }();
    ASSERT_GE(ig_star, 0);

    std::vector<std::complex<double>> drho_g(pw_rho_.npw, std::complex<double>(0.0, 0.0));
    drho_g[ig_star] = std::complex<double>(0.3, -0.7);
    std::vector<std::complex<double>> dv;
    rho_.v_hartree_q(q_cart_, drho_g, dv);
    ASSERT_EQ(dv.size(), static_cast<size_t>(pw_rho_.npw));
    const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig_star] + q_cart_;
    const std::complex<double> expect
        = ModuleBase::e2 * ModuleBase::FOUR_PI / (pw_rho_.tpiba2 * (w * w)) * drho_g[ig_star];
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        if (ig == ig_star)
        {
            EXPECT_NEAR(dv[ig].real(), expect.real(), 1.0e-10);
            EXPECT_NEAR(dv[ig].imag(), expect.imag(), 1.0e-10);
        }
        else
        {
            EXPECT_EQ(dv[ig], std::complex<double>(0.0, 0.0));
        }
    }

    // |G+q| = 0 (ig = -q) is skipped like v_hartree skips ig_gge0
    const ModuleBase::Vector3<double> q_minus = pw_rho_.gcar[ig_star] * (-1.0);
    std::vector<std::complex<double>> dv0;
    rho_.v_hartree_q(q_minus, drho_g, dv0);
    EXPECT_EQ(dv0[ig_star], std::complex<double>(0.0, 0.0));

    // a wrong-size drho clears the output instead of aliasing it
    std::vector<std::complex<double>> short_input(3, std::complex<double>(1.0, 1.0));
    rho_.v_hartree_q(q_cart_, short_input, dv0);
    EXPECT_TRUE(dv0.empty());
}

TEST_F(DFPTRhoSerialTest, MixDrhoKerkerFirstStepIsPreconditionedScaledOutput)
{
    // first step from the zero input: mixed[ig] = beta * f[ig] * out[ig]
    // with the Kerker screen f[ig] = |G+q|^2 / (|G+q|^2 + a^2) built from
    // gcar + q_frac * G (1/lat0^2 units), the v_hartree_q convention
    ModuleDFPT::DFPT_Rho rho_k;
    double w2_min = 0.0;
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + q_cart_;
        const double w2 = w * w;
        if (w2 > 1.0e-12 && (w2_min == 0.0 || w2 < w2_min))
        {
            w2_min = w2;
        }
    }
    ASSERT_GT(w2_min, 0.0);
    const double a2 = 4.0 * w2_min;
    rho_k.init({1, pw_rho_.nrxx, &pw_rho_, &pw_wfc_, G_, std::string("rho_kerker"), 0.7, a2});

    std::vector<std::complex<double>> out(pw_rho_.npw);
    int n_small = 0;
    int n_large = 0;
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        out[ig] = 0.01 * std::complex<double>(std::cos(0.7 * ig), std::sin(0.5 * ig));
        const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + q_cart_;
        const double w2 = w * w;
        if (w2 < a2)
        {
            ++n_small;
        }
        if (w2 > 100.0 * a2)
        {
            ++n_large;
        }
    }
    // sanity: the screen actually varies across the basis
    ASSERT_GT(n_small, 0);
    ASSERT_GT(n_large, 0);

    data_.set_drho_g(0, 0, out);
    rho_k.mix_drho(0, data_);
    const std::vector<std::complex<double>> mixed = data_.get_drho_g(0, 0);
    ASSERT_EQ(mixed.size(), out.size());
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + q_cart_;
        const double w2 = w * w;
        const double f = (w2 < 1.0e-12) ? 0.0 : w2 / (w2 + a2);
        const std::complex<double> ref = 0.7 * f * out[ig];
        EXPECT_NEAR(mixed[ig].real(), ref.real(), 1.0e-12);
        EXPECT_NEAR(mixed[ig].imag(), ref.imag(), 1.0e-12);
    }
    EXPECT_NEAR(rho_k.get_residual(0, data_), 1.0, 1.0e-12);
}

TEST_F(DFPTRhoSerialTest, MixDrhoKerkerStabilizesStiffModelProblem)
{
    // model SCF problem: out(g) = D(g) in(g) + s(g) with the measured
    // diamond-smoke Coulomb-stiffness eigenvalue D = lambda ~ -2.2 on the
    // smallest |G+q| shell and D = 0.3 elsewhere; the fixed point is a
    // fixed target pattern t(g), s = (1 - D) t. Plain mixing must satisfy
    // beta < 2 / (1 + |lambda|) = 0.625, so beta = 0.7 diverges
    // (amplification |1 - beta (1 - D)| = 1.24), while the Kerker screen
    // damps the stiff shell amplification below 1 and converges.
    const int npw = pw_rho_.npw;

    double w2_min = 0.0;
    for (int ig = 0; ig < npw; ++ig)
    {
        const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + q_cart_;
        const double w2 = w * w;
        if (w2 > 1.0e-12 && (w2_min == 0.0 || w2 < w2_min))
        {
            w2_min = w2;
        }
    }
    ASSERT_GT(w2_min, 0.0);

    std::vector<double> stiff(npw, 0.3);
    for (int ig = 0; ig < npw; ++ig)
    {
        const ModuleBase::Vector3<double> w = pw_rho_.gcar[ig] + q_cart_;
        const double w2 = w * w;
        if (w2 > 1.0e-12 && w2 < 1.5 * w2_min)
        {
            stiff[ig] = -2.2;
        }
    }

    std::vector<std::complex<double>> target(npw);
    for (int ig = 0; ig < npw; ++ig)
    {
        target[ig] = 0.01 * std::complex<double>(std::cos(0.3 * ig), std::sin(0.9 * ig));
    }

    auto model_out = [&](const std::vector<std::complex<double>>& in) {
        std::vector<std::complex<double>> o(npw);
        for (int ig = 0; ig < npw; ++ig)
        {
            o[ig] = stiff[ig] * in[ig] + (1.0 - stiff[ig]) * target[ig];
        }
        return o;
    };

    // plain beta = 0.7 on the stiff model diverges
    ModuleDFPT::DFPT_Rho rho_p;
    rho_p.init({1, pw_rho_.nrxx, &pw_rho_, &pw_wfc_, G_, std::string("rho_plain"), 0.7, 0.0});
    data_.set_drho_g(0, 0, std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));
    for (int it = 0; it < 40; ++it)
    {
        data_.set_drho_g(0, 0, model_out(data_.get_drho_g(0, 0)));
        rho_p.mix_drho(0, data_);
    }
    const double residual_plain = rho_p.get_residual(0, data_);
    EXPECT_GT(residual_plain, 1.0);

    // kerker beta = 0.7 with a^2 = 9 w2_min (f ~ 0.1 on the stiff shell)
    // converges to the target
    ModuleDFPT::DFPT_Rho rho_k;
    rho_k.init({1, pw_rho_.nrxx, &pw_rho_, &pw_wfc_, G_, std::string("rho_kerker"), 0.7, 9.0 * w2_min});
    data_.set_drho_g(0, 0, std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));
    for (int it = 0; it < 300; ++it)
    {
        data_.set_drho_g(0, 0, model_out(data_.get_drho_g(0, 0)));
        rho_k.mix_drho(0, data_);
    }
    const double residual_kerker = rho_k.get_residual(0, data_);
    EXPECT_LT(residual_kerker, 1.0e-8);
    const std::vector<std::complex<double>> final_in = data_.get_drho_g(0, 0);
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(final_in[ig].real(), target[ig].real(), 1.0e-10);
        EXPECT_NEAR(final_in[ig].imag(), target[ig].imag(), 1.0e-10);
    }
}

TEST_F(DFPTRhoSerialTest, BroydenConvergesDiagonalSystemFast)
{
    // synthetic fixed point: out = D in + (1 - D) t with a diagonal D
    // (0.5 and 0.3 on two G shells, 0 elsewhere) and a fixed target t,
    // so the response map has a source and the relative residual
    // ||out - in|| / ||out|| genuinely decays. Plain beta = 0.4
    // contracts every populated shell by |1 - beta (1 - D)| = 0.8 (d = 0.5,
    // the binding shell; 0.72 for d = 0.3), i.e. roughly 100 plain
    // iterations to 1e-10; Modified Broyden on the few distinct
    // eigenvalues should converge in a handful of iterations.
    SetupBases(k_d_, q_d_, 2);
    ModuleDFPT::DFPT_Rho::Config cfg;
    cfg.nspin = 1;
    cfg.nrxx = pw_rho_.nrxx;
    cfg.pw_rho = &pw_rho_;
    cfg.pw_wfc = &pw_wfc_;
    cfg.recip_matrix = G_;
    cfg.mix_type = "rho_broyden";
    cfg.mix_beta = 0.4;
    cfg.kerker_a2 = 1.0;
    cfg.broyden_ndim = 8;
    ModuleDFPT::DFPT_Rho rho;
    rho.init(cfg);

    const int npw = pw_rho_.npw;
    std::vector<std::complex<double>> target(npw, std::complex<double>(0.0, 0.0));
    target[0] = std::complex<double>(0.9, -0.2);
    target[1] = std::complex<double>(0.5, 0.2);
    std::vector<double> d(npw, 0.0);
    d[0] = 0.5;
    d[1] = 0.3;
    auto model_out = [&](const std::vector<std::complex<double>>& in) {
        std::vector<std::complex<double>> o(npw);
        for (int ig = 0; ig < npw; ++ig)
        {
            o[ig] = d[ig] * in[ig] + (1.0 - d[ig]) * target[ig];
        }
        return o;
    };

    // run mixing iterations from the zero input; each step the model
    // rebuilds the output from the mixed input
    int n_iter = 0;
    std::vector<std::complex<double>> rho_out = model_out(std::vector<std::complex<double>>(npw, std::complex<double>(0.0, 0.0)));
    for (int iter = 0; iter < 40; ++iter)
    {
        data_.set_drho_g(0, 0, rho_out);
        rho.mix_drho(0, data_);
        ++n_iter;
        const double res = rho.get_residual(0, data_);
        if (res < 1.0e-10)
        {
            break;
        }
        rho_out = model_out(data_.get_drho_g(0, 0));
    }
    // Broyden on a 3-eigenvalue diagonal system: converges dramatically
    // faster than plain mixing (which needs roughly 100 iterations at
    // 0.8/step on the binding shell)
    EXPECT_LT(n_iter, 12);
    const double final_res = rho.get_residual(0, data_);
    EXPECT_LT(final_res, 1.0e-8);
    // the converged input is the fixed point of the model map
    const std::vector<std::complex<double>> final_in = data_.get_drho_g(0, 0);
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(final_in[ig].real(), target[ig].real(), 1.0e-8);
        EXPECT_NEAR(final_in[ig].imag(), target[ig].imag(), 1.0e-8);
    }
}

TEST_F(DFPTRhoSerialTest, BroydenResetMixingAllowsFreshQAndSameQCycle)
{
    // regression (PR review): entering a fresh q after q = 0 must reset
    // the shared Broyden mixer's address binding even though the new q
    // has no history entry yet -- otherwise the first push of the newly
    // allocated Mixing_Data exits through the one-object-one-data guard
    // (broyden_mixing.cpp) and kills the process. The same-q cycle must
    // also restart from a zero input density.
    SetupBases(k_d_, q_d_, 2);
    ModuleDFPT::DFPT_Rho::Config cfg;
    cfg.nspin = 1;
    cfg.nrxx = pw_rho_.nrxx;
    cfg.pw_rho = &pw_rho_;
    cfg.pw_wfc = &pw_wfc_;
    cfg.recip_matrix = G_;
    cfg.mix_type = "rho_broyden";
    cfg.mix_beta = 0.4;
    cfg.kerker_a2 = 1.0;
    cfg.broyden_ndim = 8;
    ModuleDFPT::DFPT_Rho rho;
    rho.init(cfg);

    const int npw = pw_rho_.npw;
    std::vector<std::complex<double>> out(npw);
    for (int ig = 0; ig < npw; ++ig)
    {
        out[ig] = 0.01 * std::complex<double>(std::cos(0.3 * ig), std::sin(0.9 * ig));
    }

    // two pushes at q = 0: the mixer binds the q = 0 Mixing_Data
    data_.set_drho_g(0, 0, out);
    rho.mix_drho(0, data_);
    data_.set_drho_g(0, 0, out);
    rho.mix_drho(0, data_);

    // entering q = 1 for the first time (no history entry exists): with
    // the buggy guard the mixer stayed bound to q = 0 and this died
    rho.reset_mixing(1);
    data_.set_drho_g(1, 0, out);
    rho.mix_drho(1, data_);
    // fresh q starts from the zero input: first step is beta * out
    const std::vector<std::complex<double>> mixed_q1 = data_.get_drho_g(1, 0);
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(mixed_q1[ig].real(), 0.4 * out[ig].real(), 1.0e-12);
        EXPECT_NEAR(mixed_q1[ig].imag(), 0.4 * out[ig].imag(), 1.0e-12);
    }

    // the same-q cycle also restarts from a zero input density
    rho.reset_mixing(0);
    data_.set_drho_g(0, 0, out);
    rho.mix_drho(0, data_);
    const std::vector<std::complex<double>> mixed_q0 = data_.get_drho_g(0, 0);
    for (int ig = 0; ig < npw; ++ig)
    {
        EXPECT_NEAR(mixed_q0[ig].real(), 0.4 * out[ig].real(), 1.0e-12);
        EXPECT_NEAR(mixed_q0[ig].imag(), 0.4 * out[ig].imag(), 1.0e-12);
    }
}

// ---------------------------------------------------------------------------
// nlcc phase regression: the first-order core response must localize at the
// ATOM site, not its mirror. This is the anti-blind-spot test for the
// structure-factor sign: exp(+i(G+q)tau) mirrors the response to -tau, a bug
// every historical validation missed because the nlcc atom sat at the origin
// (bn_cubic) where both conventions coincide. The wurtzite bn_hex comparison
// against QE dyn.G pinned the sign: the mirror of B1 lands exactly on an N
// site of the other layer, collapsing the interlayer shear to -399 cm^-1.
// ---------------------------------------------------------------------------
TEST_F(DFPTRhoSerialTest, DrhocLocalizesAtTheAtomNotItsMirror)
{
    MakeNlccTwoAtomCell();
    SetupBases(k_d_, q_d_, 2);
    ModuleDFPT::DFPT_Rho::Config cfg;
    cfg.nspin = 1;
    cfg.nrxx = pw_rho_.nrxx;
    cfg.pw_rho = &pw_rho_;
    cfg.pw_wfc = &pw_wfc_;
    cfg.recip_matrix = G_;
    cfg.mix_type = "rho_plain";
    cfg.mix_beta = 0.4;
    cfg.kerker_a2 = 1.0;
    cfg.broyden_ndim = 8;
    ModuleDFPT::DFPT_Rho rho;
    rho.init(cfg);

    // atom 1 carries the nlcc at the non-origin tau_nlcc_; all three
    // displacement directions must localize at the atom
    const int iat = 1;
    const int it = ucell_.iat2it[iat];
    ASSERT_TRUE(ucell_.atoms[it].ncpp.nlcc);
    // tau is stored in lat0-cartesian units (fixture convention: cell edge
    // a_ in lat0 units); the true fractional position divides by a_
    const double inv_a = 1.0 / a_;
    const ModuleBase::Vector3<double> tau_frac(ucell_.atoms[it].tau[ucell_.iat2ia[iat]].x * inv_a,
                                               ucell_.atoms[it].tau[ucell_.iat2ia[iat]].y * inv_a,
                                               ucell_.atoms[it].tau[ucell_.iat2ia[iat]].z * inv_a);
    // mirror position folded back into the cell
    ModuleBase::Vector3<double> mirror_frac(1.0 - tau_frac.x, 1.0 - tau_frac.y, 1.0 - tau_frac.z);

    const int nrxx = pw_rho_.nrxx;
    std::vector<std::complex<double>> drhoc_g;
    std::vector<std::complex<double>> drhoc_r(nrxx);
    for (int dir = 0; dir < 3; ++dir)
    {
        // q = 0 first: the amplitude is the plain derivative profile, its
        // |.|^2 centroid is at the atom (up to the derivative skew)
        const ModuleBase::Vector3<double> q_zero(0.0, 0.0, 0.0);
        rho.compute_drhoc_g(ucell_, iat, dir, q_zero, drhoc_g);
        ASSERT_EQ(static_cast<int>(drhoc_g.size()), pw_rho_.npw);
        pw_rho_.recip2real(drhoc_g.data(), drhoc_r.data());
        double wsum = 0.0;
        double cx = 0.0;
        double cy = 0.0;
        double cz = 0.0;
        for (int ir = 0; ir < nrxx; ++ir)
        {
            const double wgt = std::norm(drhoc_r[ir]);
            // fractional coordinates of the grid point via the FFT triple
            const int iz = ir % pw_rho_.nz;
            const int iy = (ir / pw_rho_.nz) % pw_rho_.ny;
            const int ix = ir / (pw_rho_.ny * pw_rho_.nz);
            // wrap-around: grid index -> fractional in [0, 1)
            const double fx = static_cast<double>(ix) / pw_rho_.nx;
            const double fy = static_cast<double>(iy) / pw_rho_.ny;
            const double fz = static_cast<double>(iz) / pw_rho_.nz;
            // weight by |field|^2, folding the +1 wrap for the periodic image
            // closest to tau: the field is localized, so the unwrapped
            // centroid near tau is what we need
            double wx = fx;
            double wy = fy;
            double wz = fz;
            if (std::abs(fx - tau_frac.x) > 0.5)
            {
                wx += (fx < tau_frac.x) ? 1.0 : -1.0;
            }
            if (std::abs(fy - tau_frac.y) > 0.5)
            {
                wy += (fy < tau_frac.y) ? 1.0 : -1.0;
            }
            if (std::abs(fz - tau_frac.z) > 0.5)
            {
                wz += (fz < tau_frac.z) ? 1.0 : -1.0;
            }
            wsum += wgt;
            cx += wgt * wx;
            cy += wgt * wy;
            cz += wgt * wz;
        }
        ASSERT_GT(wsum, 0.0);
        const ModuleBase::Vector3<double> com(cx / wsum, cy / wsum, cz / wsum);
        // minimum-image distance: wrap each component into [-0.5, 0.5)
        auto mindist2 = [](const ModuleBase::Vector3<double>& a,
                           const ModuleBase::Vector3<double>& b) {
            double dx = a.x - b.x; dx -= std::round(dx);
            double dy = a.y - b.y; dy -= std::round(dy);
            double dz = a.z - b.z; dz -= std::round(dz);
            return dx * dx + dy * dy + dz * dz;
        };
        const double dist_atom = mindist2(com, tau_frac);
        const double dist_mirror = mindist2(com, mirror_frac);
        // the centroid of a derivative field sits at the atom up to the
        // derivative skew; the mirror is a full half-cell away in the
        // minimum-image sense, far beyond any tolerance that would admit
        // the reflected profile
        EXPECT_LT(dist_atom, 0.01) << "dir=" << dir;
        EXPECT_GT(dist_mirror, 25.0 * dist_atom) << "dir=" << dir;
    }

    // non-zero q variant: the q-phase exp(-i q.tau) must not mirror either.
    // Use a generic incommensurate q (the fixture default q_d_ is already
    // generic); the real-space field of the q-shifted amplitude still
    // localizes at the atom (the e^{iqr} carrier is smooth on the sigma
    // scale of the core).
    const ModuleBase::Vector3<double> q2 = q_d_ * 0.7 + ModuleBase::Vector3<double>(0.11, -0.07, 0.05);
    const ModuleBase::Vector3<double> q2_cart = q2 * ucell_.G;
    rho.compute_drhoc_g(ucell_, iat, 1, q2_cart, drhoc_g);
    ASSERT_EQ(static_cast<int>(drhoc_g.size()), pw_rho_.npw);
    pw_rho_.recip2real(drhoc_g.data(), drhoc_r.data());
    double wsum = 0.0;
    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;
    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double wgt = std::norm(drhoc_r[ir]);
        const int iz = ir % pw_rho_.nz;
        const int iy = (ir / pw_rho_.nz) % pw_rho_.ny;
        const int ix = ir / (pw_rho_.ny * pw_rho_.nz);
        double wx = static_cast<double>(ix) / pw_rho_.nx;
        double wy = static_cast<double>(iy) / pw_rho_.ny;
        double wz = static_cast<double>(iz) / pw_rho_.nz;
        if (std::abs(wx - tau_frac.x) > 0.5)
        {
            wx += (wx < tau_frac.x) ? 1.0 : -1.0;
        }
        if (std::abs(wy - tau_frac.y) > 0.5)
        {
            wy += (wy < tau_frac.y) ? 1.0 : -1.0;
        }
        if (std::abs(wz - tau_frac.z) > 0.5)
        {
            wz += (wz < tau_frac.z) ? 1.0 : -1.0;
        }
        wsum += wgt;
        cx += wgt * wx;
        cy += wgt * wy;
        cz += wgt * wz;
    }
    ASSERT_GT(wsum, 0.0);
    const ModuleBase::Vector3<double> com2(cx / wsum, cy / wsum, cz / wsum);
    // same minimum-image comparison as above
    auto mindist2 = [](const ModuleBase::Vector3<double>& a,
                       const ModuleBase::Vector3<double>& b) {
        double dx = a.x - b.x; dx -= std::round(dx);
        double dy = a.y - b.y; dy -= std::round(dy);
        double dz = a.z - b.z; dz -= std::round(dz);
        return dx * dx + dy * dy + dz * dz;
    };
    const double d2a = mindist2(com2, tau_frac);
    const double d2m = mindist2(com2, mirror_frac);
    EXPECT_LT(d2a, 0.02);
    EXPECT_GT(d2m, 25.0 * d2a);

    // phase-sensitive check at generic q: the centroid via |field|^2 is
    // invariant under a global phase rotation, so flipping the sign of
    // exp(i q.tau) alone would leave it unchanged. Instead, verify the
    // G-space coefficient phase directly: drhoc_g(G) must carry
    // exp(-i(G+q).tau) (the QE addcore convention), not the conjugate.
    // Compare the q=0 and q=q2 coefficients: their ratio at each G must
    // be exp(-i q.cart . tau) * (ratio of the radial/Bessel factor at
    // |G| vs |G+q|). For a single G with a dominant coefficient this
    // isolates the structure-factor phase.
    const ModuleBase::Vector3<double> q_zero_cart(0.0, 0.0, 0.0);
    std::vector<std::complex<double>> drhoc_g0;
    rho.compute_drhoc_g(ucell_, iat, 1, q_zero_cart, drhoc_g0);
    ASSERT_EQ(static_cast<int>(drhoc_g0.size()), pw_rho_.npw);
    // find the G with the largest |G|=0-adjacent coefficient (where the
    // radial factor is comparable) and verify the phase relation
    int ig_best = -1;
    double best_mag = 0.0;
    for (int ig = 0; ig < pw_rho_.npw; ++ig)
    {
        const double m = std::norm(drhoc_g0[ig]);
        if (m > best_mag)
        {
            best_mag = m;
            ig_best = ig;
        }
    }
    ASSERT_GE(ig_best, 0);
    // the expected extra phase from q is exp(-i q_cart . tau_cart):
    // drhoc_q(G) = -i|G+q| tpiba rho~(|G+q|) exp(-i(G+q).tau)
    // drhoc_0(G) = -i|G|   tpiba rho~(|G|)   exp(-i G. tau)
    // ratio = [|G+q|/|G| · rho~(|G+q|)/rho~(|G|)] · exp(-i q. tau)
    // The bracket is positive real; the phase is exactly -q_cart . tau_cart.
    // Use the full TWO_PI-scaled phase (the tau and q_cart are already in
    // consistent units: tau in lat0-cartesian, q_cart in 2pi/lat0 units,
    // so the dot product directly gives the phase in radians).
    const std::complex<double> ratio = drhoc_g[ig_best] / drhoc_g0[ig_best];
    const ModuleBase::Vector3<double> tau_cart_alat = ucell_.atoms[it].tau[ucell_.iat2ia[iat]];
    const double expected_phase_arg = -ModuleBase::TWO_PI * (q2_cart * tau_cart_alat) / ucell_.lat0;
    const double cos_arg = std::cos(expected_phase_arg);
    const double sin_arg = std::sin(expected_phase_arg);
    const std::complex<double> expected_phase(cos_arg, sin_arg);
    // tight tolerance: the exact phase relation is exp(-i q.tau), so both
    // omission (arg=0) and reversal (arg=+q.tau) must be rejected; the
    // only residual is the positive-real radial-ratio bracket (which has
    // zero phase contribution) and floating-point noise
    // rule 14: named locals for computed function arguments
    const std::complex<double> unit_phase_ratio = ratio / expected_phase;
    const double phase_diff = std::arg(unit_phase_ratio);
    const double sin_diff = std::sin(phase_diff);
    const double cos_diff = std::cos(phase_diff);
    const double phase_diff_wrapped = std::atan2(sin_diff, cos_diff);
    EXPECT_NEAR(phase_diff_wrapped, 0.0, 0.05) << "phase of drhoc_q/drhoc_0 ratio vs expected exp(-i q.tau)";
}
