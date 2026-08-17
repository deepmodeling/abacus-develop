#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <functional>
#include <vector>

// serial unit test of the first-order density response (C3).
// Everything here runs without __MPI: the plane-wave bases are built through
// the real serial initgrids/initparameters/setuptransform path on a shared
// FFT grid, exactly like the production setup_pwrho/setup_pwwfc sequence.

#define private public
#include "source_cell/qlist.h"
#undef private

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

namespace {

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

class DFPTRhoSerialTest : public testing::Test
{
  protected:
    const double lat0_ = 1.8897261254578281;
    const double ecutwfc_ = 2.5; // Ry
    const double rho_mult_ = 9.0;

    ModuleBase::Matrix3 latvec_;
    ModulePW::PW_Basis pw_rho_;
    ModulePW::PW_Basis_K pw_wfc_;
    ModuleCell::QList qlist_;
    ModuleDFPT::DFPT_PW_Data data_;
    ModuleDFPT::DFPT_Rho rho_;

    const ModuleBase::Vector3<double> q_d_{0.13, 0.0, 0.07};
    const ModuleBase::Vector3<double> k_d_{-0.13, 0.0, -0.07};
    ModuleBase::Vector3<double> q_cart_;
    ModuleBase::Matrix3 G_;

    static const int nbands_ = 2;

    void SetUp() override
    {
        latvec_ = ModuleBase::Matrix3(10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0);
        G_ = latvec_.Inverse().Transpose();

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
        q_cart_ = q_d_ * G_;

        data_.init(&qlist_, 1, nbands_, pw_wfc_.npwk_max, pw_rho_.nrxx, 1, 1, nullptr);
        rho_.init(1, pw_rho_.nrxx, &pw_rho_, &pw_wfc_, G_, "plain", 0.4);
    }

    void FillRandomStates(psi::Psi<std::complex<double>>& psi,
                          std::vector<std::vector<std::vector<std::complex<double>>>>& dpsi)
    {
        psi::Psi<std::complex<double>> p(1, nbands_, pw_wfc_.npwk_max, pw_wfc_.npwk[0], true);
        ModuleDFPT::DFPT_KQ_Basis kq;
        kq.init(&pw_wfc_, q_cart_, 0);
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
    kq.init(&pw_wfc_, q_cart_, 0);
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
        const ModuleBase::Vector3<double> delta =
            ModuleBase::Vector3<double>(mx, my, mz) * G_;
        // A_Delta = sum_G c*_G d_{G+Delta}, brute-forced over the two lists
        std::complex<double> aref(0.0, 0.0);
        for (int jgl = 0; jgl < kq.get_npwk(); ++jgl)
        {
            const ModuleBase::Vector3<double> gq = kq.get_gcar(jgl);
            for (size_t j = 0; j < glist.size(); ++j)
            {
                if (std::abs(glist[j].x - (gq.x - delta.x)) < 1.0e-6 &&
                    std::abs(glist[j].y - (gq.y - delta.y)) < 1.0e-6 &&
                    std::abs(glist[j].z - (gq.z - delta.z)) < 1.0e-6)
                {
                    aref += std::conj(clist[j]) * dvec[jgl];
                    break;
                }
            }
        }
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
    kq.init(&pw_wfc_, q_cart_, 0);
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
        const ModuleBase::Vector3<double> r_cart =
            ModuleBase::Vector3<double>(fx, fy, fz) * latvec_;
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
        const double ref = 2.0 * (std::conj(u) * du * eq).real();
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
    rho0.init(1, pw_rho_.nrxx, &pw_rho_, &pw_wfc0, G_, "plain", 0.4);

    psi::Psi<std::complex<double>> psi(1, nbands_, pw_wfc0.npwk_max, pw_wfc0.npwk[0], true);
    ModuleDFPT::DFPT_KQ_Basis kq0;
    kq0.init(&pw_wfc0, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0);
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
    const int ig_star = [this]()
    {
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
        = ModuleBase::e2 * ModuleBase::FOUR_PI / (pw_rho_.tpiba2 * (w * w))
          * drho_g[ig_star];
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
