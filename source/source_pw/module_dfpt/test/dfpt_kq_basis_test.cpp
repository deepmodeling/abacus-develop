#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_pw/module_dfpt/dfpt_kq_basis.h"

/************************************************
 *  unit test of DFPT_KQ_Basis (C0)
 ***********************************************/

/**
 * - Tested Functions:
 *   - DFPT_KQ_Basis::init() - enumeration of the local k+q plane-wave
 *     basis from the union of the ground-state wavefunction G grid and the
 *     denser charge G grid (the wavefunction reservoir alone does not
 *     cover q-shifted balls: its radius is (sqrt(gk_ecut)+max|k_mesh|)^2).
 *   - Accessors get_npwk / get_ig_rho / get_gcar / get_gpluskq / get_gk2 /
 *     get_kplusq.
 *
 * Both bases are hand-built with public members only (no FFT setup
 * needed): complex (gamma_only=false) grids on a cubic lattice. Every
 * selection is cross-checked against an independent brute-force count
 * over the full FFT grid. The TruncatedWfcReservoirCompletedByRhoGrid
 * case reproduces the q!=0 truncation defect the dual-reservoir
 * enumeration fixes.
 */

namespace {

bool VecLess(const ModuleBase::Vector3<double>& a, const ModuleBase::Vector3<double>& b)
{
    if (a.x != b.x)
    {
        return a.x < b.x;
    }
    if (a.y != b.y)
    {
        return a.y < b.y;
    }
    return a.z < b.z;
}

// mirror the wrap used by cal_GplusK_cartesian / collect_local_pw
int WrapIndex(int i, int n)
{
    if (i >= n / 2 + 1)
    {
        return i - n;
    }
    return i;
}

class DFPTKQBasisTest : public testing::Test
{
  protected:
    ModulePW::PW_Basis_K pw_;
    ModulePW::PW_Basis prho_;
    const double lat0_ = 1.8897261254578281;
    const double ecutwfc_ = 520.0; // ball radius ~ 2 G shells
    double tpiba2_ = 0.0;
    double gk_ecut_ = 0.0;
    double ggecut_ = 0.0;
    ModuleBase::Matrix3 G_;
    const int nx_ = 7, ny_ = 7, nz_ = 7;

    void ResetBasis()
    {
        delete[] pw_.kvec_c;
        delete[] pw_.ig2isz;
        delete[] pw_.is2fftixy;
        pw_ = ModulePW::PW_Basis_K();
        pw_.nx = 0;
        pw_.ny = 0;
        pw_.nz = 0;
        pw_.fftny = 0;
        pw_.npw = 0;
        pw_.nst = 0;
        delete[] prho_.ig2isz;
        delete[] prho_.is2fftixy;
        prho_ = ModulePW::PW_Basis();
        prho_.nx = 0;
        prho_.ny = 0;
        prho_.nz = 0;
        prho_.fftny = 0;
        prho_.npw = 0;
        prho_.nst = 0;
    }

    // fill the (stick, z) layout the real distribution code produces for a
    // given ball; fftny = ny (complex basis). nxyz must be preset.
    void FillSticks(ModulePW::PW_Basis& pw, double ggecut)
    {
        std::vector<int> is2fftixy;
        std::vector<int> ig2isz;
        for (int ix0 = 0; ix0 < pw.nx; ++ix0)
        {
            const int wix = WrapIndex(ix0, pw.nx);
            for (int iy0 = 0; iy0 < pw.ny; ++iy0)
            {
                const int wiy = WrapIndex(iy0, pw.ny);
                std::vector<int> stick;
                for (int iz0 = 0; iz0 < pw.nz; ++iz0)
                {
                    const int wiz = WrapIndex(iz0, pw.nz);
                    ModuleBase::Vector3<double> f(wix, wiy, wiz);
                    const ModuleBase::Vector3<double> g = f * G_;
                    if (g * g <= ggecut)
                    {
                        stick.push_back(iz0);
                    }
                }
                if (!stick.empty())
                {
                    const int is = static_cast<int>(is2fftixy.size());
                    is2fftixy.push_back(iy0 + ix0 * pw.fftny);
                    for (size_t s = 0; s < stick.size(); ++s)
                    {
                        ig2isz.push_back(is * pw.nz + stick[s]);
                    }
                }
            }
        }
        pw.nst = static_cast<int>(is2fftixy.size());
        pw.npw = static_cast<int>(ig2isz.size());
        pw.is2fftixy = new int[pw.nst];
        pw.ig2isz = new int[pw.npw];
        for (int i = 0; i < pw.nst; ++i)
        {
            pw.is2fftixy[i] = is2fftixy[i];
        }
        for (int i = 0; i < pw.npw; ++i)
        {
            pw.ig2isz[i] = ig2isz[i];
        }
    }

    // shared field setup for a cubic complex basis. The wavefunction grid
    // uses the production radius (sqrt(gk_ecut) + max|k_mesh|)^2; the
    // charge grid uses the production 4*ecutwfc (= 4*gk_ecut ball, no
    // k list) that covers every q-shifted ball for q in the first BZ.
    void BuildBase(const std::vector<ModuleBase::Vector3<double>>& kvec_c)
    {
        tpiba2_ = ModuleBase::TWO_PI * ModuleBase::TWO_PI / (lat0_ * lat0_);
        gk_ecut_ = ecutwfc_ / tpiba2_;
        const double b = ModuleBase::TWO_PI / lat0_;
        G_.e11 = b;
        G_.e12 = 0;
        G_.e13 = 0;
        G_.e21 = 0;
        G_.e22 = b;
        G_.e23 = 0;
        G_.e31 = 0;
        G_.e32 = 0;
        G_.e33 = b;

        double kmaxmod = 0.0;
        for (size_t i = 0; i < kvec_c.size(); ++i)
        {
            kmaxmod = std::max(kmaxmod, std::sqrt(kvec_c[i] * kvec_c[i]));
        }
        ggecut_ = std::pow(std::sqrt(gk_ecut_) + kmaxmod, 2);

        pw_.nx = nx_;
        pw_.ny = ny_;
        pw_.nz = nz_;
        pw_.nxyz = nx_ * ny_ * nz_;
        pw_.fftny = ny_; // gamma_only = false
        pw_.gamma_only = false;
        pw_.G = G_;
        pw_.ggecut = ggecut_;
        pw_.gk_ecut = gk_ecut_;
        pw_.nks = static_cast<int>(kvec_c.size());
        pw_.kvec_c = new ModuleBase::Vector3<double>[pw_.nks];
        for (int i = 0; i < pw_.nks; ++i)
        {
            pw_.kvec_c[i] = kvec_c[i];
        }
        FillSticks(pw_, ggecut_);

        prho_.nx = nx_;
        prho_.ny = ny_;
        prho_.nz = nz_;
        prho_.nxyz = nx_ * ny_ * nz_;
        prho_.fftny = ny_;
        prho_.gamma_only = false;
        prho_.G = G_;
        prho_.ggecut = 4.0 * gk_ecut_;
        FillSticks(prho_, prho_.ggecut);
    }

    // independent reference: brute-force count/collect over the whole FFT
    // grid for a given shifted center, sorted for set comparison.
    std::vector<ModuleBase::Vector3<double>> ReferenceSelection(const ModuleBase::Vector3<double>& center)
    {
        std::vector<ModuleBase::Vector3<double>> out;
        for (int ix0 = 0; ix0 < nx_; ++ix0)
        {
            const int wix = WrapIndex(ix0, nx_);
            for (int iy0 = 0; iy0 < ny_; ++iy0)
            {
                const int wiy = WrapIndex(iy0, ny_);
                for (int iz0 = 0; iz0 < nz_; ++iz0)
                {
                    const int wiz = WrapIndex(iz0, nz_);
                    ModuleBase::Vector3<double> f(wix, wiy, wiz);
                    ModuleBase::Vector3<double> g = f * G_;
                    ModuleBase::Vector3<double> gp = g + center;
                    if (gp * gp <= gk_ecut_)
                    {
                        out.push_back(gp);
                    }
                }
            }
        }
        std::sort(out.begin(), out.end(), VecLess);
        return out;
    }

    std::vector<ModuleBase::Vector3<double>> KqSet(const ModuleDFPT::DFPT_KQ_Basis& kq)
    {
        std::vector<ModuleBase::Vector3<double>> out = kq.get_gcar_all();
        for (size_t i = 0; i < out.size(); ++i)
        {
            out[i] = kq.get_gpluskq(i);
        }
        std::sort(out.begin(), out.end(), VecLess);
        return out;
    }
};

TEST_F(DFPTKQBasisTest, GammaQ0ReproducesWfcGrid)
{
    // single Gamma k: kmaxmod = 0, so the wavefunction grid ball equals the
    // Gamma ball and the q=0 selection must reproduce it verbatim.
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0)});

    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_, &prho_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0);
    ASSERT_TRUE(kq.is_valid());
    EXPECT_EQ(kq.get_npwk(), pw_.npw);

    // every selected vector lies inside the cutoff and on the brute-force set
    const std::vector<ModuleBase::Vector3<double>> ref = ReferenceSelection(
        ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    EXPECT_EQ(static_cast<int>(ref.size()), pw_.npw);
    const std::vector<ModuleBase::Vector3<double>> sel = KqSet(kq);
    EXPECT_EQ(sel.size(), ref.size());
    for (size_t i = 0; i < ref.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(sel[i].x, ref[i].x);
        EXPECT_DOUBLE_EQ(sel[i].y, ref[i].y);
        EXPECT_DOUBLE_EQ(sel[i].z, ref[i].z);
    }

    // rho-grid indices are unique and in range
    std::vector<int> igs;
    for (int igl = 0; igl < kq.get_npwk(); ++igl)
    {
        const int ig = kq.get_ig_rho(igl);
        EXPECT_GE(ig, 0);
        EXPECT_LT(ig, prho_.npw);
        igs.push_back(ig);
    }
    std::sort(igs.begin(), igs.end());
    for (size_t i = 1; i < igs.size(); ++i)
    {
        EXPECT_NE(igs[i], igs[i - 1]);
    }
}

TEST_F(DFPTKQBasisTest, ShiftedCenterSelectsAsymmetricSphere)
{
    // k = (0,0,0.5b): the |G+k|^2 cut keeps an asymmetric shell
    const double b = ModuleBase::TWO_PI / lat0_;
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
               ModuleBase::Vector3<double>(0.0, 0.0, 0.5 * b)});

    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_, &prho_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 1);
    const ModuleBase::Vector3<double> center = kq.get_kplusq();
    EXPECT_NEAR(center.x, 0.0, 1e-12);
    EXPECT_NEAR(center.y, 0.0, 1e-12);
    EXPECT_NEAR(center.z, 0.5 * b, 1e-12);

    // independent brute force on the full FFT grid
    const std::vector<ModuleBase::Vector3<double>> ref = ReferenceSelection(center);
    const std::vector<ModuleBase::Vector3<double>> sel = KqSet(kq);
    EXPECT_EQ(sel.size(), ref.size());
    for (size_t i = 0; i < ref.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(sel[i].x, ref[i].x);
        EXPECT_DOUBLE_EQ(sel[i].y, ref[i].y);
        EXPECT_DOUBLE_EQ(sel[i].z, ref[i].z);
    }
}

TEST_F(DFPTKQBasisTest, TruncatedWfcReservoirCompletedByRhoGrid)
{
    // regression of the q!=0 truncation defect: the wavefunction grid is
    // sized for the k-mesh only, so a q-shifted ball can poke outside it
    // while staying inside the 4*ecutwfc charge ball. The dual-reservoir
    // enumeration must still return the full brute-force set.
    const double b = ModuleBase::TWO_PI / lat0_;
    const ModuleBase::Vector3<double> k1(0.0, 0.0, 0.5 * b);
    const ModuleBase::Vector3<double> q(0.5 * b, 0.5 * b, 0.5 * b);
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0), k1});

    // the k+q ball must reach beyond the wavefunction grid ball, otherwise
    // this test would not exercise the completing reservoir
    const double need = std::sqrt((k1 + q) * (k1 + q)) + std::sqrt(gk_ecut_);
    EXPECT_GT(need, std::sqrt(ggecut_));
    EXPECT_LE(need, std::sqrt(prho_.ggecut));

    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_, &prho_, q, 1);
    const std::vector<ModuleBase::Vector3<double>> ref = ReferenceSelection(k1 + q);
    const std::vector<ModuleBase::Vector3<double>> sel = KqSet(kq);
    ASSERT_EQ(sel.size(), ref.size());
    for (size_t i = 0; i < ref.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(sel[i].x, ref[i].x);
        EXPECT_DOUBLE_EQ(sel[i].y, ref[i].y);
        EXPECT_DOUBLE_EQ(sel[i].z, ref[i].z);
    }
    for (int igl = 0; igl < kq.get_npwk(); ++igl)
    {
        EXPECT_LE(kq.get_gk2(igl), gk_ecut_ + 1e-12);
        EXPECT_GE(kq.get_ig_rho(igl), 0);
    }

    // cross-check every returned rho-grid index: the rho G it points to
    // must carry the same integer triplet as the k+q entry
    for (int igl = 0; igl < kq.get_npwk(); ++igl)
    {
        const int ig = kq.get_ig_rho(igl);
        ASSERT_GE(ig, 0);
        const int isz = prho_.ig2isz[ig];
        const int iz = WrapIndex(isz % prho_.nz, prho_.nz);
        const int ixy = prho_.is2fftixy[isz / prho_.nz];
        const int ix = WrapIndex(ixy / prho_.fftny, prho_.nx);
        const int iy = WrapIndex(ixy % prho_.fftny, prho_.ny);
        const ModuleBase::Vector3<double> want = kq.get_gcar(igl);
        EXPECT_NEAR(ix * b, want.x, 1e-10);
        EXPECT_NEAR(iy * b, want.y, 1e-10);
        EXPECT_NEAR(iz * b, want.z, 1e-10);
    }
}

TEST_F(DFPTKQBasisTest, TranslationInvarianceOfKQ)
{
    // the k+q basis depends only on the sum k+q: (ik=0, q=k1) must agree
    // with (ik=1, q=0), both centered at k1
    const ModuleBase::Vector3<double> k1(0.0, 0.0, 0.5 * ModuleBase::TWO_PI / lat0_);
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0), k1});

    ModuleDFPT::DFPT_KQ_Basis a, b;
    a.init(&pw_, &prho_, k1, 0);
    b.init(&pw_, &prho_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 1);
    ASSERT_EQ(a.get_npwk(), b.get_npwk());
    for (int igl = 0; igl < a.get_npwk(); ++igl)
    {
        EXPECT_EQ(a.get_ig_rho(igl), b.get_ig_rho(igl));
        EXPECT_DOUBLE_EQ(a.get_gk2(igl), b.get_gk2(igl));
        EXPECT_DOUBLE_EQ(a.get_gpluskq(igl).x, b.get_gpluskq(igl).x);
        EXPECT_DOUBLE_EQ(a.get_gpluskq(igl).y, b.get_gpluskq(igl).y);
        EXPECT_DOUBLE_EQ(a.get_gpluskq(igl).z, b.get_gpluskq(igl).z);
    }

    // shifting back by -k1 recovers the Gamma basis
    ModuleDFPT::DFPT_KQ_Basis c;
    c.init(&pw_, &prho_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0) - k1, 1);
    for (int igl = 0; igl < c.get_npwk(); ++igl)
    {
        EXPECT_DOUBLE_EQ(c.get_gk2(igl), c.get_gcar(igl) * c.get_gcar(igl));
    }
}

TEST_F(DFPTKQBasisTest, InvalidOrMismatchedBaseIsRejected)
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    // null providers: valid-but-empty basis
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(nullptr, nullptr, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0);
    EXPECT_FALSE(kq.is_valid());
    EXPECT_EQ(kq.get_npwk(), 0);

    // gamma_only base is rejected (DFPT needs the full complex G ball)
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0)});
    pw_.gamma_only = true;
    pw_.fftny = ny_ / 2 + 1;
    ModuleDFPT::DFPT_KQ_Basis kq2;
    EXPECT_EXIT(kq2.init(&pw_, &prho_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0),
                ::testing::ExitedWithCode(1),
                "");

    // mismatched FFT grid dimensions between the two bases are rejected
    pw_.gamma_only = false;
    pw_.fftny = ny_;
    prho_.nx = 9;
    prho_.nxyz = 9 * prho_.ny * prho_.nz;
    ModuleDFPT::DFPT_KQ_Basis kq3;
    EXPECT_EXIT(kq3.init(&pw_, &prho_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0),
                ::testing::ExitedWithCode(1),
                "");
}

} // namespace
