#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_pw/module_dfpt/dfpt_kq_basis.h"

/************************************************
 *  unit test of DFPT_KQ_Basis (C0)
 ***********************************************/

/**
 * - Tested Functions:
 *   - DFPT_KQ_Basis::init() - enumeration of the local k+q plane-wave
 *     basis from an initialized ground-state k-basis by re-filtering the
 *     shared G grid at the shifted center k+q.
 *   - Accessors get_npwk / get_ig / get_ig2isz / get_gcar /
 *     get_gpluskq / get_gk2 / get_kplusq.
 *
 * The ground-state k-basis is hand-built with public members only (no FFT
 * setup needed): a complex (gamma_only=false) basis on a cubic lattice with
 * an ecutwfc ball large enough to contain several G shells. Every selection
 * is cross-checked against an independent brute-force count over the full
 * FFT grid, and the k+q -> k+q translation invariance is verified.
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
    const double lat0_ = 1.8897261254578281;
    const double ecutwfc_ = 130.0; // gamma ball reaches the first G shell
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
    }

    // shared field setup for a cubic complex basis; builds the FFT-grid G
    // set with |G|^2 <= ggecut and fills ig2isz / is2fftixy in the same
    // (stick, z) layout the real distribution code produces.
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

        // stick layout: one stick per (ix, iy) pair that has at least one
        // qualifying z plane; is2fftixy[is] = iy + ix * fftny.
        std::vector<int> is2fftixy;
        std::vector<int> ig2isz;
        for (int ix0 = 0; ix0 < nx_; ++ix0)
        {
            const int wix = WrapIndex(ix0, nx_);
            for (int iy0 = 0; iy0 < ny_; ++iy0)
            {
                const int wiy = WrapIndex(iy0, ny_);
                std::vector<int> stick;
                for (int iz0 = 0; iz0 < nz_; ++iz0)
                {
                    const int wiz = WrapIndex(iz0, nz_);
                    ModuleBase::Vector3<double> f(wix, wiy, wiz);
                    const ModuleBase::Vector3<double> g = f * G_;
                    if (g * g <= ggecut_)
                    {
                        stick.push_back(iz0);
                    }
                }
                if (!stick.empty())
                {
                    const int is = static_cast<int>(is2fftixy.size());
                    is2fftixy.push_back(iy0 + ix0 * pw_.fftny);
                    for (size_t s = 0; s < stick.size(); ++s)
                    {
                        ig2isz.push_back(is * nz_ + stick[s]);
                    }
                }
            }
        }
        pw_.nst = static_cast<int>(is2fftixy.size());
        pw_.npw = static_cast<int>(ig2isz.size());
        pw_.is2fftixy = new int[pw_.nst];
        pw_.ig2isz = new int[pw_.npw];
        for (int i = 0; i < pw_.nst; ++i)
        {
            pw_.is2fftixy[i] = is2fftixy[i];
        }
        for (int i = 0; i < pw_.npw; ++i)
        {
            pw_.ig2isz[i] = ig2isz[i];
        }
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

TEST_F(DFPTKQBasisTest, GammaQ0ReproducesBaseOrdering)
{
    // single Gamma k: kmaxmod = 0, so the shared grid ball equals the Gamma
    // ball and the q=0 selection must reproduce the base basis verbatim.
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0)});
    ASSERT_EQ(pw_.npw, 7);

    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0);
    ASSERT_TRUE(kq.is_valid());
    EXPECT_EQ(kq.get_npwk(), pw_.npw);

    for (int igl = 0; igl < kq.get_npwk(); ++igl)
    {
        EXPECT_EQ(kq.get_ig(igl), igl); // ordering preserved
        EXPECT_EQ(kq.get_ig2isz(igl), pw_.ig2isz[igl]);
        const ModuleBase::Vector3<double> gcar = kq.get_gcar(igl);
        EXPECT_DOUBLE_EQ(kq.get_gk2(igl), gcar * gcar);
        const ModuleBase::Vector3<double> gp = kq.get_gpluskq(igl);
        EXPECT_DOUBLE_EQ(gp.x, gcar.x);
        EXPECT_DOUBLE_EQ(gp.y, gcar.y);
        EXPECT_DOUBLE_EQ(gp.z, gcar.z);
    }

    // every selected vector lies inside the cutoff
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
}

TEST_F(DFPTKQBasisTest, ShiftedCenterSelectsAsymmetricSphere)
{
    // k = (0,0,0.5): only G=(0,0,0) and G=(0,0,-1) survive the |G+k|^2 cut
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
               ModuleBase::Vector3<double>(0.0, 0.0, 0.5 * ModuleBase::TWO_PI / lat0_)});

    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(&pw_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 1);
    const ModuleBase::Vector3<double> center = kq.get_kplusq();
    EXPECT_NEAR(center.x, 0.0, 1e-12);
    EXPECT_NEAR(center.y, 0.0, 1e-12);
    EXPECT_NEAR(center.z, 0.5 * ModuleBase::TWO_PI / lat0_, 1e-12);

    // independent brute force on the full FFT grid
    const std::vector<ModuleBase::Vector3<double>> ref = ReferenceSelection(center);
    ASSERT_EQ(ref.size(), 2u); // G=(0,0,0) and G=(0,0,-1)
    const std::vector<ModuleBase::Vector3<double>> sel = KqSet(kq);
    EXPECT_EQ(sel.size(), ref.size());
    for (size_t i = 0; i < ref.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(sel[i].x, ref[i].x);
        EXPECT_DOUBLE_EQ(sel[i].y, ref[i].y);
        EXPECT_DOUBLE_EQ(sel[i].z, ref[i].z);
    }

    // indices must be unique (each selected vector corresponds to exactly
    // one underlying G vector)
    std::vector<int> igs;
    for (int igl = 0; igl < kq.get_npwk(); ++igl)
    {
        igs.push_back(kq.get_ig(igl));
    }
    std::sort(igs.begin(), igs.end());
    for (size_t i = 1; i < igs.size(); ++i)
    {
        EXPECT_NE(igs[i], igs[i - 1]);
    }
}

TEST_F(DFPTKQBasisTest, TranslationInvarianceOfKQ)
{
    // the k+q basis depends only on the sum k+q: (ik=0, q=k1) must agree
    // with (ik=1, q=0), both centered at k1
    const ModuleBase::Vector3<double> k1(0.0, 0.0, 0.5 * ModuleBase::TWO_PI / lat0_);
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0), k1});

    ModuleDFPT::DFPT_KQ_Basis a, b;
    a.init(&pw_, k1, 0);
    b.init(&pw_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 1);
    ASSERT_EQ(a.get_npwk(), b.get_npwk());
    EXPECT_EQ(a.get_npwk(), b.get_npwk());
    for (int igl = 0; igl < a.get_npwk(); ++igl)
    {
        EXPECT_EQ(a.get_ig(igl), b.get_ig(igl));
        EXPECT_DOUBLE_EQ(a.get_gk2(igl), b.get_gk2(igl));
        EXPECT_DOUBLE_EQ(a.get_gpluskq(igl).x, b.get_gpluskq(igl).x);
        EXPECT_DOUBLE_EQ(a.get_gpluskq(igl).y, b.get_gpluskq(igl).y);
        EXPECT_DOUBLE_EQ(a.get_gpluskq(igl).z, b.get_gpluskq(igl).z);
    }

    // shifting back by -k1 recovers the Gamma basis
    ModuleDFPT::DFPT_KQ_Basis c;
    c.init(&pw_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0) - k1, 1);
    EXPECT_EQ(c.get_npwk(), 7);
    for (int igl = 0; igl < c.get_npwk(); ++igl)
    {
        EXPECT_DOUBLE_EQ(c.get_gk2(igl), c.get_gcar(igl) * c.get_gcar(igl));
    }
}

TEST_F(DFPTKQBasisTest, NonzeroQMatchesBruteForce)
{
    const ModuleBase::Vector3<double> k1(0.0, 0.0, 0.5 * ModuleBase::TWO_PI / lat0_);
    const ModuleBase::Vector3<double> q(0.5 * ModuleBase::TWO_PI / lat0_, 0.0, 0.25 * ModuleBase::TWO_PI / lat0_);
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0), k1});

    const ModuleBase::Vector3<double> centers[2] = {q, k1 + q};
    for (int ik = 0; ik < 2; ++ik)
    {
        ModuleDFPT::DFPT_KQ_Basis kq;
        kq.init(&pw_, q, ik);
        const std::vector<ModuleBase::Vector3<double>> ref = ReferenceSelection(centers[ik]);
        const std::vector<ModuleBase::Vector3<double>> sel = KqSet(kq);
        ASSERT_EQ(sel.size(), ref.size());
        for (size_t i = 0; i < ref.size(); ++i)
        {
            EXPECT_DOUBLE_EQ(sel[i].x, ref[i].x);
            EXPECT_DOUBLE_EQ(sel[i].y, ref[i].y);
            EXPECT_DOUBLE_EQ(sel[i].z, ref[i].z);
        }
        // cutoff consistency: every retained plane wave is below the cut
        for (int igl = 0; igl < kq.get_npwk(); ++igl)
        {
            EXPECT_LE(kq.get_gk2(igl), gk_ecut_ + 1e-12);
        }
    }
}

TEST_F(DFPTKQBasisTest, InvalidOrGammaOnlyBaseIsRejected)
{
    // null provider: valid-but-empty basis
    ModuleDFPT::DFPT_KQ_Basis kq;
    kq.init(nullptr, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0);
    EXPECT_FALSE(kq.is_valid());
    EXPECT_EQ(kq.get_npwk(), 0);
    EXPECT_TRUE(kq.get_igl2ig().empty());

    // gamma_only base is rejected (DFPT needs the full complex G ball)
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0)});
    pw_.gamma_only = true;
    pw_.fftny = ny_ / 2 + 1;
    ModuleDFPT::DFPT_KQ_Basis kq2;
    EXPECT_EXIT(kq2.init(&pw_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0),
                ::testing::ExitedWithCode(1),
                "");
}

} // namespace