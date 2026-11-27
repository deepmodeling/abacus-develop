#include "source_pw/module_dfpt/dfpt_kq_basis.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <iostream>
#include <vector>

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

namespace
{

bool VecLess(const ModuleBase::Vector3<double> &a, const ModuleBase::Vector3<double> &b)
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
    void FillSticks(ModulePW::PW_Basis &pw, double ggecut)
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

    // per-k plane-wave data (npwk/gcar/gk2/igl2isz_k) in the production
    // PW_Basis_K layout: the DFPT Gamma special case reads exactly these
    // arrays, so the fixture must fill them for every k of the mesh. The
    // arrays are allocated once for the whole mesh; each call fills only
    // its own k slot, so BuildBase must call this once per k AFTER nst/
    // npw of the layout exist and npwk_max is known (two-pass below)
    void FillSticksK(
        ModulePW::PW_Basis_K &pw,
        int ik_fill,
        const ModuleBase::Vector3<double> &kvec,
        double ggecut,
        int npwk_max)
    {
        std::vector<int> isz_list;
        std::vector<ModuleBase::Vector3<double>> gcar_list;
        std::vector<double> gk2_list;
        for (int is = 0; is < pw.nst; ++is)
        {
            const int ixy = pw.is2fftixy[is];
            const int ix = ixy / pw.fftny;
            const int iy = ixy % pw.fftny;
            const int wix = WrapIndex(ix, pw.nx);
            const int wiy = WrapIndex(iy, pw.ny);
            for (int iz0 = 0; iz0 < pw.nz; ++iz0)
            {
                // the stick members in ig2isz order (z ascending); recover
                // the member index by matching isz = is*nz + iz0
                const int isz = is * pw.nz + iz0;
                const int wiz = WrapIndex(iz0, pw.nz);
                ModuleBase::Vector3<double> f(wix, wiy, wiz);
                const ModuleBase::Vector3<double> g = f * G_;
                const ModuleBase::Vector3<double> gk = g + kvec;
                if (gk * gk <= ggecut)
                {
                    isz_list.push_back(isz);
                    gcar_list.push_back(g);
                    gk2_list.push_back(gk * gk);
                }
            }
        }
        pw.npwk[ik_fill] = static_cast<int>(isz_list.size());
        for (int i = 0; i < pw.npwk[ik_fill]; ++i)
        {
            pw.gcar[ik_fill * npwk_max + i] = gcar_list[i];
            pw.gk2[ik_fill * npwk_max + i] = gk2_list[i];
            pw.igl2isz_k[ik_fill * npwk_max + i] = isz_list[i];
        }
    }

    // shared field setup for a cubic complex basis. The wavefunction grid
    // uses the production radius (sqrt(gk_ecut) + max|k_mesh|)^2; the
    // charge grid uses the production 4*ecutwfc (= 4*gk_ecut ball, no
    // k list) that covers every q-shifted ball for q in the first BZ.
    void BuildBase(const std::vector<ModuleBase::Vector3<double>> &kvec_c)
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
        // two-pass per-k fill: first pass determines npwk_max, second pass
        // allocates the arrays once and fills each k slot
        {
            int npwk_max = 0;
            for (size_t i = 0; i < kvec_c.size(); ++i)
            {
                int cnt = 0;
                for (int is = 0; is < pw_.nst; ++is)
                {
                    const int ixy = pw_.is2fftixy[is];
                    const int ix = ixy / pw_.fftny;
                    const int iy = ixy % pw_.fftny;
                    const int wix = WrapIndex(ix, pw_.nx);
                    const int wiy = WrapIndex(iy, pw_.ny);
                    for (int iz0 = 0; iz0 < pw_.nz; ++iz0)
                    {
                        const int wiz = WrapIndex(iz0, pw_.nz);
                        ModuleBase::Vector3<double> f(wix, wiy, wiz);
                        const ModuleBase::Vector3<double> g = f * G_;
                        const ModuleBase::Vector3<double> gk = g + kvec_c[i];
                        if (gk * gk <= gk_ecut_)
                        {
                            ++cnt;
                        }
                    }
                }
                npwk_max = std::max(npwk_max, cnt);
            }
            const int nks = pw_.nks;
            pw_.npwk = new int[nks];
            pw_.npwk_max = npwk_max;
            pw_.gcar = new ModuleBase::Vector3<double>[npwk_max * nks];
            pw_.gk2 = new double[npwk_max * nks];
            pw_.igl2isz_k = new int[npwk_max * nks];
            for (int i = 0; i < nks; ++i)
            {
                FillSticksK(pw_, i, kvec_c[i], gk_ecut_, npwk_max);
            }
        }

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
    std::vector<ModuleBase::Vector3<double>> ReferenceSelection(const ModuleBase::Vector3<double> &center)
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

    std::vector<ModuleBase::Vector3<double>> KqSet(const ModuleDFPT::DFPT_KQ_Basis &kq)
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
    const std::vector<ModuleBase::Vector3<double>> ref = ReferenceSelection(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
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
    BuildBase({ModuleBase::Vector3<double>(0.0, 0.0, 0.0), ModuleBase::Vector3<double>(0.0, 0.0, 0.5 * b)});

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
    EXPECT_EXIT(
        kq2.init(&pw_, &prho_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0), ::testing::ExitedWithCode(1), "");

    // mismatched FFT grid dimensions between the two bases are rejected
    pw_.gamma_only = false;
    pw_.fftny = ny_;
    prho_.nx = 9;
    prho_.nxyz = 9 * prho_.ny * prho_.nz;
    ModuleDFPT::DFPT_KQ_Basis kq3;
    EXPECT_EXIT(
        kq3.init(&pw_, &prho_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 0), ::testing::ExitedWithCode(1), "");
}

TEST_F(DFPTKQBasisTest, GammaQOnMultiKMeshFollowsGroundStateOrder)
{
    // order lock for the q = 0 legs: on a multi-k production mesh, the
    // Gamma ball of EVERY k must be enumerated in that k's ground-state
    // pw order (stick layout), because the position/field response builds
    // its right-hand side in that order and passes it to the operator and
    // to the occupied-state projectors without remapping. A basis that
    // returns the same set in a different order silently corrupts the
    // dielectric and Born tensors through dfpt_compute_q0.
    ResetBasis();
    BuildBase(
        {ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
            ModuleBase::Vector3<double>(0.0, 0.5 * ModuleBase::TWO_PI / lat0_, 0.5 * ModuleBase::TWO_PI / lat0_)});

    for (int ik = 0; ik < 2; ++ik)
    {
        ModuleDFPT::DFPT_KQ_Basis kq;
        kq.init(&pw_, &prho_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), ik);
        ASSERT_TRUE(kq.is_valid());
        // each k's Gamma ball is that k's ground-state ball, enumerated in
        // the stick layout order (ig2isz traversal with the |g+k|^2
        // membership filter), so the order matches what dfpt_pw_q0.cpp
        // reads
        const int npwk_expected = pw_.npwk[ik];
        EXPECT_EQ(kq.get_npwk(), npwk_expected);
        const int hw_x = static_cast<int>(pw_.nx / 2) + 1;
        const int hw_y = static_cast<int>(pw_.ny / 2) + 1;
        const int hw_z = static_cast<int>(pw_.nz / 2) + 1;
        int igl_out = 0;
        for (int igl = 0; igl < pw_.npw; ++igl)
        {
            const int isz = pw_.ig2isz[igl];
            const int iz = isz % pw_.nz;
            const int is = isz / pw_.nz;
            const int ixy = pw_.is2fftixy[is];
            const int ix = ixy / pw_.fftny;
            const int iy = ixy % pw_.fftny;
            const int wx = (ix >= hw_x) ? ix - pw_.nx : ix;
            const int wy = (iy >= hw_y) ? iy - pw_.ny : iy;
            const int wz = (iz >= hw_z) ? iz - pw_.nz : iz;
            const ModuleBase::Vector3<double> gfrac(wx, wy, wz);
            const ModuleBase::Vector3<double> gcart = gfrac * G_;
            const ModuleBase::Vector3<double> gplusk = gcart + pw_.kvec_c[ik];
            if (gplusk * gplusk > pw_.gk_ecut)
            {
                continue;
            }
            ASSERT_LT(igl_out, kq.get_npwk());
            // sel is the k+q ball in cartesian coordinates with the
            // gplusk convention (compare against g + k, not g alone)
            const ModuleBase::Vector3<double> sel = kq.get_gpluskq(igl_out);
            EXPECT_DOUBLE_EQ(sel.x, gplusk.x);
            EXPECT_DOUBLE_EQ(sel.y, gplusk.y);
            EXPECT_DOUBLE_EQ(sel.z, gplusk.z);
            ++igl_out;
        }
        EXPECT_EQ(igl_out, npwk_expected);
    }
}

} // namespace
