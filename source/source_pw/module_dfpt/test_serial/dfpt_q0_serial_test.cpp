#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>

// serial unit test of the q -> 0 response (C6): the position operator in
// the velocity (commutator) form, the dielectric tensor and the Born
// charges. Runs without __MPI on the shared FFT grid like the other DFPT
// serial tests; all references are closed-form or operator finite
// differences, no ground-state solver is involved.

#define private public
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/pseudo.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_cell/magnetism.h"
#include "source_pw/module_pwdft/stru_fac.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#include "source_pw/module_dfpt/dfpt_q0.h"
#undef private

#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_psi/psi.h"
#include "dfpt_serial_fixture.h"

// ctor/dtor stubs for the cell/spepot/stru_fac link closures live in the
// shared test/dfpt_test_mocks.cpp compiled into every DFPT test binary.

/************************************************
 *  serial unit test of DFPT_Q0 (C6)
 ***********************************************/

/**
 * - Tested Functions:
 *   - DFPT_Pert::build_vkb_dk: the three analytic derivative terms (atomic
 *     phase, radial chain rule, harmonic direction chain) against a central
 *     finite difference of build_vkb on a generic shifted gk list.
 *   - DFPT_Q0::pos_matrix: the kinetic velocity term, the -i commutator
 *     factor, the tpiba scaling and the Hermitian structure against
 *     closed-form plane-wave-combination states; the nonlocal contraction
 *     against an operator finite difference of <u_m|V_nl(k)|u_n>; exactly
 *     degenerate pairs are skipped.
 *   - DFPT_Q0::compute_eps: the SCF contraction (16 pi / Omega, wg,
 *     occupied-band sum) on synthetic pos_resp/dpsi_efield stashes;
 *     the empty-band rows must be skipped and indices/conj pinned by
 *     distinct per-direction complexes.
 *   - DFPT_Q0::compute_born: elementwise against the closed-form
 *     <u_v|dV/dtau_b|u_m><u_m|r_a|u_v> sums at q = 0 (Coulomb local part),
 *     the ionic Z delta_ab, and the dpsi-slot backup/restore.
 */

class DFPTQ0SerialTest : public DFPTSerialBase
{
  protected:
    Structure_Factor sf_;
    ModuleDFPT::DFPT_Pert pert_;
    ModuleDFPT::DFPT_Q0 q0_;

    // this fixture is Gamma-only: shadow the generic default k of the base
    const ModuleBase::Vector3<double> k_d_{0.0, 0.0, 0.0};
    const ModuleBase::Vector3<double> gx_{0.1, 0.0, 0.0}; // 1/lat0 units
    const ModuleBase::Vector3<double> gy_{0.0, 0.1, 0.0};

    void SetUp() override
    {
        SetUpCell();
        SetupBases(k_d_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 4);
        pert_.init(ucell_, &pw_rho_, &pw_wfc_, sf_);
        q0_.init(ucell_, &pw_rho_, &pw_wfc_, &pert_);
    }

    // wfc-basis index of the reciprocal vector (ix, iy, iz)/a at Gamma
    int IgOf(int ix, int iy, int iz) const
    {
        const int npwk = pw_wfc_.npwk[0];
        for (int ig = 0; ig < npwk; ++ig)
        {
            const ModuleBase::Vector3<double> g = pw_wfc_.getgpluskcar(0, ig);
            if (std::llround(g.x * a_) == ix && std::llround(g.y * a_) == iy
                && std::llround(g.z * a_) == iz)
            {
                return ig;
            }
        }
        return -1;
    }
};

// ---------------------------------------------------------------------------
// build_vkb_dk against a finite difference of build_vkb
// ---------------------------------------------------------------------------

TEST_F(DFPTQ0SerialTest, BuildVkbDkMatchesFiniteDifference)
{
    MakeNCAtom();
    // generic shifted list without any |g| = 0 entry (the direction chain is
    // singular at the origin for l >= 1 rows)
    const int ng = 8;
    std::vector<ModuleBase::Vector3<double>> gk(ng);
    gk[0] = ModuleBase::Vector3<double>(0.1, 0.0, 0.0) + k_d_;
    gk[1] = ModuleBase::Vector3<double>(0.0, 0.1, 0.0) + k_d_;
    gk[2] = ModuleBase::Vector3<double>(0.07, 0.13, 0.05) + k_d_;
    gk[3] = ModuleBase::Vector3<double>(-0.11, 0.23, -0.17) + k_d_;
    gk[4] = ModuleBase::Vector3<double>(0.29, -0.19, 0.31) + k_d_;
    gk[5] = ModuleBase::Vector3<double>(-0.05, 0.0, 0.11) + k_d_;
    gk[6] = ModuleBase::Vector3<double>(0.13, -0.07, 0.03) + k_d_;
    gk[7] = ModuleBase::Vector3<double>(0.21, 0.11, -0.13) + k_d_;

    std::vector<std::vector<std::complex<double>>> vkb;
    pert_.build_vkb(0, 0, gk, vkb);
    const int nh = ucell_.atoms[0].ncpp.nh;
    ASSERT_EQ(vkb.size(), static_cast<size_t>(nh));

    const double eps = 1.0e-5; // gcar units
    for (int d = 0; d < 3; ++d)
    {
        ModuleBase::Vector3<double> shift(0.0, 0.0, 0.0);
        shift[d] = eps;
        std::vector<ModuleBase::Vector3<double>> gk_p(ng), gk_m(ng);
        for (int i = 0; i < ng; ++i)
        {
            gk_p[i] = gk[i] + shift;
            gk_m[i] = gk[i] - shift;
        }
        std::vector<std::vector<std::complex<double>>> vkb_p, vkb_m;
        pert_.build_vkb(0, 0, gk_p, vkb_p);
        pert_.build_vkb(0, 0, gk_m, vkb_m);

        std::vector<std::vector<std::complex<double>>> dvkb;
        pert_.build_vkb_dk(0, 0, d, gk, vkb, dvkb);
        ASSERT_EQ(dvkb.size(), static_cast<size_t>(nh));
        for (int mu = 0; mu < nh; ++mu)
        {
            for (int i = 0; i < ng; ++i)
            {
                const std::complex<double> fd = (vkb_p[mu][i] - vkb_m[mu][i]) / (2.0 * eps);
                const double scale = std::max(1.0, std::abs(fd));
                EXPECT_NEAR(dvkb[mu][i].real(), fd.real(), 1.0e-5 * scale)
                    << "mu=" << mu << " i=" << i << " d=" << d;
                EXPECT_NEAR(dvkb[mu][i].imag(), fd.imag(), 1.0e-5 * scale)
                    << "mu=" << mu << " i=" << i << " d=" << d;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// pos_matrix: kinetic part, -i factor, Hermiticity, degenerate skip
// ---------------------------------------------------------------------------

TEST_F(DFPTQ0SerialTest, PosMatrixKineticAndDegenerateSkip)
{
    const int npwk = pw_wfc_.npwk[0];
    const int ig0 = IgOf(0, 0, 0);
    const int igx = IgOf(1, 0, 0);
    const int igy = IgOf(0, 1, 0);
    ASSERT_GE(ig0, 0);
    ASSERT_GE(igx, 0);
    ASSERT_GE(igy, 0);

    const double e1 = ucell_.tpiba2 * (gx_ * gx_); // |Gx|^2 in Ry
    // b0 = sqrt(0.8)|G0> + sqrt(0.2)|Gx>   (eps = 0.2 e1)
    // b1 = sqrt(0.2)|G0> - sqrt(0.8)|Gx>   (eps = 0.8 e1)
    // b2 = (|Gx> + |Gy>)/sqrt(2)           (eps = e1, degenerate with b3)
    // b3 = (|Gx> - |Gy>)/sqrt(2)           (eps = e1)
    psi::Psi<std::complex<double>> psi(1, 4, npwk, npwk, true);
    psi.zero_out();
    psi(0, 0, ig0) = std::sqrt(0.8);
    psi(0, 0, igx) = std::sqrt(0.2);
    psi(0, 1, ig0) = std::sqrt(0.2);
    psi(0, 1, igx) = -std::sqrt(0.8);
    psi(0, 2, igx) = std::complex<double>(1.0 / std::sqrt(2.0), 0.0);
    psi(0, 2, igy) = std::complex<double>(1.0 / std::sqrt(2.0), 0.0);
    psi(0, 3, igx) = std::complex<double>(1.0 / std::sqrt(2.0), 0.0);
    psi(0, 3, igy) = std::complex<double>(-1.0 / std::sqrt(2.0), 0.0);

    ModuleBase::matrix eig(1, 4);
    eig(0, 0) = 0.2 * e1;
    eig(0, 1) = 0.8 * e1;
    eig(0, 2) = e1;
    eig(0, 3) = e1;

    std::vector<std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>> r_mat;
    q0_.pos_matrix(psi, eig, r_mat);
    ASSERT_EQ(r_mat.size(), 1u);
    ASSERT_EQ(r_mat[0].size(), 4u);

    // p_01^d = 2 tpiba^2 <b0| (k+G)_d |b1> = 2 tpiba^2 (-sqrt(0.16)) Gx_d
    for (int d = 0; d < 3; ++d)
    {
        const std::complex<double> p01(2.0 * ucell_.tpiba2 * (-std::sqrt(0.16)) * gx_[d], 0.0);
        const std::complex<double> expect
            = std::complex<double>(0.0, -1.0) * p01 / (ucell_.tpiba * (0.2 * e1 - 0.8 * e1));
        EXPECT_NEAR(r_mat[0][0][1][d].real(), expect.real(), 1.0e-10);
        EXPECT_NEAR(r_mat[0][0][1][d].imag(), expect.imag(), 1.0e-10);
        // Hermiticity: r_10 = conj(r_01)
        EXPECT_NEAR(r_mat[0][1][0][d].real(), expect.real(), 1.0e-10);
        EXPECT_NEAR(r_mat[0][1][0][d].imag(), -expect.imag(), 1.0e-10);
        // diagonal vanishes
        EXPECT_EQ(r_mat[0][0][0][d], std::complex<double>(0.0, 0.0));
        // exactly degenerate pairs are skipped
        EXPECT_EQ(r_mat[0][2][3][d], std::complex<double>(0.0, 0.0));
        EXPECT_EQ(r_mat[0][3][2][d], std::complex<double>(0.0, 0.0));
    }
}

// ---------------------------------------------------------------------------
// pos_matrix: nonlocal velocity against an operator finite difference
// ---------------------------------------------------------------------------

TEST_F(DFPTQ0SerialTest, PosMatrixNonlocalMatchesOperatorFiniteDifference)
{
    MakeNCAtom();
    const int npwk = pw_wfc_.npwk[0];
    const int ig0 = IgOf(0, 0, 0);
    ASSERT_GE(ig0, 0);

    // deterministic pseudo-random orthonormal bands with the |G| = 0
    // component forced to zero: the projector derivative is direction
    // singular exactly at g = 0 (l >= 1 rows), so that single column is
    // excluded from both sides of the comparison
    const int nb = 4;
    psi::Psi<std::complex<double>> psi(1, nb, npwk, npwk, true);
    psi.zero_out();
    unsigned seed = 20260817u;
    auto rnd = [&]()
    {
        seed = seed * 1664525u + 1013904223u;
        return ((seed >> 8) & 0xffffff) / 16777216.0 * 2.0 - 1.0;
    };
    std::vector<std::vector<std::complex<double>>> c(nb, std::vector<std::complex<double>>(npwk));
    for (int b = 0; b < nb; ++b)
    {
        for (int ig = 0; ig < npwk; ++ig)
        {
            c[b][ig] = (ig == ig0) ? std::complex<double>(0.0, 0.0)
                                   : std::complex<double>(rnd(), rnd());
        }
    }
    // Gram-Schmidt, skipping the zero column keeps the norm from column 1 on
    for (int b = 0; b < nb; ++b)
    {
        for (int p = 0; p < b; ++p)
        {
            std::complex<double> dot(0.0, 0.0);
            for (int ig = 0; ig < npwk; ++ig)
            {
                dot += std::conj(c[p][ig]) * c[b][ig];
            }
            for (int ig = 0; ig < npwk; ++ig)
            {
                c[b][ig] -= dot * c[p][ig];
            }
        }
        double nrm = 0.0;
        for (int ig = 0; ig < npwk; ++ig)
        {
            nrm += std::norm(c[b][ig]);
        }
        nrm = std::sqrt(nrm);
        for (int ig = 0; ig < npwk; ++ig)
        {
            c[b][ig] /= nrm;
            psi(0, b, ig) = c[b][ig];
        }
    }

    // non-degenerate fake eigenvalues (pos_matrix only uses them as divisors)
    ModuleBase::matrix eig(1, nb);
    for (int b = 0; b < nb; ++b)
    {
        eig(0, b) = 0.31 + 0.17 * b;
    }

    std::vector<std::vector<std::vector<ModuleBase::Vector3<std::complex<double>>>>> r_mat;
    q0_.pos_matrix(psi, eig, r_mat);

    // finite-difference reference of the full velocity operator
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

    std::vector<ModuleBase::Vector3<double>> gk(npwk);
    for (int ig = 0; ig < npwk; ++ig)
    {
        gk[ig] = pw_wfc_.getgpluskcar(0, ig);
    }
    const double eps = 1.0e-5;

    // becp with the |G| = 0 column dropped on a shifted list
    auto vnl_matrix = [&](const std::vector<ModuleBase::Vector3<double>>& glist,
                          std::vector<std::vector<std::complex<double>>>& mmat)
    {
        std::vector<std::vector<std::complex<double>>> vkb;
        pert_.build_vkb(0, 0, glist, vkb);
        std::vector<std::vector<std::complex<double>>> becp(nb);
        for (int b = 0; b < nb; ++b)
        {
            becp[b].assign(nh, std::complex<double>(0.0, 0.0));
            for (int mu = 0; mu < nh; ++mu)
            {
                for (int ig = 0; ig < npwk; ++ig)
                {
                    if (ig == ig0)
                    {
                        continue; // the singular column, see above
                    }
                    becp[b][mu] += std::conj(vkb[mu][ig]) * psi(0, b, ig);
                }
            }
        }
        mmat.assign(nb, std::vector<std::complex<double>>(nb, std::complex<double>(0.0, 0.0)));
        for (int m = 0; m < nb; ++m)
        {
            for (int n = 0; n < nb; ++n)
            {
                for (int mu = 0; mu < nh; ++mu)
                {
                    std::complex<double> dc(0.0, 0.0);
                    for (int nu = 0; nu < nh; ++nu)
                    {
                        if (row_m[mu] == row_m[nu])
                        {
                            dc += p.dion(row_ib[mu], row_ib[nu]) * becp[n][nu];
                        }
                    }
                    mmat[m][n] += std::conj(becp[m][mu]) * dc;
                }
            }
        }
    };

    for (int d = 0; d < 3; ++d)
    {
        ModuleBase::Vector3<double> shift(0.0, 0.0, 0.0);
        shift[d] = eps;
        std::vector<ModuleBase::Vector3<double>> gk_p(npwk), gk_m(npwk);
        for (int ig = 0; ig < npwk; ++ig)
        {
            gk_p[ig] = gk[ig] + shift;
            gk_m[ig] = gk[ig] - shift;
        }
        std::vector<std::vector<std::complex<double>>> mm_p, mm_m;
        vnl_matrix(gk_p, mm_p);
        vnl_matrix(gk_m, mm_m);

        for (int m = 0; m < nb; ++m)
        {
            for (int n = 0; n < nb; ++n)
            {
                if (m == n)
                {
                    continue;
                }
                const double de = eig(0, m) - eig(0, n);
                // recover p from r: r = -i p / (tpiba de)
                const std::complex<double> p_r
                    = std::complex<double>(0.0, 1.0) * ucell_.tpiba * de * r_mat[0][m][n][d];
                // analytic kinetic + finite-difference nonlocal
                std::complex<double> p_kin(0.0, 0.0);
                for (int ig = 0; ig < npwk; ++ig)
                {
                    p_kin += 2.0 * ucell_.tpiba2 * gk[ig][d] * std::conj(psi(0, m, ig))
                             * psi(0, n, ig);
                }
                const std::complex<double> p_nl = (mm_p[m][n] - mm_m[m][n]) / (2.0 * eps);
                const double scale = std::max(1.0, std::abs(p_kin) + std::abs(p_nl));
                EXPECT_NEAR(p_r.real(), (p_kin + p_nl).real(), 1.0e-6 * scale)
                    << "m=" << m << " n=" << n << " d=" << d;
                EXPECT_NEAR(p_r.imag(), (p_kin + p_nl).imag(), 1.0e-6 * scale)
                    << "m=" << m << " n=" << n << " d=" << d;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// compute_eps (SCF contraction) against synthetic response stashes:
// eps(a,b) = delta_ab - (16 pi/Omega) sum_k wg sum_occ Re<Y^a_v|dpsi^E,b_v>
// with distinct complexes per direction catching transposed indices and
// conj placement; the empty-band rows stay unsolved and must be skipped
// ---------------------------------------------------------------------------

TEST_F(DFPTQ0SerialTest, ComputeEpsScfSyntheticStash)
{
    const int npwk = pw_wfc_.npwk[0];
    const int ig0 = IgOf(0, 0, 0);
    const int igx = IgOf(1, 0, 0);
    ASSERT_GE(ig0, 0);
    ASSERT_GE(igx, 0);

    ModuleBase::matrix wg(1, 2);
    wg(0, 0) = 2.0;
    wg(0, 1) = 0.0;

    // synthetic bare position legs Y^a_{0,0} = P_c x_a|psi_0>
    const std::complex<double> gam[3] = {std::complex<double>(0.15, -0.3),
                                         std::complex<double>(0.4, 0.05),
                                         std::complex<double>(-0.35, 0.2)};
    const std::complex<double> del[3] = {std::complex<double>(-0.25, 0.45),
                                         std::complex<double>(0.1, -0.1),
                                         std::complex<double>(0.3, 0.25)};
    for (int a = 0; a < 3; ++a)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> y(
            1, std::vector<std::vector<std::complex<double>>>(2));
        y[0][0].assign(npwk, std::complex<double>(0.0, 0.0));
        y[0][0][ig0] = gam[a];
        y[0][0][igx] = del[a];
        data_.set_pos_resp(a, y);
    }

    // synthetic converged E-field responses dpsi^E,b_{0,0}
    const std::complex<double> mue[3] = {std::complex<double>(0.3, 0.2),
                                         std::complex<double>(-0.1, 0.4),
                                         std::complex<double>(0.25, -0.15)};
    const std::complex<double> nue[3] = {std::complex<double>(0.2, -0.35),
                                         std::complex<double>(0.45, 0.1),
                                         std::complex<double>(-0.2, -0.05)};
    for (int b = 0; b < 3; ++b)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> e(
            1, std::vector<std::vector<std::complex<double>>>(2));
        e[0][0].assign(npwk, std::complex<double>(0.0, 0.0));
        e[0][0][ig0] = mue[b];
        e[0][0][igx] = nue[b];
        data_.set_dpsi_efield(b, e);
    }

    q0_.compute_eps(wg, data_);
    const ModuleBase::matrix eps = data_.get_dielectric();

    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            // <Y^a|dE^b> = conj(gam_a) mue_b + conj(del_a) nue_b over the
            // shared G support, wg-weighted with the 16 pi/Omega prefactor
            const std::complex<double> dot = std::conj(gam[a]) * mue[b]
                                             + std::conj(del[a]) * nue[b];
            const double expect = ((a == b) ? 1.0 : 0.0)
                                  - 16.0 * ModuleBase::PI / ucell_.omega
                                        * wg(0, 0) * dot.real();
            EXPECT_NEAR(eps(a, b), expect, 1.0e-12) << "a=" << a << " b=" << b;
        }
    }
}

// ---------------------------------------------------------------------------
// compute_born against the closed-form screened-leg product (v4 QE anchor):
// Z*(a,idir) = zion delta - 2 sum wg Re<dpsi^{idir}(scf)|Y^a> with synthetic
// stashed displacement responses and position legs
// ---------------------------------------------------------------------------

TEST_F(DFPTQ0SerialTest, ComputeBornTwoLevelAnalytic)
{
    const int npwk = pw_wfc_.npwk[0];
    const int ig0 = IgOf(0, 0, 0);
    const int igx = IgOf(1, 0, 0);
    ASSERT_GE(ig0, 0);
    ASSERT_GE(igx, 0);

    const double e1 = ucell_.tpiba2 * (gx_ * gx_);
    psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
    psi.zero_out();
    psi(0, 0, ig0) = std::sqrt(0.6);
    psi(0, 0, igx) = std::sqrt(0.4);
    psi(0, 1, ig0) = std::sqrt(0.2);

    ModuleBase::matrix wg(1, 2);
    wg(0, 0) = 2.0;
    wg(0, 1) = 0.0;
    ModuleBase::matrix eig(1, 2);
    eig(0, 0) = 0.4 * e1;
    eig(0, 1) = e1;

    // sentinel dpsi in the q = 0 slot: compute_born never touches the dpsi
    // slots, the sentinel must survive verbatim
    const std::vector<std::complex<double>> sentinel(npwk, std::complex<double>(0.5, -0.25));
    data_.set_dpsi(0, 0, 0, sentinel);

    // synthetic converged displacement responses dpsi(scf)/du_{0,idir} for
    // the occupied band (G0/Gx components, distinct complexes per idir
    // catch transposed indices); the empty-band row stays unsolved
    const std::complex<double> alpha[3] = {std::complex<double>(0.3, 0.2),
                                           std::complex<double>(-0.1, 0.4),
                                           std::complex<double>(0.25, -0.15)};
    const std::complex<double> beta[3] = {std::complex<double>(0.2, -0.35),
                                          std::complex<double>(0.45, 0.1),
                                          std::complex<double>(-0.2, -0.05)};
    for (int idir = 0; idir < 3; ++idir)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> disp(
            1, std::vector<std::vector<std::complex<double>>>(2));
        disp[0][0].assign(npwk, std::complex<double>(0.0, 0.0));
        disp[0][0][ig0] = alpha[idir];
        disp[0][0][igx] = beta[idir];
        data_.set_dpsi_disp(0, idir, disp);
    }

    // synthetic solved position legs Y^a_{0,0} = P_c x_a|psi_0>
    const std::complex<double> gam[3] = {std::complex<double>(0.15, -0.3),
                                         std::complex<double>(0.4, 0.05),
                                         std::complex<double>(-0.35, 0.2)};
    const std::complex<double> del[3] = {std::complex<double>(-0.25, 0.45),
                                         std::complex<double>(0.1, -0.1),
                                         std::complex<double>(0.3, 0.25)};
    for (int a = 0; a < 3; ++a)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> y(
            1, std::vector<std::vector<std::complex<double>>>(2));
        y[0][0].assign(npwk, std::complex<double>(0.0, 0.0));
        y[0][0][ig0] = gam[a];
        y[0][0][igx] = del[a];
        data_.set_pos_resp(a, y);
    }

    q0_.compute_born(psi, wg, eig, data_);
    const ModuleBase::matrix zstar = data_.get_born(0);

    const double zion = ucell_.atoms[0].ncpp.zv;
    for (int idir = 0; idir < 3; ++idir)
    {
        for (int a = 0; a < 3; ++a)
        {
            // <dpsi^{idir}|Y^a> = conj(alpha)gam + conj(beta)del over the
            // shared G support, wg-weighted with the -2 spin prefactor
            const std::complex<double> dot = std::conj(alpha[idir]) * gam[a]
                                             + std::conj(beta[idir]) * del[a];
            const double expect = ((a == idir) ? zion : 0.0)
                                  - 2.0 * wg(0, 0) * dot.real();
            EXPECT_NEAR(zstar(a, idir), expect, 1.0e-12) << "a=" << a << " idir=" << idir;
        }
    }

    // the q = 0 dpsi slot is untouched
    const std::vector<std::complex<double>> after = data_.get_dpsi(0, 0, 0);
    ASSERT_EQ(after.size(), sentinel.size());
    for (size_t i = 0; i < after.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(after[i].real(), sentinel[i].real());
        EXPECT_DOUBLE_EQ(after[i].imag(), sentinel[i].imag());
    }
}

// ---------------------------------------------------------------------------
// build_stars + rotate_tensor: C3 group on a simple-cubic 3-atom orbit cell.
// Atoms sit on the orbit tau = (0.1,0.2,0.3) under the cyclic permutation
// (x,y,z)->(z,x,y), so every group operation maps atom i -> i+1 (mod 3) and
// the star of k = (1/4,0,0) has the three members
// {(1/4,0,0),(0,1/4,0),(0,0,1/4)}. An anisotropic trace-6 tensor star-
// averages to 2*delta_ab, and the members carry the cyclic atom map.
// ---------------------------------------------------------------------------

TEST_F(DFPTQ0SerialTest, StarRotationCyclicGroup)
{
    // rebuild the cell as 3 atoms on the C3 orbit of (0.1,0.2,0.3)
    ucell_.nat = 3;
    ucell_.atoms[0].na = 3;
    ucell_.atoms[0].tau.resize(3);
    ucell_.atoms[0].taud.resize(3);
    const double t0[3] = {0.1, 0.2, 0.3};
    for (int i = 0; i < 3; ++i)
    {
        ucell_.atoms[0].taud[i] = ModuleBase::Vector3<double>(t0[i], t0[(i + 1) % 3], t0[(i + 2) % 3]);
        ucell_.atoms[0].tau[i] = ucell_.atoms[0].taud[i] * a_;
    }
    delete[] ucell_.iat2it;
    delete[] ucell_.iat2ia;
    ucell_.iat2it = new int[3];
    ucell_.iat2ia = new int[3];
    for (int i = 0; i < 3; ++i)
    {
        ucell_.iat2it[i] = 0;
        ucell_.iat2ia[i] = i;
    }

    // C3 about (111): direct-space row-convention matrices g with tau' = tau*g
    // (the cycle (x,y,z)->(y,z,x) sends atom i -> i+1 mod 3); g2 = g*g. In
    // reciprocal space kgmatrix = G*g*G^-1 = g for a cubic cell.
    const ModuleBase::Matrix3 g1(0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    const ModuleBase::Matrix3 g2 = g1 * g1;
    ucell_.symm.nrotk = 3;
    ucell_.symm.gmatrix[0] = ModuleBase::Matrix3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    ucell_.symm.gmatrix[1] = g1;
    ucell_.symm.gmatrix[2] = g2;
    for (int j = 0; j < 3; ++j)
    {
        ucell_.symm.kgmatrix[j] = ucell_.symm.gmatrix[j];
        ucell_.symm.gtrans[j] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    }

    // single reduced k point (1/4,0,0) on its own wfc basis
    ModulePW::PW_Basis_K kwfc;
    const ModuleBase::Vector3<double> klist[1]
        = {ModuleBase::Vector3<double>(0.25, 0.0, 0.0)};
    kwfc.initgrids(lat0_, latvec_, pw_rho_.nx, pw_rho_.ny, pw_rho_.nz);
    kwfc.initparameters(false, ecutwfc_, 1, klist);
    kwfc.fft_bundle.initfftmode(0);
    kwfc.setuptransform();
    kwfc.collect_local_pw();
    q0_.init(ucell_, &pw_rho_, &kwfc, &pert_);

    q0_.build_stars(1);
    ASSERT_EQ(q0_.stars_[0].size(), (size_t)3);

    // star-average an anisotropic trace-6 tensor over the members
    ModuleBase::matrix chi(3, 3, true);
    chi(0, 0) = 1.0;
    chi(1, 1) = 2.0;
    chi(2, 2) = 3.0;
    double avg[9] = {0.0};
    for (size_t im = 0; im < q0_.stars_[0].size(); ++im)
    {
        double rot[9];
        q0_.rotate_tensor(q0_.stars_[0][im].cart, chi, rot);
        for (int i = 0; i < 9; ++i)
        {
            avg[i] += rot[i] / 3.0;
        }
    }
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            EXPECT_NEAR(avg[3 * a + b], (a == b) ? 2.0 : 0.0, 1.0e-12) << "a=" << a << " b=" << b;
        }
    }

    // every member carries a cyclic atom map: the identity member is the
    // pre-filled one with an empty map, the two rotation members carry the
    // g and g*g maps i -> (i+1)%3 and i -> (i+2)%3
    int seen_shift[3] = {1, 0, 0}; // identity shift 0 already seen
    for (size_t im = 0; im < q0_.stars_[0].size(); ++im)
    {
        const std::vector<int>& amap = q0_.stars_[0][im].atom_map;
        if (amap.empty())
        {
            // identity member: rotation must be the identity
            EXPECT_DOUBLE_EQ(q0_.stars_[0][im].cart.e11, 1.0);
            EXPECT_DOUBLE_EQ(q0_.stars_[0][im].cart.e12, 0.0);
            continue;
        }
        ASSERT_EQ(amap.size(), (size_t)3);
        int shift = -1;
        for (int s = 0; s < 3; ++s)
        {
            bool ok = true;
            for (int i = 0; i < 3; ++i)
            {
                ok = ok && amap[i] == (i + s) % 3;
            }
            if (ok)
            {
                shift = s;
                break;
            }
        }
        ASSERT_GE(shift, 0);
        seen_shift[shift] = 1;
    }
    for (int s = 0; s < 3; ++s)
    {
        EXPECT_EQ(seen_shift[s], 1) << "atom-map shift " << s << " missing";
    }

    // with the point group unavailable the stars degenerate to identity
    ucell_.symm.nrotk = 0;
    q0_.build_stars(1);
    ASSERT_EQ(q0_.stars_[0].size(), (size_t)1);
    double rot[9];
    q0_.rotate_tensor(q0_.stars_[0][0].cart, chi, rot);
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            EXPECT_NEAR(rot[3 * a + b], chi(a, b), 1.0e-12);
        }
    }
}
