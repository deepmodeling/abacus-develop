// Serial unit tests of the DFPT_Q0 core (dfpt_q0.cpp): build_vkb_dk
// against a central finite difference of build_vkb, compute_eps (the SCF
// contraction on synthetic response stashes), compute_born (the
// closed-form screened-leg product with the dpsi-slot sentinel), and the
// build_stars + rotate_tensor cyclic-group regression.

#include "dfpt_test_tool.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_dfpt/dfpt_pw_data.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_cell/stru_fac.h"
#include "source_psi/psi.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#include "source_pw/module_dfpt/dfpt_q0.h"

// per-file self-contained fixture (gtest style): owns its objects and
// calls the dfpt_test free-function utilities; Gamma-only
class DFPTQ0SerialTest : public ::testing::Test
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
    ModuleDFPT::DFPT_Q0 q0_;

    const ModuleBase::Vector3<double> k_d_{0.0, 0.0, 0.0};
    const ModuleBase::Vector3<double> gx_{0.1, 0.0, 0.0}; // 1/lat0 units

    void SetUp() override
    {
        dfpt_test::setup_cell(ucell_, latvec_);
        dfpt_test::setup_bases(pw_rho_, pw_wfc_, qlist_, data_, ucell_,
                               k_d_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 4, q_cart_);
        pert_.init(ucell_, &pw_rho_, &pw_wfc_, sf_);
        q0_.init(ucell_, &pw_rho_, &pw_wfc_, &pert_);
    }

    void TearDown() override
    {
        delete[] ucell_.atoms;
        ucell_.atoms = nullptr;
    }

    // wfc-basis index of the reciprocal vector (ix, iy, iz)/a at Gamma
    int IgOf(int ix, int iy, int iz) const
    {
        const int npwk = pw_wfc_.npwk[0];
        for (int ig = 0; ig < npwk; ++ig)
        {
            const ModuleBase::Vector3<double> g = pw_wfc_.getgpluskcar(0, ig);
            if (std::llround(g.x * dfpt_test::cubic_edge) == ix
                && std::llround(g.y * dfpt_test::cubic_edge) == iy
                && std::llround(g.z * dfpt_test::cubic_edge) == iz)
            {
                return ig;
            }
        }
        return -1;
    }
};

#include <cmath>
#include <complex>
#include <vector>

TEST_F(DFPTQ0SerialTest, BuildVkbDkMatchesFiniteDifference)
{
    dfpt_test::make_nc_atom(ucell_);
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
                EXPECT_NEAR(dvkb[mu][i].real(), fd.real(), 1.0e-5 * scale) << "mu=" << mu << " i=" << i << " d=" << d;
                EXPECT_NEAR(dvkb[mu][i].imag(), fd.imag(), 1.0e-5 * scale) << "mu=" << mu << " i=" << i << " d=" << d;
            }
        }
    }
}

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
    const std::complex<double> gam[3]
        = {std::complex<double>(0.15, -0.3), std::complex<double>(0.4, 0.05), std::complex<double>(-0.35, 0.2)};
    const std::complex<double> del[3]
        = {std::complex<double>(-0.25, 0.45), std::complex<double>(0.1, -0.1), std::complex<double>(0.3, 0.25)};
    for (int a = 0; a < 3; ++a)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> y(
            1,
            std::vector<std::vector<std::complex<double>>>(2));
        y[0][0].assign(npwk, std::complex<double>(0.0, 0.0));
        y[0][0][ig0] = gam[a];
        y[0][0][igx] = del[a];
        data_.set_pos_resp(a, y);
    }

    // synthetic converged E-field responses dpsi^E,b_{0,0}
    const std::complex<double> mue[3]
        = {std::complex<double>(0.3, 0.2), std::complex<double>(-0.1, 0.4), std::complex<double>(0.25, -0.15)};
    const std::complex<double> nue[3]
        = {std::complex<double>(0.2, -0.35), std::complex<double>(0.45, 0.1), std::complex<double>(-0.2, -0.05)};
    for (int b = 0; b < 3; ++b)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> e(
            1,
            std::vector<std::vector<std::complex<double>>>(2));
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
            const std::complex<double> dot = std::conj(gam[a]) * mue[b] + std::conj(del[a]) * nue[b];
            const double expect = ((a == b) ? 1.0 : 0.0) - 16.0 * ModuleBase::PI / ucell_.omega * wg(0, 0) * dot.real();
            EXPECT_NEAR(eps(a, b), expect, 1.0e-12) << "a=" << a << " b=" << b;
        }
    }
}

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
    const std::complex<double> alpha[3]
        = {std::complex<double>(0.3, 0.2), std::complex<double>(-0.1, 0.4), std::complex<double>(0.25, -0.15)};
    const std::complex<double> beta[3]
        = {std::complex<double>(0.2, -0.35), std::complex<double>(0.45, 0.1), std::complex<double>(-0.2, -0.05)};
    for (int idir = 0; idir < 3; ++idir)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> disp(
            1,
            std::vector<std::vector<std::complex<double>>>(2));
        disp[0][0].assign(npwk, std::complex<double>(0.0, 0.0));
        disp[0][0][ig0] = alpha[idir];
        disp[0][0][igx] = beta[idir];
        data_.set_dpsi_disp(0, idir, disp);
    }

    // synthetic solved position legs Y^a_{0,0} = P_c x_a|psi_0>
    const std::complex<double> gam[3]
        = {std::complex<double>(0.15, -0.3), std::complex<double>(0.4, 0.05), std::complex<double>(-0.35, 0.2)};
    const std::complex<double> del[3]
        = {std::complex<double>(-0.25, 0.45), std::complex<double>(0.1, -0.1), std::complex<double>(0.3, 0.25)};
    for (int a = 0; a < 3; ++a)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> y(
            1,
            std::vector<std::vector<std::complex<double>>>(2));
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
            const std::complex<double> dot = std::conj(alpha[idir]) * gam[a] + std::conj(beta[idir]) * del[a];
            const double expect = ((a == idir) ? zion : 0.0) - 2.0 * wg(0, 0) * dot.real();
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
        ucell_.atoms[0].tau[i] = ucell_.atoms[0].taud[i] * dfpt_test::cubic_edge;
    }
    ucell_.iat2it.resize(3);
    ucell_.iat2ia.resize(3);
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
    const ModuleBase::Vector3<double> klist[1] = {ModuleBase::Vector3<double>(0.25, 0.0, 0.0)};
    kwfc.initgrids(dfpt_test::lat0, latvec_, pw_rho_.nx, pw_rho_.ny, pw_rho_.nz);
    kwfc.initparameters(false, dfpt_test::ecutwfc, 1, klist);
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
