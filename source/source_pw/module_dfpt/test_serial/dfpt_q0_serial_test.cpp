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
#include "source_lcao/module_dftu/dftu.h"
#include "source_psi/psi.h"

// test-support ctor/dtor stubs (see dfpt_pert_serial_test.cpp)
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
 *   - DFPT_Q0::compute_eps: the full prefactor chain (8 pi / Omega, wg,
 *     1/(eps_c - eps_v)) on a two-level toy system with a complex excited
 *     state.
 *   - DFPT_Q0::compute_born: elementwise against the closed-form
 *     <u_v|dV/dtau_b|u_m><u_m|r_a|u_v> sums at q = 0 (Coulomb local part),
 *     the ionic Z delta_ab, and the dpsi-slot backup/restore.
 */

class DFPTQ0SerialTest : public testing::Test
{
  protected:
    const double lat0_ = 1.8897261254578281;
    const double ecutwfc_ = 2.5;
    const double rho_mult_ = 9.0;
    const double a_ = 10.0; // cubic edge in lat0 units

    ModuleBase::Matrix3 latvec_;
    UnitCell ucell_;
    ModulePW::PW_Basis pw_rho_;
    ModulePW::PW_Basis_K pw_wfc_;
    Structure_Factor sf_;
    ModuleDFPT::DFPT_Pert pert_;
    ModuleDFPT::DFPT_Q0 q0_;
    ModuleCell::QList qlist_;
    ModuleDFPT::DFPT_PW_Data data_;

    const ModuleBase::Vector3<double> k_d_{0.0, 0.0, 0.0}; // Gamma only
    const ModuleBase::Vector3<double> tau_{1.1, 2.3, 0.7};  // lat0 units
    const ModuleBase::Vector3<double> gx_{0.1, 0.0, 0.0};   // 1/lat0 units
    const ModuleBase::Vector3<double> gy_{0.0, 0.1, 0.0};

    void SetUp() override
    {
        latvec_ = ModuleBase::Matrix3(a_, 0.0, 0.0, 0.0, a_, 0.0, 0.0, 0.0, a_);
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
        ucell_.omega = a_ * a_ * a_ * lat0_ * lat0_ * lat0_;
        ucell_.iat2it = new int[1];
        ucell_.iat2ia = new int[1];
        ucell_.iat2it[0] = 0;
        ucell_.iat2ia[0] = 0;
        MakeCoulombAtom();

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
        qlist_.kvec_d.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0));

        data_.init(&qlist_, 1, 4, pw_wfc_.npwk_max, pw_rho_.nrxx, 1, 1, nullptr);
        pert_.init(ucell_, &pw_rho_, &pw_wfc_, sf_);
        q0_.init(ucell_, &pw_rho_, &pw_wfc_, &pert_);
    }

    void TearDown() override
    {
        delete[] ucell_.atoms;
        ucell_.atoms = nullptr;
        delete[] ucell_.iat2it;
        ucell_.iat2it = nullptr;
        delete[] ucell_.iat2ia;
        ucell_.iat2ia = nullptr;
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
        at.mass = 12.0;
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
            p.betar(0, i) = std::exp(-std::pow(r - 1.0, 2.0) / (2.0 * 0.3 * 0.3));
            p.betar(1, i) = std::exp(-std::pow(r - 1.2, 2.0) / (2.0 * 0.35 * 0.35));
        }
        p.dion.create(2, 2);
        p.dion(0, 0) = 0.8;
        p.dion(0, 1) = 0.15;
        p.dion(1, 0) = -0.25;
        p.dion(1, 1) = 1.1;
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

    double VlocCoulomb(double g2_bohr) const
    {
        return -ucell_.atoms[0].ncpp.zv * ModuleBase::e2 * ModuleBase::FOUR_PI / ucell_.omega
               / g2_bohr;
    }

    // analytic dVloc/dtau_dir coefficient at displacement vector w (1/lat0),
    // GS structure-factor phase convention exp(-i 2pi w.tau)
    std::complex<double> AnalyticDVloc(int dir, const ModuleBase::Vector3<double>& w) const
    {
        const double w2 = w * w;
        if (w2 < 1.0e-12)
        {
            return std::complex<double>(0.0, 0.0);
        }
        const double arg = -ModuleBase::TWO_PI * (w * tau_);
        return std::complex<double>(0.0, -1.0) * (ucell_.tpiba * w[dir])
               * VlocCoulomb(w2 * ucell_.tpiba2)
               * std::complex<double>(std::cos(arg), std::sin(arg));
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
// compute_eps on a two-level system with a complex excited state
// ---------------------------------------------------------------------------

TEST_F(DFPTQ0SerialTest, ComputeEpsTwoLevelAnalytic)
{
    const int npwk = pw_wfc_.npwk[0];
    const int ig0 = IgOf(0, 0, 0);
    const int igx = IgOf(1, 0, 0);
    const int igy = IgOf(0, 1, 0);
    ASSERT_GE(ig0, 0);
    ASSERT_GE(igx, 0);
    ASSERT_GE(igy, 0);

    const double e1 = ucell_.tpiba2 * (gx_ * gx_);
    // v = sqrt(0.6)|G0> + sqrt(0.4)|Gx>          (eps = 0.4 e1)
    // c = sqrt(0.2)|G0> - sqrt(0.3)|Gx> + i sqrt(0.5)|Gy>  (eps = e1)
    // (the relative i phase keeps the reference sensitive to conj placement)
    psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
    psi.zero_out();
    psi(0, 0, ig0) = std::sqrt(0.6);
    psi(0, 0, igx) = std::sqrt(0.4);
    psi(0, 1, ig0) = std::sqrt(0.2);
    psi(0, 1, igx) = -std::sqrt(0.3);
    psi(0, 1, igy) = std::complex<double>(0.0, std::sqrt(0.5));

    ModuleBase::matrix wg(1, 2);
    wg(0, 0) = 2.0;
    wg(0, 1) = 0.0;
    ModuleBase::matrix eig(1, 2);
    eig(0, 0) = 0.4 * e1;
    eig(0, 1) = e1;

    q0_.compute_eps(psi, wg, eig, data_);
    const ModuleBase::matrix eps = data_.get_dielectric();

    // closed-form velocity elements between the two states (kinetic operator
    // diagonal in G: G0 pairs G0 with G = 0, the Gy component of c pairs
    // with the vanishing Gy component of v)
    std::complex<double> p01[3];
    for (int d = 0; d < 3; ++d)
    {
        p01[d] = 2.0 * ucell_.tpiba2
                 * (std::sqrt(0.6) * std::sqrt(0.2) * 0.0
                    + std::sqrt(0.4) * (-std::sqrt(0.3)) * gx_[d]);
    }
    const double de = eig(0, 1) - eig(0, 0);
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            const std::complex<double> r_vc
                = std::complex<double>(0.0, -1.0) * p01[a] / (ucell_.tpiba * (eig(0, 0) - eig(0, 1)));
            const std::complex<double> r_cv
                = std::complex<double>(0.0, -1.0) * std::conj(p01[b])
                      / (ucell_.tpiba * (eig(0, 1) - eig(0, 0)));
            const double expect = ((a == b) ? 1.0 : 0.0)
                                  + 8.0 * ModuleBase::PI / ucell_.omega * wg(0, 0)
                                        * (r_vc * r_cv).real() / de;
            EXPECT_NEAR(eps(a, b), expect, 1.0e-10) << "a=" << a << " b=" << b;
        }
    }
}

// ---------------------------------------------------------------------------
// compute_born against the closed-form q = 0 sums (Coulomb local dV)
// ---------------------------------------------------------------------------

TEST_F(DFPTQ0SerialTest, ComputeBornAnalyticCoulomb)
{
    const int npwk = pw_wfc_.npwk[0];
    const int ig0 = IgOf(0, 0, 0);
    const int igx = IgOf(1, 0, 0);
    const int igy = IgOf(0, 1, 0);
    ASSERT_GE(ig0, 0);
    ASSERT_GE(igx, 0);
    ASSERT_GE(igy, 0);

    const double e1 = ucell_.tpiba2 * (gx_ * gx_);
    psi::Psi<std::complex<double>> psi(1, 2, npwk, npwk, true);
    psi.zero_out();
    psi(0, 0, ig0) = std::sqrt(0.6);
    psi(0, 0, igx) = std::sqrt(0.4);
    psi(0, 1, ig0) = std::sqrt(0.2);
    psi(0, 1, igx) = -std::sqrt(0.3);
    psi(0, 1, igy) = std::complex<double>(0.0, std::sqrt(0.5));

    ModuleBase::matrix wg(1, 2);
    wg(0, 0) = 2.0;
    wg(0, 1) = 0.0;
    ModuleBase::matrix eig(1, 2);
    eig(0, 0) = 0.4 * e1;
    eig(0, 1) = e1;

    // sentinel dpsi in the q = 0 slot: compute_born must restore it
    const std::vector<std::complex<double>> sentinel(npwk, std::complex<double>(0.5, -0.25));
    data_.set_dpsi(0, 0, 0, sentinel);

    q0_.compute_born(psi, wg, eig, data_);
    const ModuleBase::matrix zstar = data_.get_born(0);

    // closed-form dV matrix elements: supp(v) = {G0, Gx}, supp(m) = {G0, Gx, Gy}
    const std::complex<double> cv[3] = {psi(0, 0, ig0), psi(0, 0, igx), std::complex<double>(0.0, 0.0)};
    const ModuleBase::Vector3<double> gv[3] = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0), gx_,
                                               gy_};
    const std::complex<double> cm[3] = {psi(0, 1, ig0), psi(0, 1, igx), psi(0, 1, igy)};
    const double de = eig(0, 1) - eig(0, 0);
    for (int idir = 0; idir < 3; ++idir)
    {
        // dv_m0 = <u_m|dV/dtau_idir|u_v> = sum_{G'' in supp(m)} cc_m(G'')
        //         sum_{G' in supp(v)} c_v(G') AnalyticDVloc(G'' - G')
        std::complex<double> dv10(0.0, 0.0);
        for (int im = 0; im < 3; ++im)
        {
            if (cm[im] == std::complex<double>(0.0, 0.0))
            {
                continue;
            }
            for (int iv = 0; iv < 2; ++iv)
            {
                dv10 += std::conj(cm[im]) * cv[iv] * AnalyticDVloc(idir, gv[im] - gv[iv]);
            }
        }
        for (int a = 0; a < 3; ++a)
        {
            // p_01^a = 2 tpiba^2 <v|(k+G)_a|m>, kinetic operator diagonal in G:
            // only shared components pair (G0 = 0 drops out, v has no Gy)
            std::complex<double> p01_a(0.0, 0.0);
            for (int g = 0; g < 3; ++g)
            {
                p01_a += 2.0 * ucell_.tpiba2 * std::conj(cv[g]) * gv[g][a] * cm[g];
            }
            const std::complex<double> r_10
                = std::complex<double>(0.0, -1.0) * std::conj(p01_a)
                      / (ucell_.tpiba * (eig(0, 1) - eig(0, 0)));
            // ionic Z sits on the (a == idir) diagonal only
            const double zion = (a == idir) ? ucell_.atoms[0].ncpp.zv : 0.0;
            const double expect
                = zion - 4.0 * wg(0, 0) * (std::conj(dv10) * r_10).real() / de;
            EXPECT_NEAR(zstar(a, idir), expect, 1.0e-9) << "a=" << a << " idir=" << idir;
        }
    }

    // the q = 0 dpsi slot is restored
    const std::vector<std::complex<double>> after = data_.get_dpsi(0, 0, 0);
    ASSERT_EQ(after.size(), sentinel.size());
    for (size_t i = 0; i < after.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(after[i].real(), sentinel[i].real());
        EXPECT_DOUBLE_EQ(after[i].imag(), sentinel[i].imag());
    }
}
