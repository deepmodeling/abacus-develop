// Serial unit tests of DFPT_Q0::pos_matrix (dfpt_q0_pos.cpp): the
// position operator in the velocity (commutator) form. The kinetic
// velocity term, the -i commutator factor, the tpiba scaling and the
// Hermitian structure against closed-form plane-wave-combination states
// (exactly degenerate pairs are skipped); the nonlocal contraction
// against an operator finite difference of <u_m|V_nl(k)|u_n>.

#include "dfpt_serial_fixture.h"
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

// per-file fixture (self-contained, gtest style): this fixture is
// Gamma-only (k_d_ shadows the generic default k of the base)
class DFPTQ0SerialTest : public DFPTSerialBase
{
  protected:
    Structure_Factor sf_;
    ModuleDFPT::DFPT_Pert pert_;
    ModuleDFPT::DFPT_Q0 q0_;

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
            if (std::llround(g.x * a_) == ix && std::llround(g.y * a_) == iy && std::llround(g.z * a_) == iz)
            {
                return ig;
            }
        }
        return -1;
    }
};

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>

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
    auto rnd = [&]() {
        seed = seed * 1664525u + 1013904223u;
        return ((seed >> 8) & 0xffffff) / 16777216.0 * 2.0 - 1.0;
    };
    std::vector<std::vector<std::complex<double>>> c(nb, std::vector<std::complex<double>>(npwk));
    for (int b = 0; b < nb; ++b)
    {
        for (int ig = 0; ig < npwk; ++ig)
        {
            c[b][ig] = (ig == ig0) ? std::complex<double>(0.0, 0.0) : std::complex<double>(rnd(), rnd());
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
                          std::vector<std::vector<std::complex<double>>>& mmat) {
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
                const std::complex<double> p_r = std::complex<double>(0.0, 1.0) * ucell_.tpiba * de * r_mat[0][m][n][d];
                // analytic kinetic + finite-difference nonlocal
                std::complex<double> p_kin(0.0, 0.0);
                for (int ig = 0; ig < npwk; ++ig)
                {
                    p_kin += 2.0 * ucell_.tpiba2 * gk[ig][d] * std::conj(psi(0, m, ig)) * psi(0, n, ig);
                }
                const std::complex<double> p_nl = (mm_p[m][n] - mm_m[m][n]) / (2.0 * eps);
                const double scale = std::max(1.0, std::abs(p_kin) + std::abs(p_nl));
                EXPECT_NEAR(p_r.real(), (p_kin + p_nl).real(), 1.0e-6 * scale) << "m=" << m << " n=" << n << " d=" << d;
                EXPECT_NEAR(p_r.imag(), (p_kin + p_nl).imag(), 1.0e-6 * scale) << "m=" << m << " n=" << n << " d=" << d;
            }
        }
    }
}
