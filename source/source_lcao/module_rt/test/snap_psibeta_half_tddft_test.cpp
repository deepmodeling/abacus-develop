#include "source_lcao/module_rt/snap_psibeta_half_tddft.h"

#include "source_base/ylm.h"
#include "source_basis/module_nao/radial_collection.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/setup_nonlocal.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <gtest/gtest.h>
#include <vector>

InfoNonlocal::InfoNonlocal()
{
    this->Beta = new Numerical_Nonlocal[1];
    this->nproj = nullptr;
    this->nprojmax = 0;
    this->rcutmax_Beta = 0.0;
}

InfoNonlocal::~InfoNonlocal()
{
    delete[] this->Beta;
    delete[] this->nproj;
}

namespace
{
struct ComparisonStats
{
    double max_real_diff = 0.0;
    double max_imag_abs = 0.0;
    double max_reference_abs = 0.0;
};

class SnapPsibetaHalfTddftTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        ModuleBase::Ylm::set_coefficients();

        const std::string root = "../../../../../";
        const std::string orb_file = "tests/PP_ORB/C_gga_8au_100Ry_2s2p1d.orb";
        const std::string full_orb_file = root + orb_file;
        const std::string orbital_files[1] = {orb_file};

        std::ofstream ofs("snap_psibeta_half_tddft_test.log");
        orb.init(ofs, 1, root, orbital_files, "", 2, 100.0, 0.01, 0.01, 30.0, false, 0, false, false, 0);

        build_fake_beta_projectors();

        orb_radials.build(1, &full_orb_file, 'o');
        beta_radials.build(1, info_nl.Beta);

        const double rmax = std::max(orb_radials.rcut_max(), beta_radials.rcut_max());
        const double cutoff = 2.0 * rmax;
        const int nr = static_cast<int>(rmax / 0.01) + 1;

        orb_radials.set_uniform_grid(true, nr, cutoff, 'i', true);
        beta_radials.set_uniform_grid(true, nr, cutoff, 'i', true);
        overlap_orb_beta.tabulate(orb_radials, beta_radials, 'S', nr, cutoff);
    }

    void build_fake_beta_projectors()
    {
        const int nproj = 2;
        std::vector<Numerical_Nonlocal_Lm> beta_lm(nproj);

        for (int iproj = 0; iproj < nproj; ++iproj)
        {
            const int l = iproj;
            const auto& phi_ln = orb.Phi[0].PhiLN(l, 0);
            beta_lm[iproj].set_NL_proj("C",
                                       0,
                                       l,
                                       phi_ln.getNr(),
                                       phi_ln.getRab(),
                                       phi_ln.getRadial(),
                                       phi_ln.getPsi_r(),
                                       orb.get_kmesh(),
                                       orb.get_dk(),
                                       orb.get_dr_uniform());
        }

        info_nl.nproj = new int[1];
        info_nl.nproj[0] = nproj;
        info_nl.nprojmax = nproj;
        info_nl.Beta[0].set_type_info(0, "C", "NC", 1, nproj, beta_lm.data());
        info_nl.rcutmax_Beta = info_nl.Beta[0].get_rcut_max();
    }

    static int abacus_m_to_m(const int m)
    {
        return (m % 2 == 0) ? -m / 2 : (m + 1) / 2;
    }

    ComparisonStats compare_zero_vector_potential(const int lebedev_grid_points)
    {
        const ModuleBase::Vector3<double> R0(0.1, -0.2, 0.3);
        const ModuleBase::Vector3<double> R1(0.4, 0.2, -0.1);
        const ModuleBase::Vector3<double> zero_A(0.0, 0.0, 0.0);
        module_rt::SnapIntegrationOptions options;
        options.lebedev_grid_points = lebedev_grid_points;

        ComparisonStats stats;

        for (int L1 = 0; L1 <= orb.Phi[0].getLmax(); ++L1)
        {
            for (int N1 = 0; N1 < orb.Phi[0].getNchi(L1); ++N1)
            {
                for (int m1 = 0; m1 < 2 * L1 + 1; ++m1)
                {
                    std::vector<std::vector<std::complex<double>>> grid_nlm;
                    module_rt::snap_psibeta_half_tddft(orb,
                                                       info_nl,
                                                       grid_nlm,
                                                       R1,
                                                       0,
                                                       L1,
                                                       m1,
                                                       N1,
                                                       R0,
                                                       0,
                                                       zero_A,
                                                       false,
                                                       options);

                    std::vector<std::vector<double>> tci_nlm;
                    overlap_orb_beta.snap(0, L1, N1, abacus_m_to_m(m1), 0, R0 - R1, false, tci_nlm);

                    EXPECT_FALSE(grid_nlm.empty());
                    EXPECT_FALSE(tci_nlm.empty());
                    if (grid_nlm.empty() || tci_nlm.empty())
                    {
                        continue;
                    }
                    EXPECT_EQ(grid_nlm[0].size(), tci_nlm[0].size());
                    if (grid_nlm[0].size() != tci_nlm[0].size())
                    {
                        continue;
                    }

                    for (size_t i = 0; i < grid_nlm[0].size(); ++i)
                    {
                        stats.max_real_diff
                            = std::max(stats.max_real_diff, std::abs(grid_nlm[0][i].real() - tci_nlm[0][i]));
                        stats.max_imag_abs = std::max(stats.max_imag_abs, std::abs(grid_nlm[0][i].imag()));
                        stats.max_reference_abs = std::max(stats.max_reference_abs, std::abs(tci_nlm[0][i]));
                    }
                }
            }
        }

        return stats;
    }

    LCAO_Orbitals orb;
    InfoNonlocal info_nl;
    RadialCollection orb_radials;
    RadialCollection beta_radials;
    TwoCenterIntegrator overlap_orb_beta;
};
} // namespace

TEST_F(SnapPsibetaHalfTddftTest, ZeroVectorPotentialMatchesTwoCenterIntegral)
{
    const double real_tolerance = 5.0e-8;
    const double imag_tolerance = 1.0e-12;
    const ComparisonStats stats = compare_zero_vector_potential(110);

    EXPECT_LT(stats.max_real_diff, real_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_imag_abs, imag_tolerance) << "max reference abs = " << stats.max_reference_abs;
}

TEST_F(SnapPsibetaHalfTddftTest, ZeroVectorPotentialHighOrderGridMatchesTwoCenterIntegral)
{
    const double real_tolerance = 5.0e-8;
    const double imag_tolerance = 1.0e-12;
    const ComparisonStats stats = compare_zero_vector_potential(590);

    EXPECT_LT(stats.max_real_diff, real_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_imag_abs, imag_tolerance) << "max reference abs = " << stats.max_reference_abs;
}

TEST_F(SnapPsibetaHalfTddftTest, ZeroVectorPotentialPositionMomentsAreReal)
{
    const ModuleBase::Vector3<double> R0(-0.3, 0.2, 0.1);
    const ModuleBase::Vector3<double> R1(0.2, -0.1, 0.4);
    const ModuleBase::Vector3<double> zero_A(0.0, 0.0, 0.0);
    const double tolerance = 1.0e-12;

    std::vector<std::vector<std::complex<double>>> nlm;
    module_rt::snap_psibeta_half_tddft(orb, info_nl, nlm, R1, 0, 1, 1, 0, R0, 0, zero_A, true);

    ASSERT_EQ(nlm.size(), 4);
    for (const auto& dim: nlm)
    {
        ASSERT_EQ(dim.size(), 4);
        for (const std::complex<double>& value: dim)
        {
            EXPECT_NEAR(value.imag(), 0.0, tolerance);
        }
    }
}
