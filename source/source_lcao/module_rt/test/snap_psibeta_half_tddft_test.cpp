#include "source_lcao/module_rt/snap_psibeta_half_tddft.h"

#include "source_base/ylm.h"
#include "source_cell/read_pp.h"
#include "source_cell/setup_nonlocal.h"
#include "source_cell/unitcell.h"
#include "source_io/module_hs/cal_r_overlap_R.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <gtest/gtest.h>
#include <vector>

SepPot::SepPot() = default;

SepPot::~SepPot()
{
    delete[] r;
    delete[] rv;
}

Sep_Cell::Sep_Cell() noexcept : ntype(0), omega(0.0), tpiba2(0.0)
{
}

Sep_Cell::~Sep_Cell() noexcept = default;

Magnetism::Magnetism()
{
    tot_mag = 0.0;
    abs_mag = 0.0;
}

Magnetism::~Magnetism()
{
    delete[] start_mag;
}

UnitCell::UnitCell()
{
    itia2iat.create(1, 1);
}

UnitCell::~UnitCell()
{
    if (set_atom_flag)
    {
        delete[] atoms;
    }
}

namespace
{
struct ComparisonStats
{
    double max_overlap_diff = 0.0;
    double max_position_diff = 0.0;
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
        const std::string orb_file = "tests/PP_ORB/Ti_gga_10au_100Ry_4s2p2d1f.orb";
        const std::string orbital_files[1] = {orb_file};

        std::ofstream ofs("snap_psibeta_half_tddft_test.log");
        orb.init(ofs, 1, root, orbital_files, "", 3, 100.0, 0.01, 0.01, 30.0, false, 0, false, false, 0);

        ASSERT_EQ(orb.Phi[0].getLmax(), 3);
        ASSERT_EQ(orb.Phi[0].getNchi(0), 4);
        ASSERT_EQ(orb.Phi[0].getNchi(1), 2);
        ASSERT_EQ(orb.Phi[0].getNchi(2), 2);
        ASSERT_EQ(orb.Phi[0].getNchi(3), 1);

        build_ti_beta_projectors(root);
        initialize_r_overlap_reference();
    }

    void build_ti_beta_projectors(const std::string& root)
    {
        ucell.ntype = 1;
        ucell.nat = 1;
        ucell.atoms = new Atom[1];
        ucell.set_atom_flag = true;

        Atom& atom = ucell.atoms[0];
        atom.label = "Ti";
        atom.type = 0;
        atom.na = 1;
        atom.nwl = orb.Phi[0].getLmax();
        atom.l_nchi.resize(atom.nwl + 1);
        atom.nw = 0;
        for (int L = 0; L <= atom.nwl; ++L)
        {
            atom.l_nchi[L] = orb.Phi[0].getNchi(L);
            atom.nw += (2 * L + 1) * atom.l_nchi[L];
        }
        atom.tau.resize(1);
        atom.tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);

        Pseudopot_upf pseudo_reader;
        std::string pseudo_type = "auto";
        const int pseudo_error = pseudo_reader.init_pseudo_reader(root + "tests/PP_ORB/Ti_ONCV_PBE-1.0.upf", pseudo_type, atom.ncpp);
        ASSERT_EQ(pseudo_error, 0);
        ASSERT_EQ(pseudo_type, "upf201");
        ASSERT_EQ(atom.ncpp.psd, "Ti");
        ASSERT_EQ(atom.ncpp.pp_type, "NC");
        ASSERT_EQ(atom.ncpp.nbeta, 6);
        ASSERT_EQ(atom.ncpp.lll, std::vector<int>({0, 0, 1, 1, 2, 2}));
        pseudo_reader.complete_default(atom.ncpp);
        ASSERT_EQ(atom.ncpp.nh, 18);
        ASSERT_EQ(atom.ncpp.jjj.size(), 6);

        ucell.infoNL.nproj = new int[1];
        std::ofstream log("snap_psibeta_half_tddft_nonlocal.log");
        ucell.infoNL.Set_NonLocal(0, &atom, ucell.infoNL.nproj[0], orb.get_kmesh(), orb.get_dk(), orb.get_dr_uniform(), log);

        ASSERT_EQ(ucell.infoNL.nproj[0], 6);
        ucell.infoNL.nprojmax = ucell.infoNL.nproj[0];
        ucell.infoNL.rcutmax_Beta = ucell.infoNL.Beta[0].get_rcut_max();
    }

    void initialize_r_overlap_reference()
    {
        r_calculator.init_nonlocal(ucell, pv, orb);
    }

    ComparisonStats compare_zero_vector_potential(const int radial_grid_num, const int lebedev_grid_points)
    {
        const ModuleBase::Vector3<double> R0(0.1, -0.2, 0.3);
        const ModuleBase::Vector3<double> R1(0.4, 0.2, -0.1);
        const ModuleBase::Vector3<double> zero_A(0.0, 0.0, 0.0);
        module_rt::SnapIntegrationOptions options;
        options.radial_grid_num = radial_grid_num;
        options.lebedev_grid_points = lebedev_grid_points;

        ComparisonStats stats;

        for (int L1 = 0; L1 <= orb.Phi[0].getLmax(); ++L1)
        {
            for (int N1 = 0; N1 < orb.Phi[0].getNchi(L1); ++N1)
            {
                for (int m1 = 0; m1 < 2 * L1 + 1; ++m1)
                {
                    std::vector<std::vector<std::complex<double>>> grid_nlm;
                    module_rt::snap_psibeta_half_tddft(orb, ucell.infoNL, grid_nlm, R1, 0, L1, m1, N1, R0, 0, zero_A, true, options);

                    std::vector<std::vector<double>> reference_nlm;
                    r_calculator.get_psi_r_beta(ucell, reference_nlm, R1, 0, L1, m1, N1, R0, 0);

                    EXPECT_EQ(grid_nlm.size(), 4);
                    EXPECT_EQ(reference_nlm.size(), 4);
                    if (grid_nlm.size() != 4 || reference_nlm.size() != 4)
                    {
                        continue;
                    }

                    bool sizes_match = true;
                    for (size_t dim = 0; dim < grid_nlm.size(); ++dim)
                    {
                        EXPECT_EQ(grid_nlm[dim].size(), reference_nlm[dim].size());
                        sizes_match = sizes_match && (grid_nlm[dim].size() == reference_nlm[dim].size());
                    }
                    if (!sizes_match)
                    {
                        continue;
                    }

                    for (size_t dim = 0; dim < grid_nlm.size(); ++dim)
                    {
                        for (size_t i = 0; i < grid_nlm[dim].size(); ++i)
                        {
                            const double real_diff = std::abs(grid_nlm[dim][i].real() - reference_nlm[dim][i]);
                            if (dim == 0)
                            {
                                stats.max_overlap_diff = std::max(stats.max_overlap_diff, real_diff);
                            }
                            else
                            {
                                stats.max_position_diff = std::max(stats.max_position_diff, real_diff);
                            }
                            stats.max_imag_abs = std::max(stats.max_imag_abs, std::abs(grid_nlm[dim][i].imag()));
                            stats.max_reference_abs = std::max(stats.max_reference_abs, std::abs(reference_nlm[dim][i]));
                        }
                    }
                }
            }
        }

        return stats;
    }

    LCAO_Orbitals orb;
    UnitCell ucell;
    Parallel_Orbitals pv;
    cal_r_overlap_R r_calculator;
};
} // namespace

TEST_F(SnapPsibetaHalfTddftTest, ZeroVectorPotentialMatchesTwoCenterIntegral)
{
    const double overlap_tolerance = 4.0e-7;
    const double position_tolerance = 6.0e-7;
    const double imag_tolerance = 1.0e-12;
    const ComparisonStats stats = compare_zero_vector_potential(140, 110);

    EXPECT_LT(stats.max_overlap_diff, overlap_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_position_diff, position_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_imag_abs, imag_tolerance) << "max reference abs = " << stats.max_reference_abs;
}

TEST_F(SnapPsibetaHalfTddftTest, ZeroVectorPotentialDenseRadialGridMatchesTwoCenterIntegral)
{
    const double overlap_tolerance = 3.0e-7;
    const double position_tolerance = 5.0e-7;
    const double imag_tolerance = 1.0e-12;
    const ComparisonStats stats = compare_zero_vector_potential(280, 110);

    EXPECT_LT(stats.max_overlap_diff, overlap_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_position_diff, position_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_imag_abs, imag_tolerance) << "max reference abs = " << stats.max_reference_abs;
}

TEST_F(SnapPsibetaHalfTddftTest, ZeroVectorPotentialHighOrderGridMatchesTwoCenterIntegral)
{
    const double overlap_tolerance = 4.0e-7;
    const double position_tolerance = 6.0e-7;
    const double imag_tolerance = 1.0e-12;
    const ComparisonStats stats = compare_zero_vector_potential(140, 590);

    EXPECT_LT(stats.max_overlap_diff, overlap_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_position_diff, position_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_imag_abs, imag_tolerance) << "max reference abs = " << stats.max_reference_abs;
}
