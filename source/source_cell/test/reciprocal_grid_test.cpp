/**
 * @file reciprocal_grid_test.cpp
 * @brief Unit tests for ModuleCell::ReciprocalGrid base class.
 *
 * Covers the spin-free shared functionality: Monkhorst-Pack mesh generation,
 * the Monkhorst-Pack coordinate formula, direct/Cartesian conversion, weight
 * normalization and the star (IBZ) reduction primitive.
 */
#include "gtest/gtest.h"

#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/magnetism.h"
#include "source_cell/pseudo.h"
#include "source_cell/reciprocal_grid.h"
#include "source_cell/unitcell.h"

#include <cmath>
#include <vector>

// Linker stubs: the symmetry library referenced by build_star_ops needs these
// symbols; the real definitions live in the cell_info object library which is
// not linked into this test. Mirror of klist_test.cpp.
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
SepPot::SepPot() {}
SepPot::~SepPot() {}
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
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

/**
 * @brief Minimal concrete subclass exposing the pure-virtual hook.
 */
class TestGrid : public ModuleCell::ReciprocalGrid
{
  public:
    void reduce_by_symmetry(const UnitCell&,
                            const ModuleSymmetry::Symmetry&,
                            bool,
                            std::string&,
                            bool&) override
    {
    }
};

class ReciprocalGridTest : public testing::Test
{
  protected:
    TestGrid grid;
};

TEST_F(ReciprocalGridTest, Construct)
{
    EXPECT_EQ(grid.nks, 0);
    EXPECT_EQ(grid.nkstot, 0);
    EXPECT_EQ(grid.nkstot_full, 0);
    EXPECT_FALSE(grid.kc_done);
    EXPECT_FALSE(grid.kd_done);
    EXPECT_FALSE(grid.is_mp);
}

TEST_F(ReciprocalGridTest, MPFormula)
{
    // k_type=1 (MP, without Gamma)
    EXPECT_DOUBLE_EQ(grid.Monkhorst_Pack_formula(1, 0.0, 1, 4), (0.0 + 2.0 - 4.0 - 1.0) / 8.0);
    EXPECT_DOUBLE_EQ(grid.Monkhorst_Pack_formula(1, 0.5, 2, 4), (0.5 + 4.0 - 4.0 - 1.0) / 8.0);
    // k_type=0 (Gamma-centered)
    EXPECT_DOUBLE_EQ(grid.Monkhorst_Pack_formula(0, 0.0, 1, 4), 0.0);
    EXPECT_DOUBLE_EQ(grid.Monkhorst_Pack_formula(0, 0.0, 4, 4), 3.0 / 4.0);
}

TEST_F(ReciprocalGridTest, MonkhorstPackGeneration)
{
    const int nmp[3] = {2, 3, 4};
    const double offset[3] = {0.0, 0.0, 0.0};
    grid.Monkhorst_Pack(nmp, offset, 0);

    EXPECT_TRUE(grid.is_mp == false); // is_mp is not touched by Monkhorst_Pack
    EXPECT_EQ(grid.nkstot, 24);
    EXPECT_TRUE(grid.kd_done);
    EXPECT_EQ(grid.kvec_d.size(), 24);
    EXPECT_EQ(grid.wk.size(), 24);

    const double weight = 1.0 / 24.0;
    for (int i = 0; i < grid.nkstot; ++i)
    {
        EXPECT_DOUBLE_EQ(grid.wk[i], weight);
    }

    // Gamma-centered: first point is (0,0,0)
    EXPECT_DOUBLE_EQ(grid.kvec_d[0].x, 0.0);
    EXPECT_DOUBLE_EQ(grid.kvec_d[0].y, 0.0);
    EXPECT_DOUBLE_EQ(grid.kvec_d[0].z, 0.0);

    // Last point (x=2,y=3,z=4): (0.5, 2/3, 0.75)
    EXPECT_DOUBLE_EQ(grid.kvec_d[23].x, 0.5);
    EXPECT_DOUBLE_EQ(grid.kvec_d[23].y, 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(grid.kvec_d[23].z, 0.75);
}

TEST_F(ReciprocalGridTest, D2CConversion)
{
    const ModuleBase::Matrix3 G(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    grid.kvec_d.resize(1);
    grid.kvec_d[0] = ModuleBase::Vector3<double>(0.2, 0.3, 0.4);
    grid.kvec_c.resize(1);

    grid.kvec_d2c(G);
    EXPECT_DOUBLE_EQ(grid.kvec_c[0].x, 0.2);
    EXPECT_DOUBLE_EQ(grid.kvec_c[0].y, 0.3);
    EXPECT_DOUBLE_EQ(grid.kvec_c[0].z, 0.4);
}

TEST_F(ReciprocalGridTest, D2CCleansNumericalNoise)
{
    const ModuleBase::Matrix3 G(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    grid.kvec_d.resize(1);
    grid.kvec_d[0] = ModuleBase::Vector3<double>(1.0e-12, 0.5, -1.0e-11);
    grid.kvec_c.resize(1);

    grid.kvec_d2c(G);
    EXPECT_DOUBLE_EQ(grid.kvec_d[0].x, 0.0);
    EXPECT_DOUBLE_EQ(grid.kvec_d[0].z, 0.0);
    EXPECT_DOUBLE_EQ(grid.kvec_c[0].x, 0.0);
    EXPECT_DOUBLE_EQ(grid.kvec_c[0].z, 0.0);
}

TEST_F(ReciprocalGridTest, D2CC2DRoundTrip)
{
    // non-trivial lattice
    const ModuleBase::Matrix3 G(2.0, 0.0, 1.0, 0.0, 3.0, 0.5, 1.0, 0.0, 4.0);
    // kvec_c2d uses R^T; round-trip requires R^T = G^{-1}
    const ModuleBase::Matrix3 R = G.Inverse().Transpose();

    grid.kvec_d.resize(3);
    grid.kvec_d[0] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3);
    grid.kvec_d[1] = ModuleBase::Vector3<double>(0.5, 0.5, 0.5);
    grid.kvec_d[2] = ModuleBase::Vector3<double>(-0.25, 0.75, 0.125);
    grid.kvec_c.resize(3);

    const ModuleBase::Vector3<double> d_ref[3] = {ModuleBase::Vector3<double>(0.1, 0.2, 0.3),
                                                  ModuleBase::Vector3<double>(0.5, 0.5, 0.5),
                                                  ModuleBase::Vector3<double>(-0.25, 0.75, 0.125)};

    grid.kvec_d2c(G);
    grid.kvec_c2d(R);

    for (int i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(grid.kvec_d[i].x, d_ref[i].x, 1e-12);
        EXPECT_NEAR(grid.kvec_d[i].y, d_ref[i].y, 1e-12);
        EXPECT_NEAR(grid.kvec_d[i].z, d_ref[i].z, 1e-12);
    }
}

TEST_F(ReciprocalGridTest, NormalizeWk)
{
    const int nmp[3] = {2, 2, 2};
    const double offset[3] = {0.0, 0.0, 0.0};
    grid.Monkhorst_Pack(nmp, offset, 0);

    grid.normalize_wk(1);
    double sum = 0.0;
    for (int i = 0; i < grid.nkstot; ++i)
    {
        sum += grid.wk[i];
    }
    EXPECT_NEAR(sum, 1.0, 1e-12);

    grid.normalize_wk(2);
    sum = 0.0;
    for (int i = 0; i < grid.nkstot; ++i)
    {
        sum += grid.wk[i];
    }
    EXPECT_NEAR(sum, 2.0, 1e-12);
}

TEST_F(ReciprocalGridTest, ReduceIbzNonMp)
{
    // two points related by inversion: they must fold into one, with the
    // combined weight. This exercises the -q (time-reversal) folding path.
    const ModuleBase::Matrix3 G(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    const ModuleBase::Matrix3 inv(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0);
    const ModuleBase::Matrix3 ind(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    grid.is_mp = false;
    grid.nkstot = 2;
    grid.nkstot_full = 2;
    grid.kvec_d.resize(2);
    grid.kvec_d[0] = ModuleBase::Vector3<double>(0.25, 0.25, 0.25);
    grid.kvec_d[1] = ModuleBase::Vector3<double>(-0.25, -0.25, -0.25);
    grid.wk.resize(2);
    grid.wk[0] = 0.5;
    grid.wk[1] = 0.5;

    ModuleBase::Matrix3 ops[2] = {ind, inv};
    std::vector<ModuleBase::Vector3<double>> vec_ibz;
    std::vector<double> wk_ibz;
    std::vector<int> ibz_index;
    std::vector<int> ibz2bz;
    grid.reduce_ibz(ops, 2, G, G, nullptr, 1e-6, vec_ibz, wk_ibz, ibz_index, ibz2bz);

    EXPECT_EQ(vec_ibz.size(), 1);
    EXPECT_DOUBLE_EQ(wk_ibz[0], 1.0);
    EXPECT_EQ(ibz_index[0], 0);
    EXPECT_EQ(ibz_index[1], 0);
    EXPECT_EQ(ibz2bz[0], 0);
}

TEST_F(ReciprocalGridTest, ReduceIbzKeepsDistinctPoints)
{
    // two points NOT related by any operation: both survive
    const ModuleBase::Matrix3 G(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    const ModuleBase::Matrix3 ind(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    grid.is_mp = false;
    grid.nkstot = 2;
    grid.nkstot_full = 2;
    grid.kvec_d.resize(2);
    grid.kvec_d[0] = ModuleBase::Vector3<double>(0.25, 0.25, 0.25);
    grid.kvec_d[1] = ModuleBase::Vector3<double>(0.50, 0.50, 0.50);
    grid.wk.resize(2);
    grid.wk[0] = 0.5;
    grid.wk[1] = 0.5;

    ModuleBase::Matrix3 ops[1] = {ind};
    std::vector<ModuleBase::Vector3<double>> vec_ibz;
    std::vector<double> wk_ibz;
    std::vector<int> ibz_index;
    std::vector<int> ibz2bz;
    grid.reduce_ibz(ops, 1, G, G, nullptr, 1e-6, vec_ibz, wk_ibz, ibz_index, ibz2bz);

    EXPECT_EQ(vec_ibz.size(), 2);
    EXPECT_DOUBLE_EQ(wk_ibz[0], 0.5);
    EXPECT_DOUBLE_EQ(wk_ibz[1], 0.5);
    EXPECT_EQ(ibz_index[0], 0);
    EXPECT_EQ(ibz_index[1], 1);
}

TEST_F(ReciprocalGridTest, ReduceIbzMpKLattice)
{
    // Monkhorst-Pack path: the {0, 0.5}^3 gamma-centered mesh folded by the
    // closed group {I, C3, C3^2} (order-3 rotations about (1,1,1)) yields
    // Gamma + 3 X + 3 M + R, with the k-lattice consistency asserts active.
    const ModuleBase::Matrix3 G(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    const ModuleBase::Matrix3 ind(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    // row-vector convention: (a,b,c) * c3 = (c,a,b); * c3sq = (b,c,a)
    const ModuleBase::Matrix3 c3(0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
    const ModuleBase::Matrix3 c3sq(0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    const int nmp[3] = {2, 2, 2};
    const double offset[3] = {0.0, 0.0, 0.0};
    grid.Monkhorst_Pack(nmp, offset, 0); // sets nkstot=8, wk=1/8, kd_done
    grid.is_mp = true;
    grid.nkstot_full = grid.nkstot;

    // k-lattice basis of the 2x2x2 mesh: G/2 along each reciprocal axis.
    // In this diagonal frame the k-lattice rotations equal the reciprocal ones.
    const ModuleBase::Matrix3 k_lattice(0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5);
    const ModuleBase::Matrix3 ops[3] = {ind, c3, c3sq};
    const std::vector<ModuleBase::Matrix3> kkmatrix(ops, ops + 3);

    std::vector<ModuleBase::Vector3<double>> vec_ibz;
    std::vector<double> wk_ibz;
    std::vector<int> ibz_index;
    std::vector<int> ibz2bz;
    grid.reduce_ibz(ops, 3, G, k_lattice, kkmatrix.data(), 1e-6, vec_ibz, wk_ibz, ibz_index, ibz2bz);

    ASSERT_EQ(vec_ibz.size(), 4);
    // every mesh point is mapped to an irreducible point
    for (int i = 0; i < grid.nkstot; ++i)
    {
        EXPECT_GE(ibz_index[i], 0);
    }
    // Gamma first, then representatives of the X, M, R stars
    EXPECT_DOUBLE_EQ(vec_ibz[0].x, 0.0);
    EXPECT_DOUBLE_EQ(vec_ibz[0].y, 0.0);
    EXPECT_DOUBLE_EQ(vec_ibz[0].z, 0.0);
    // stars: Gamma(1) + X(3) + M(3) + R(1) -> weights 1/8, 3/8, 3/8, 1/8
    EXPECT_DOUBLE_EQ(wk_ibz[0], 0.125);
    EXPECT_DOUBLE_EQ(wk_ibz[1], 0.375);
    EXPECT_DOUBLE_EQ(wk_ibz[2], 0.375);
    EXPECT_DOUBLE_EQ(wk_ibz[3], 0.125);
    double sum = 0.0;
    for (size_t i = 0; i < wk_ibz.size(); ++i)
    {
        sum += wk_ibz[i];
    }
    EXPECT_NEAR(sum, 1.0, 1e-12);
}
