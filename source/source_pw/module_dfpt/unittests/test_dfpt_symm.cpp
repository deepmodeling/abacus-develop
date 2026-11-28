#include "source_pw/module_dfpt/dfpt_symm.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{

/// fcc primitive lattice with cubic side a = 1 (lat0 = 1 units), the
/// diamond-structure cell of the production Si cases
ModuleBase::Matrix3 fcc_latvec() { return ModuleBase::Matrix3(0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0); }

/// the 48 signed-permutation cartesian rotations (the Oh point group)
std::vector<ModuleBase::Matrix3> cubic_cartesian_ops()
{
    std::vector<ModuleBase::Matrix3> ops;
    const int rows[6][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}};
    for (int p = 0; p < 6; ++p)
    {
        for (int sx = 1; sx >= -1; sx -= 2)
        {
            for (int sy = 1; sy >= -1; sy -= 2)
            {
                for (int sz = 1; sz >= -1; sz -= 2)
                {
                    const int signs[3] = {sx, sy, sz};
                    double m[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                    for (int i = 0; i < 3; ++i)
                    {
                        m[i * 3 + rows[p][i]] = signs[i];
                    }
                    ops.push_back(ModuleBase::Matrix3(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]));
                }
            }
        }
    }
    return ops;
}

/// the Oh rotations expressed in the fcc primitive fractional basis
/// (integer matrices, the gmatrix convention of the symmetry module)
std::vector<ModuleBase::Matrix3> fcc_gmatrices()
{
    const ModuleBase::Matrix3 lat = fcc_latvec();
    const ModuleBase::Matrix3 inv_lat = lat.Inverse();
    const std::vector<ModuleBase::Matrix3> ops = cubic_cartesian_ops();
    std::vector<ModuleBase::Matrix3> gm;
    for (const ModuleBase::Matrix3 &r : ops)
    {
        gm.push_back(lat * r * inv_lat);
    }
    return gm;
}

} // namespace

TEST(DFPTLittleGroup, GammaKeepsAllOperations)
{
    const std::vector<ModuleBase::Matrix3> gm = fcc_gmatrices();
    const ModuleBase::Vector3<double> q(0.0, 0.0, 0.0);
    const std::vector<int> group = ModuleDFPT::little_group_of_q(fcc_latvec(), gm.data(), gm.size(), q);
    EXPECT_EQ(group.size(), gm.size());
}

TEST(DFPTLittleGroup, LPointLittleGroupIsD3dOfTwelve)
{
    // L = (pi/a)(-1, 1, 1) with a = 1: the little group including time
    // reversal is D3d (12 operations); the dot-product reciprocal test
    // must reject C2y, whose cartesian difference (1, 0, -1)*2pi/a is
    // integral but not a reciprocal-lattice vector
    const std::vector<ModuleBase::Matrix3> gm = fcc_gmatrices();
    const ModuleBase::Vector3<double> q(-M_PI, M_PI, M_PI);
    const std::vector<int> group = ModuleDFPT::little_group_of_q(fcc_latvec(), gm.data(), gm.size(), q);
    EXPECT_EQ(group.size(), 12u);
    // the identity operation is always a member
    bool has_identity = false;
    for (const int iop : group)
    {
        if (std::abs(gm[iop].e11 - 1.0) < 1.0e-12 && std::abs(gm[iop].e22 - 1.0) < 1.0e-12 &&
            std::abs(gm[iop].e33 - 1.0) < 1.0e-12)
        {
            has_identity = true;
        }
    }
    EXPECT_TRUE(has_identity);
}

TEST(DFPTLittleGroup, XPointLittleGroupIsD4hOfSixteen)
{
    // X = (0, 2pi/a, 0) with a = 1: the little group including time
    // reversal is D4h (16 operations: the maps with the y column +-y)
    const std::vector<ModuleBase::Matrix3> gm = fcc_gmatrices();
    const ModuleBase::Vector3<double> q(0.0, 2.0 * M_PI, 0.0);
    const std::vector<int> group = ModuleDFPT::little_group_of_q(fcc_latvec(), gm.data(), gm.size(), q);
    EXPECT_EQ(group.size(), 16u);
}

TEST(DFPTLittleGroup, GenericQKeepsIdentityAndInversionTimesTR)
{
    // for a generic q the unitary branch keeps only the identity, but the
    // antiunitary branch always admits inversion combined with time
    // reversal ((theta*i) q = q for every q), so the gray little group of
    // any q has at least two members
    const std::vector<ModuleBase::Matrix3> gm = fcc_gmatrices();
    const ModuleBase::Vector3<double> q(1.234, -0.567, 0.891);
    const std::vector<int> group = ModuleDFPT::little_group_of_q(fcc_latvec(), gm.data(), gm.size(), q);
    EXPECT_EQ(group.size(), 2u);
    EXPECT_EQ(group.front(), 0);
}

TEST(DFPTLittleGroup, EmptyInputIsEmptyOutput)
{
    const ModuleBase::Vector3<double> q(0.0, 0.0, 0.0);
    const std::vector<int> group = ModuleDFPT::little_group_of_q(fcc_latvec(), nullptr, 0, q);
    EXPECT_TRUE(group.empty());
}

namespace
{

/// the 2x2x2 Gamma-centered mesh in fractional coordinates (the Si cases)
std::vector<ModuleBase::Vector3<double>> mesh_2x2x2()
{
    std::vector<ModuleBase::Vector3<double>> mesh;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                mesh.push_back(ModuleBase::Vector3<double>(i * 0.5, j * 0.5, k * 0.5));
            }
        }
    }
    return mesh;
}

/// sorted orbit sizes of a pair list, plus the total (must equal the mesh)
std::vector<int> orbit_sizes(const std::vector<ModuleDFPT::KQPair> &pairs)
{
    std::vector<int> sizes;
    for (const ModuleDFPT::KQPair &p : pairs)
    {
        sizes.push_back(p.orbit_size);
    }
    std::sort(sizes.begin(), sizes.end(), std::greater<int>());
    return sizes;
}

} // namespace

TEST(DFPTKqPairs, GammaReproducesTheStandardIbzReduction)
{
    // q = 0 with the full group: the pairs degenerate to (k, k) and the
    // orbit reduction must reproduce the production k-mesh IBZ of the
    // 2x2x2 fcc mesh: Gamma (1) + L star (4) + X star (3)
    const std::vector<ModuleBase::Matrix3> gm = fcc_gmatrices();
    const ModuleBase::Vector3<double> q(0.0, 0.0, 0.0);
    const std::vector<ModuleDFPT::KQPair> pairs =
        ModuleDFPT::kq_pair_reps(fcc_latvec(), gm.data(), gm.size(), q, mesh_2x2x2());
    ASSERT_EQ(pairs.size(), 3u);
    const std::vector<int> sizes = orbit_sizes(pairs);
    EXPECT_EQ(sizes[0], 4);
    EXPECT_EQ(sizes[1], 3);
    EXPECT_EQ(sizes[2], 1);
    // the Gamma representative is present with the singleton orbit
    bool has_gamma = false;
    for (const ModuleDFPT::KQPair &p : pairs)
    {
        if (p.k_frac.x == 0.0 && p.k_frac.y == 0.0 && p.k_frac.z == 0.0)
        {
            has_gamma = true;
            EXPECT_EQ(p.orbit_size, 1);
        }
        // kq = k at q = 0
        EXPECT_DOUBLE_EQ(p.kq_frac.x, p.k_frac.x);
        EXPECT_DOUBLE_EQ(p.kq_frac.y, p.k_frac.y);
        EXPECT_DOUBLE_EQ(p.kq_frac.z, p.k_frac.z);
    }
    EXPECT_TRUE(has_gamma);
}

TEST(DFPTKqPairs, LPointReducesUnderTheD3dUnitaryGroup)
{
    // q = L = (1/2, 0, 0) fractional: the unitary little co-group is D3d
    // (12 operations) and the 8-point mesh reduces to 4 pairs with orbit
    // sizes 3, 3, 1, 1 (independently derived by brute force)
    const std::vector<ModuleBase::Matrix3> gm = fcc_gmatrices();
    const ModuleBase::Vector3<double> q(0.5, 0.0, 0.0);
    const std::vector<ModuleDFPT::KQPair> pairs =
        ModuleDFPT::kq_pair_reps(fcc_latvec(), gm.data(), gm.size(), q, mesh_2x2x2());
    ASSERT_EQ(pairs.size(), 4u);
    const std::vector<int> sizes = orbit_sizes(pairs);
    EXPECT_EQ(sizes[0], 3);
    EXPECT_EQ(sizes[1], 3);
    EXPECT_EQ(sizes[2], 1);
    EXPECT_EQ(sizes[3], 1);
    // every pair carries the k + q arithmetic
    for (const ModuleDFPT::KQPair &p : pairs)
    {
        ModuleBase::Vector3<double> expect = p.k_frac + q;
        for (int c = 0; c < 3; ++c)
        {
            double v = (c == 0) ? expect.x : ((c == 1) ? expect.y : expect.z);
            v -= std::floor(v + 1.0e-9);
            if (std::abs(v - 1.0) < 1.0e-9)
            {
                v = 0.0;
            }
            const double got = (c == 0) ? p.kq_frac.x : ((c == 1) ? p.kq_frac.y : p.kq_frac.z);
            EXPECT_NEAR(got, v, 1.0e-9);
        }
    }
}

TEST(DFPTKqPairs, XPointReducesUnderTheD4hUnitaryGroup)
{
    // q = X = (1/2, 1/2, 0) fractional: the unitary little co-group is
    // D4h (16 operations) and the mesh reduces to 4 pairs with orbit
    // sizes 4, 2, 1, 1
    const std::vector<ModuleBase::Matrix3> gm = fcc_gmatrices();
    const ModuleBase::Vector3<double> q(0.5, 0.5, 0.0);
    const std::vector<ModuleDFPT::KQPair> pairs =
        ModuleDFPT::kq_pair_reps(fcc_latvec(), gm.data(), gm.size(), q, mesh_2x2x2());
    ASSERT_EQ(pairs.size(), 4u);
    const std::vector<int> sizes = orbit_sizes(pairs);
    EXPECT_EQ(sizes[0], 4);
    EXPECT_EQ(sizes[1], 2);
    EXPECT_EQ(sizes[2], 1);
    EXPECT_EQ(sizes[3], 1);
}

TEST(DFPTKqPairs, GenericQKeepsTheFullMesh)
{
    // a generic q has the trivial unitary little co-group: every mesh
    // point is its own representative with weight 1
    const std::vector<ModuleBase::Matrix3> gm = fcc_gmatrices();
    const ModuleBase::Vector3<double> q(0.123, 0.045, 0.311);
    const std::vector<ModuleDFPT::KQPair> pairs =
        ModuleDFPT::kq_pair_reps(fcc_latvec(), gm.data(), gm.size(), q, mesh_2x2x2());
    ASSERT_EQ(pairs.size(), mesh_2x2x2().size());
    for (const ModuleDFPT::KQPair &p : pairs)
    {
        EXPECT_EQ(p.orbit_size, 1);
    }
}

TEST(DFPTKqPairs, EmptyInputsAreEmptyOutputs)
{
    const std::vector<ModuleBase::Matrix3> gm = fcc_gmatrices();
    const ModuleBase::Vector3<double> q(0.0, 0.0, 0.0);
    EXPECT_TRUE(ModuleDFPT::kq_pair_reps(fcc_latvec(), nullptr, 0, q, mesh_2x2x2()).empty());
    EXPECT_TRUE(ModuleDFPT::kq_pair_reps(fcc_latvec(), gm.data(), gm.size(), q, {}).empty());
}
