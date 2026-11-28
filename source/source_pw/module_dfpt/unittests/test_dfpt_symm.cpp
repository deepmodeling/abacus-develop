#include "source_pw/module_dfpt/dfpt_symm.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
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

namespace
{

/// hexagonal primitive lattice with lat != lat^T (the basis-asymmetry
/// discriminator the fcc tests cannot exercise)
ModuleBase::Matrix3 hex_latvec()
{
    const double c_over_a = 1.6;
    return ModuleBase::Matrix3(1.0, 0.0, 0.0, -0.5, 0.5 * std::sqrt(3.0), 0.0, 0.0, 0.0, c_over_a);
}

/// the 12 cartesian rotations of the C6v point group (6-fold about z and
/// the vertical mirrors composed with the horizontal flip)
std::vector<ModuleBase::Matrix3> hex_cartesian_ops()
{
    std::vector<ModuleBase::Matrix3> ops;
    for (int n = 0; n < 6; ++n)
    {
        const double th = n * M_PI / 3.0;
        const double c = std::cos(th);
        const double s = std::sin(th);
        const ModuleBase::Matrix3 rz(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0);
        ops.push_back(rz);
        ops.push_back(
            ModuleBase::Matrix3(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0) *
            ModuleBase::Matrix3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0));
    }
    return ops;
}

/// the hexagonal operations in the fractional basis, in EITHER the test
/// convention lat R lat^-1 (dir = +1) or the production convention
/// lat R^T lat^-1 (dir = -1): the two differ for non-symmetric bases
std::vector<ModuleBase::Matrix3> hex_gmatrices(int dir)
{
    const ModuleBase::Matrix3 lat = hex_latvec();
    const ModuleBase::Matrix3 inv_lat = lat.Inverse();
    const std::vector<ModuleBase::Matrix3> ops = hex_cartesian_ops();
    std::vector<ModuleBase::Matrix3> gm;
    for (const ModuleBase::Matrix3 &r : ops)
    {
        if (dir > 0)
        {
            gm.push_back(lat * r * inv_lat);
        }
        else
        {
            gm.push_back(lat * r.Transpose() * inv_lat);
        }
    }
    return gm;
}

/// the 3x3x2 Gamma-centered mesh in fractional coordinates: its orbits
/// under C6v wrap by reciprocal-lattice vectors (e.g. (1/3, 1/3) <->
/// (2/3, 2/3)), which is what the mod-lattice membership must catch
std::vector<ModuleBase::Vector3<double>> mesh_3x3x2()
{
    std::vector<ModuleBase::Vector3<double>> mesh;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                mesh.push_back(ModuleBase::Vector3<double>(i / 3.0, j / 3.0, k / 2.0));
            }
        }
    }
    return mesh;
}

/// multiset key of a pair list for convention-robustness comparison
std::vector<std::string> pair_keys(const std::vector<ModuleDFPT::KQPair> &pairs)
{
    std::vector<std::string> keys;
    for (const ModuleDFPT::KQPair &p : pairs)
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6) << p.k_frac.x << "," << p.k_frac.y << "," << p.k_frac.z << ","
            << p.orbit_size;
        keys.push_back(oss.str());
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

} // namespace

TEST(DFPTKqPairs, HexGammaReducesUnderC6vWithWrappingOrbits)
{
    // q = 0, full C6v: the 18-point mesh reduces to 6 pairs with orbit
    // sizes 6, 6, 2, 2, 1, 1 (brute-force derived); the 2-member orbits
    // of the interior points only close through the reciprocal-lattice
    // wrap, so this test locks the mod-lattice membership
    const std::vector<ModuleBase::Matrix3> gm = hex_gmatrices(1);
    const ModuleBase::Vector3<double> q(0.0, 0.0, 0.0);
    const std::vector<ModuleDFPT::KQPair> pairs =
        ModuleDFPT::kq_pair_reps(hex_latvec(), gm.data(), gm.size(), q, mesh_3x3x2());
    ASSERT_EQ(pairs.size(), 6u);
    const std::vector<int> sizes = orbit_sizes(pairs);
    EXPECT_EQ(sizes[0], 6);
    EXPECT_EQ(sizes[1], 6);
    EXPECT_EQ(sizes[2], 2);
    EXPECT_EQ(sizes[3], 2);
    EXPECT_EQ(sizes[4], 1);
    EXPECT_EQ(sizes[5], 1);
}

TEST(DFPTKqPairs, HexKPointReducesUnderTheLittleGroup)
{
    // q = K = (1/3, 1/3, 0): the unitary little co-group has 6 members
    // (C3 plus the mirrors fixing K) and the mesh reduces to 10 pairs
    // with orbit sizes 3, 3, 3, 3, then six singletons
    const std::vector<ModuleBase::Matrix3> gm = hex_gmatrices(1);
    const ModuleBase::Vector3<double> q(1.0 / 3.0, 1.0 / 3.0, 0.0);
    const std::vector<ModuleDFPT::KQPair> pairs =
        ModuleDFPT::kq_pair_reps(hex_latvec(), gm.data(), gm.size(), q, mesh_3x3x2());
    ASSERT_EQ(pairs.size(), 10u);
    const std::vector<int> sizes = orbit_sizes(pairs);
    EXPECT_EQ(sizes[0], 3);
    EXPECT_EQ(sizes[1], 3);
    EXPECT_EQ(sizes[2], 3);
    EXPECT_EQ(sizes[3], 3);
    for (int i = 4; i < 10; ++i)
    {
        EXPECT_EQ(sizes[i], 1);
    }
}

TEST(DFPTKqPairs, HexMPointReducesUnderTheLittleGroup)
{
    // q = M = (1/2, 0, 0): the unitary little co-group has 4 members
    // and the mesh reduces to 10 pairs with orbit sizes 2 x 8, 1, 1
    const std::vector<ModuleBase::Matrix3> gm = hex_gmatrices(1);
    const ModuleBase::Vector3<double> q(0.5, 0.0, 0.0);
    const std::vector<ModuleDFPT::KQPair> pairs =
        ModuleDFPT::kq_pair_reps(hex_latvec(), gm.data(), gm.size(), q, mesh_3x3x2());
    ASSERT_EQ(pairs.size(), 10u);
    const std::vector<int> sizes = orbit_sizes(pairs);
    EXPECT_EQ(sizes[0], 2);
    for (int i = 0; i < 8; ++i)
    {
        EXPECT_EQ(sizes[i], 2);
    }
    EXPECT_EQ(sizes[8], 1);
    EXPECT_EQ(sizes[9], 1);
}

TEST(DFPTKqPairs, ProductionConventionInputMatchesTestConvention)
{
    // the production gmatrix form lat R^T lat^-1 and the test form
    // lat R lat^-1 differ by R <-> R^T, which is only the same matrix
    // set for symmetric bases (fcc); on the hexagonal lattice they are
    // genuinely different matrices, and the reconstruction must map
    // both onto valid rotations so the pair lists agree exactly
    const ModuleBase::Vector3<double> q(0.0, 0.0, 0.0);
    const std::vector<ModuleBase::Matrix3> gm_test = hex_gmatrices(1);
    const std::vector<ModuleBase::Matrix3> gm_prod = hex_gmatrices(-1);
    const std::vector<ModuleDFPT::KQPair> pairs_test =
        ModuleDFPT::kq_pair_reps(hex_latvec(), gm_test.data(), gm_test.size(), q, mesh_3x3x2());
    const std::vector<ModuleDFPT::KQPair> pairs_prod =
        ModuleDFPT::kq_pair_reps(hex_latvec(), gm_prod.data(), gm_prod.size(), q, mesh_3x3x2());
    EXPECT_EQ(pair_keys(pairs_test), pair_keys(pairs_prod));
}

TEST(DFPTLittleGroup, HexagonalLittleGroups)
{
    // Gamma keeps all 12 C6v operations; M = (1/2, 0, 0) is
    // time-reversal invariant so its union equals the 4-member unitary
    // co-group; K = (1/3, 1/3, 0) is not TR invariant and the
    // antiunitary branch doubles its 6 unitary members to 12
    // the little-group function takes CARTESIAN q (the fcc tests pass
    // (pi/a)-scale vectors); the hexagonal high-symmetry points are
    // converted through the reciprocal map 2pi latvec^-1
    const std::vector<ModuleBase::Matrix3> gm = hex_gmatrices(1);
    const ModuleBase::Matrix3 recip = hex_latvec().Inverse() * ModuleBase::TWO_PI;
    const ModuleBase::Vector3<double> qg(0.0, 0.0, 0.0);
    EXPECT_EQ(ModuleDFPT::little_group_of_q(hex_latvec(), gm.data(), gm.size(), qg).size(), 12u);
    const ModuleBase::Vector3<double> qm = recip * ModuleBase::Vector3<double>(0.5, 0.0, 0.0);
    EXPECT_EQ(ModuleDFPT::little_group_of_q(hex_latvec(), gm.data(), gm.size(), qm).size(), 4u);
    const ModuleBase::Vector3<double> qk = recip * ModuleBase::Vector3<double>(1.0 / 3.0, 1.0 / 3.0, 0.0);
    EXPECT_EQ(ModuleDFPT::little_group_of_q(hex_latvec(), gm.data(), gm.size(), qk).size(), 12u);
}
