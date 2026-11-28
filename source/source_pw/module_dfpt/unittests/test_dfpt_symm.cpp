#include "source_pw/module_dfpt/dfpt_symm.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"

#include "gtest/gtest.h"

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
