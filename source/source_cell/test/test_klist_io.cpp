#include "../klist_io.h"
#include <gtest/gtest.h>

TEST(KListIO, PreserveFullMeshAcrossPools)
{
    const std::vector<int> spins = {0, 0};
    const std::vector<double> weights = {1.0 / 3, 2.0 / 3};
    const std::vector<ModuleBase::Vector3<double>> irreducible = {{0, 0, 0}, {0, 0, 1.0 / 3}};
    const std::vector<ModuleBase::Vector3<double>> full = {{0, 0, 0}, {0, 0, 1.0 / 3}, {0, 0, 2.0 / 3}};
    std::vector<int> spin_buffer(2);
    std::vector<double> weight_buffer(2);
    std::vector<double> cart_buffer(6);
    std::vector<double> frac_buffer(6);
    std::vector<double> full_buffer(9, -9.0);
    KListIO::pack_kpts(spins, weights, irreducible, irreducible, full, 2,
                       spin_buffer, weight_buffer, cart_buffer, frac_buffer, full_buffer);
    EXPECT_DOUBLE_EQ(full_buffer[8], 2.0 / 3);

    std::vector<int> local_spin(1);
    std::vector<double> local_weight(1);
    std::vector<ModuleBase::Vector3<double>> local_cart(1);
    std::vector<ModuleBase::Vector3<double>> local_frac(1);
    std::vector<ModuleBase::Vector3<double>> local_full(1);
    // A second pool holds only the second SCF representative, but must retain
    // every original full-grid point in its original order.
    KListIO::unpack_kpts(spin_buffer, weight_buffer, cart_buffer, frac_buffer, full_buffer,
                         1, 1, local_spin, local_weight, local_cart, local_frac, local_full);
    EXPECT_DOUBLE_EQ(local_cart[0].z, 1.0 / 3);
    EXPECT_DOUBLE_EQ(local_weight[0], 2.0 / 3);
    ASSERT_EQ(local_full.size(), full.size());
    for (std::size_t ik = 0; ik < full.size(); ++ik)
    {
        EXPECT_DOUBLE_EQ(local_full[ik].x, full[ik].x);
        EXPECT_DOUBLE_EQ(local_full[ik].y, full[ik].y);
        EXPECT_DOUBLE_EQ(local_full[ik].z, full[ik].z);
    }
}
