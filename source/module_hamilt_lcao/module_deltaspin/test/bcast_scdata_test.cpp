#include <algorithm>
#include <string>

#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#ifdef __MPI
#include "mpi.h"
#endif

/************************************************
 *  unit test of bcast_ScData
 ***********************************************/

/**
 * - Tested functions:
 *  - SpinConstrain::bcast_ScData()
 *    - bcast the ScData from root to all other ranks
 */

K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

template <typename T>
class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<T, psi::DEVICE_CPU>& sc = SpinConstrain<T, psi::DEVICE_CPU>::getScInstance();
};

using MyTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(SpinConstrainTest, MyTypes);

#ifdef __MPI
TYPED_TEST(SpinConstrainTest, BcastScData)
{
    std::map<int, int> atomCounts = {
        {0, 1},
        {1, 5}
    };
    this->sc.set_atomCounts(atomCounts);
    std::map<int, int> orbitalCounts = {
        {0, 1},
        {1, 1}
    };
    this->sc.set_orbitalCounts(orbitalCounts);
    this->sc.set_decay_grad_switch(true);
    std::string sc_file = "./support/sc_f2.json";
    this->sc.bcast_ScData(sc_file, this->sc.get_nat(), this->sc.get_ntype());
    for (int iat = 0; iat < this->sc.get_nat(); iat++)
    {
        if (iat == 1)
        {
            EXPECT_NEAR(this->sc.get_sc_lambda().data()[iat].x, 0.1 * 7.349864435130999e-05, 1e-12);
            EXPECT_NEAR(this->sc.get_sc_lambda().data()[iat].y, 0.1 * 7.349864435130999e-05, 1e-12);
            EXPECT_NEAR(this->sc.get_sc_lambda().data()[iat].z, 0.2 * 7.349864435130999e-05, 1e-12);
        }
        else if (iat == 5)
        {
            EXPECT_DOUBLE_EQ(this->sc.get_target_mag().data()[iat].x, 0.0);
            EXPECT_DOUBLE_EQ(this->sc.get_target_mag().data()[iat].y, 1.5);
            EXPECT_DOUBLE_EQ(this->sc.get_target_mag().data()[iat].z, 0.0);
        }
    }
    EXPECT_DOUBLE_EQ(this->sc.get_decay_grad().data()[0], 0.0);
    EXPECT_NEAR(this->sc.get_decay_grad().data()[1], 0.9 * 13.605698, 1e-12);
}

int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);

    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#endif