#include "module_basis/module_nao/real_gaunt_table.h"
#include "gtest/gtest.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include "module_basis/module_ao/ORB_gaunt_table.h"

/***********************************************************
 *      Unit test of class "RealGauntTable"
 ***********************************************************/
/*!
 *  Tested functions:
 *
 *  - build
 *      - constructs the Gaunt table of real spherical harmonics
 *
 *  - lmax
 *      - gets the maximum angular momentum
 *
 *  - operator()
 *      - gets a specific Gaunt coefficient of real spherical harmonics
 *
 *                                                          */
class RealGauntTableTest : public ::testing::Test
{
  protected:
    void SetUp() { rgt.build(lmax); }
    void TearDown() {}

    int lmax = 10;      //!< maximum angular momentum
    RealGauntTable rgt; //!< object under test
    const double tol = 1e-12; //!< numerical error tolerance for individual Gaunt coefficient
};

TEST_F(RealGauntTableTest, LegacyConsistency)
{
    // this test checks whether the coefficients in RealGauntTable is consistent with those of ORB_gaunt_table
    // this test shall be removed in the future once the refactoring is finished
    ORB_gaunt_table ogt;
    ogt.init_Gaunt_CH(lmax);
    ogt.init_Gaunt(lmax);

    for (int l1 = 0; l1 <= lmax; ++l1)
    {
        for (int mm1 = 0; mm1 <= 2*l1; ++mm1)
        {
            int index1 = ogt.get_lm_index(l1, mm1);
            for (int l2 = 0; l2 <= lmax; ++l2)
            {
                for (int mm2 = 0; mm2 <= 2*l2; ++mm2)
                {
                    int index2 = ogt.get_lm_index(l2, mm2);
                    for (int l3 = 0; l3 <= 2*lmax; ++l3)
                    //for (int l3 = std::abs(l1-l2); l3 <= 2*lmax; l3 += 2)
                    {
                        for (int mm3 = 0; mm3 <= 2*l3; ++mm3)
                        {
                            int index3 = ogt.get_lm_index(l3, mm3);

                            int m1 = ogt.Index_M(mm1);
                            int m2 = ogt.Index_M(mm2);
                            int m3 = ogt.Index_M(mm3);

                            EXPECT_NEAR(rgt(l1, l2, l3, m1, m2, m3), ogt.Gaunt_Coefficients(index1, index2, index3), tol);
                        }
                    }
                }
            }
        }
    }
}

TEST_F(RealGauntTableTest, SanityCheck)
{
    EXPECT_EQ(rgt.lmax(), lmax);

    //EXPECT_NEAR(rgt(0, 0, 0, 0, 0, 0), 1.0, tol);
}

int main(int argc, char** argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
