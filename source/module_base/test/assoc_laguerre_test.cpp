#include <gtest/gtest.h>
#include <vector>
#include "module_base/assoc_laguerre.h"

/*
    Unittest for associated Laguerre polynomials
*/

class AssocLaguerreTest : public ::testing::Test
{
    protected:
        AssocLaguerreTest()
        {
        }

        ~AssocLaguerreTest()
        {
        }

        void SetUp()
        {
        }

        void TearDown()
        {
        }
};

TEST_F(AssocLaguerreTest, LaguerreTest)
{
    // reference value from scipy.special.assoc_laguerre with k = 0
    Assoc_Laguerre al;
    // 0-th order Laguerre polynomial
    int n = 0;
    std::vector<double> xs = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> ref_ys;
    for(int i = 0; i < xs.size(); i++)
    {
        ref_ys.push_back(1.0);
    }
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.laguerre(n, xs[i]), ref_ys[i], 1e-6);
    }
    // 1-st order Laguerre polynomial
    n = 1;
    ref_ys = {1.0, 0.0, -1.0, -2.0};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.laguerre(n, xs[i]), ref_ys[i], 1e-6);
    }
    // 2-nd order Laguerre polynomial
    n = 2;
    ref_ys = {1.0, -0.5, -1.0, -0.5};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.laguerre(n, xs[i]), ref_ys[i], 1e-6);
    }
    // 3-rd order Laguerre polynomial
    n = 3;
    ref_ys = {1.0, -0.6666666666666666, -0.33333333333333337, 1.0};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.laguerre(n, xs[i]), ref_ys[i], 1e-6);
    }
    // 4-th order Laguerre polynomial
    n = 4;
    ref_ys = {1.0, -0.625, 0.33333333333333337, 1.375};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.laguerre(n, xs[i]), ref_ys[i], 1e-6);
    }
    // 5-th order Laguerre polynomial
    n = 5;
    ref_ys = {1.0, -0.4666666666666667, 0.7333333333333334, 0.8500000000000001};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.laguerre(n, xs[i]), ref_ys[i], 1e-6);
    }
    // 6-th order Laguerre polynomial
    n = 6;
    ref_ys = {1.0, -0.2569444444444444, 0.8222222222222224, -0.012499999999999956};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.laguerre(n, xs[i]), ref_ys[i], 1e-6);
    }
    // 7-th order Laguerre polynomial
    n = 7;
    ref_ys = {1.0, -0.04047619047619044, 0.6634920634920637, -0.7464285714285714};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.laguerre(n, xs[i]), ref_ys[i], 1e-6);
    }
    // 8-th order Laguerre polynomial
    n = 8;
    ref_ys = {1.0, 0.1539930555555556, 0.3587301587301589, -1.1087053571428571};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.laguerre(n, xs[i]), ref_ys[i], 1e-6);
    }
    // all physical and possible n values are tested
}

TEST_F(AssocLaguerreTest, AssociateLaguerreTest)
{
    // reference value from scipy.special.assoc_laguerre
    Assoc_Laguerre al;
    // test n = 8, k = 0
    int n = 8;
    int k = 0;
    std::vector<double> xs = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> ref_ys = {1.0, 0.1539930555555556, 0.3587301587301589, -1.1087053571428571};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.associate_laguerre(n, xs[i], k), ref_ys[i], 1e-6);
    }
    // test n = 8, k = 1
    n = 8;
    k = 1;
    xs = {0.0, 1.0, 2.0, 3.0};
    ref_ys = {9.0, -1.4017609126984123, 1.5777777777777775, -0.14263392857142831};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.associate_laguerre(n, xs[i], k), ref_ys[i], 1e-6);
    }
    // test n = 8, k = 2
    n = 8;
    k = 2;
    xs = {0.0, 1.0, 2.0, 3.0};
    ref_ys = {45.0, -4.189459325396824, 0.7523809523809523, 3.1359374999999994};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.associate_laguerre(n, xs[i], k), ref_ys[i], 1e-6);
    }
    // test n = 8, k = 3
    n = 8;
    k = 3;
    xs = {0.0, 1.0, 2.0, 3.0};
    ref_ys = {165.0, 1.449231150793654, -6.384126984126984, 5.70200892857143};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.associate_laguerre(n, xs[i], k), ref_ys[i], 1e-6);
    }
    // test n = 8, k = 4
    n = 8;
    k = 4;
    xs = {0.0, 1.0, 2.0, 3.0};
    ref_ys = {495.0, 51.8809771825397, -20.098412698412695, -0.09441964285714727};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.associate_laguerre(n, xs[i], k), ref_ys[i], 1e-6);
    }
    // test n = 8, k = 5
    n = 8;
    k = 5;
    xs = {0.0, 1.0, 2.0, 3.0};
    ref_ys = {1287.0, 242.01411210317463, -21.990476190476176, -23.028348214285717};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.associate_laguerre(n, xs[i], k), ref_ys[i], 1e-6);
    }
    // test n = 8, k = 6
    n = 8;
    k = 6;
    xs = {0.0, 1.0, 2.0, 3.0};
    ref_ys = {3003.0, 777.6319692460318, 53.67301587301589, -60.999776785714275};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.associate_laguerre(n, xs[i], k), ref_ys[i], 1e-6);
    }
    // test n = 8, k = 7
    n = 8;
    k = 7;
    xs = {0.0, 1.0, 2.0, 3.0};
    ref_ys = {6435.0, 2055.2262152777776, 368.62539682539676, -75.53370535714284};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.associate_laguerre(n, xs[i], k), ref_ys[i], 1e-6);
    }
    // test n = 8, k = 8
    n = 8;
    k = 8;
    xs = {0.0, 1.0, 2.0, 3.0};
    ref_ys = {12870.0, 4777.330183531744, 1256.2666666666667, 53.21986607142856};
    for(int i = 0; i < xs.size(); i++)
    {
        EXPECT_NEAR(al.associate_laguerre(n, xs[i], k), ref_ys[i], 1e-6);
    }
    // all physical and possible n, k values are tested
}

TEST_F(AssocLaguerreTest, FactorialTest)
{
    // this test is simple, no need to compare with reference value
    Assoc_Laguerre al;
    EXPECT_DOUBLE_EQ(al.factorial(0), 1);
    EXPECT_DOUBLE_EQ(al.factorial(1), 1);
    EXPECT_DOUBLE_EQ(al.factorial(2), 2);
    EXPECT_DOUBLE_EQ(al.factorial(3), 6);
    EXPECT_DOUBLE_EQ(al.factorial(4), 24);
    EXPECT_DOUBLE_EQ(al.factorial(5), 120);
    EXPECT_DOUBLE_EQ(al.factorial(6), 720);
    EXPECT_DOUBLE_EQ(al.factorial(7), 5040);
    EXPECT_DOUBLE_EQ(al.factorial(8), 40320);
    EXPECT_DOUBLE_EQ(al.factorial(9), 362880);
    EXPECT_DOUBLE_EQ(al.factorial(10), 3628800);
}

TEST_F(AssocLaguerreTest, ValueTest)
{
    Assoc_Laguerre al;
    //test n = 1, 2, ..., 4, from 1s to 4f
    std::vector<double> xs = {0.0, 1.0, 2.0, 3.0};
    int nmax = 4;
    for(int n = 1; n < nmax; n++)
    {
        for(int l = 0; l < n; l++)
        {
            std::vector<double> ref_ys;
            for(int i = 0; i < xs.size(); i++)
            {
                ref_ys.push_back(al.associate_laguerre(n-l-1, xs[i], 2*l+1));
            }
            for(int i = 0; i < xs.size(); i++)
            {
                EXPECT_DOUBLE_EQ(al.value(n, l, xs[i]), ref_ys[i]);
            }
        }
    }
}

TEST_F(AssocLaguerreTest, GenerateVectorTest)
{
    Assoc_Laguerre al;
    //test n = 1, 2, ..., 4, from 1s to 4f
    std::vector<double> xs = {0.0, 1.0, 2.0, 3.0};
    int nmax = 4;
    for(int n = 1; n <= nmax; n++)
    {
        for(int l = 0; l < n; l++)
        {
            std::vector<double> ref_ys;
            for(int i = 0; i < xs.size(); i++)
            {
                ref_ys.push_back(al.associate_laguerre(n-l-1, xs[i], 2*l+1));
            }
            std::vector<double> ys;
            ys.resize(xs.size());
            al.generate(n, l, xs, ys);
            for(int i = 0; i < xs.size(); i++)
            {
                EXPECT_NEAR(ys[i], ref_ys[i], 1e-6);
            }
        }
    }
}

TEST_F(AssocLaguerreTest, GeneratePointerTest)
{
    Assoc_Laguerre al;
    //test n = 1, 2, ..., 4, from 1s to 4f
    std::vector<double> xs = {0.0, 1.0, 2.0, 3.0};
    int nmax = 4;
    for(int n = 1; n <= nmax; n++)
    {
        for(int l = 0; l < n; l++)
        {
            std::vector<double> ref_ys;
            for(int i = 0; i < xs.size(); i++)
            {
                ref_ys.push_back(al.associate_laguerre(n-l-1, xs[i], 2*l+1));
            }
            double* ys = new double[xs.size()];
            al.generate(n, l, xs.size(), xs.data(), ys);
            for(int i = 0; i < xs.size(); i++)
            {
                EXPECT_NEAR(ys[i], ref_ys[i], 1e-6);
            }
            delete[] ys;
        }
    }
}