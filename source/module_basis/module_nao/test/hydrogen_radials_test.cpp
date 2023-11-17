#include <gtest/gtest.h>
#include <iostream>
#include "module_basis/module_nao/hydrogen_radials.h"

#ifdef __MPI
#include <mpi.h>
#endif

/// @brief Simpson integral
/// @param x variable stored in a vector
/// @param y function value stored in a vector
/// @return integral value
/// @note used to check the normalization of the radial functions
double SimpsonIntegral(const std::vector<double> &x, const std::vector<double> &y)
{
    double result = 0.0;
    result += (y[0] + y[y.size() - 1]);
    /* x and y must have the same size and their length must be a same odd number */
    assert(x.size() == y.size());
    assert(x.size() % 2 == 1);
    double h = x[1] - x[0];
    for (int i = 1; i < x.size() - 1; i+=2)
    {
        result += (2*y[i] + 4*y[i+1]);
    }
    result *= h/3;
    return result;
}

class HydrogenRadialsTest : public ::testing::Test
{
    protected:
        virtual void SetUp()
        {
            // set up the test case
            itype_ = 1;
            charge_ = 1.0;
            nmax_ = 1;
            rcut_ = 20.0;
            dr_ = 0.01;
            rank_ = 0;
            ptr_log_ = NULL;
        }

        virtual void TearDown()
        {
            // tear down the test case
        }

        int itype_;
        double charge_;
        int nmax_;
        double rcut_;
        double dr_;
        int rank_;
        std::ofstream* ptr_log_;
};

TEST_F(HydrogenRadialsTest, BuildHydrogenTest)
{
    HydrogenRadials hr;
    hr.build(itype_, charge_, nmax_, rcut_, dr_, rank_, ptr_log_);
    int lmax = nmax_ - 1;
    // check maximal angular momentum
    EXPECT_EQ(hr.lmax(), lmax);
    // check number of zeta for each l
    for(int l=0; l<=lmax; ++l)
    {
        EXPECT_EQ(hr.nzeta(l), nmax_ - l);
    }

    int nrgrid = static_cast<int>(rcut_/dr_) + 1;
    std::vector<double> r;

    for(int ir = 0; ir < nrgrid; ++ir)
    {
        r.push_back(ir * dr_);
    }
    for(int n=1; n<=nmax_; n++)
    {
        for(int l=0; l<n; l++)
        {
            std::vector<double> rvalue;
            for(int ir=0; ir<nrgrid; ++ir)
            {
                rvalue.push_back(std::pow(hr.chi(l, n-l-1).rvalue(ir), 2) * r[ir] * r[ir]);
            }
            double norm = SimpsonIntegral(r, rvalue);
            EXPECT_NEAR(norm, 1.0, 1e-3);
            rvalue.clear();
            rvalue.shrink_to_fit();
        }
    }
}

TEST_F(HydrogenRadialsTest, BuildSiliconTest)
{
    // Si: 1s2 2s2 2p6 3s2 3p2, but also generate 3d
    charge_ = 14.0;
    nmax_ = 3;
    HydrogenRadials hr;
    hr.build(itype_, charge_, nmax_, rcut_, dr_, rank_, ptr_log_);
    int lmax = nmax_ - 1;
    // check maximal angular momentum
    EXPECT_EQ(hr.lmax(), lmax);
    // check number of zeta for each l
    for(int l=0; l<=lmax; ++l)
    {
        EXPECT_EQ(hr.nzeta(l), nmax_ - l);
    }

    int nrgrid = static_cast<int>(rcut_/dr_) + 1;
    std::vector<double> r;

    for(int ir = 0; ir < nrgrid; ++ir)
    {
        r.push_back(ir * dr_);
    }
    for(int n=1; n<=nmax_; n++)
    {
        for(int l=0; l<n; l++)
        {
            std::vector<double> rvalue;
            for(int ir=0; ir<nrgrid; ++ir)
            {
                rvalue.push_back(std::pow(hr.chi(l, n-l-1).rvalue(ir), 2) * r[ir] * r[ir]);
            }
            double norm = SimpsonIntegral(r, rvalue);
            EXPECT_NEAR(norm, 1.0, 1e-3);
            rvalue.clear();
            rvalue.shrink_to_fit();
        }
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    return result;
}