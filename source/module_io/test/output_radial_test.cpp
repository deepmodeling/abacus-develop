#include <gtest/gtest.h>
#include "../output_radial.h"
#include <cmath>

class OutputRadialTest : public ::testing::Test
{
    protected:
        virtual void SetUp()
        {
            symbol_ = "H";
            ecut_ = 100.0;
            rcut_ = 10.0;
            lmax_ = 0;
            l_nchi_ = {1};
            ngrid_ = 101;
            dr_ = 0.01;
            rgrid_ = new double[ngrid_];
            chi_ = new double[ngrid_];
            for(int i = 0; i < ngrid_; ++i)
            {
                rgrid_[i] = i * dr_;
                chi_[i] = std::exp(-rgrid_[i]);
            }
        }

        virtual void TearDown()
        {
            delete[] rgrid_;
            delete[] chi_;
        }

        std::string symbol_;
        double ecut_;
        double rcut_;
        int lmax_;
        std::vector<int> l_nchi_;
        int ngrid_;
        double dr_;
        double* rgrid_;
        double* chi_;
};

TEST_F(OutputRadialTest, initialize)
{
    ModuleIO::OutputRadial output_radial;
    output_radial.initialize("test_output_radial.dat");
    EXPECT_EQ(output_radial.symbol(), "H");
    EXPECT_EQ(output_radial.ecut(), 100.0);
    EXPECT_EQ(output_radial.rcut(), 10.0);
    EXPECT_EQ(output_radial.lmax(), 0);
    EXPECT_EQ(output_radial.ngrid(), 0);
    EXPECT_EQ(output_radial.dr(), 0.01);
    EXPECT_EQ(output_radial.current_l(), 0);
    EXPECT_EQ(output_radial.current_chi(), 0);
    output_radial.finalize();
}

TEST_F(OutputRadialTest, configure)
{
    ModuleIO::OutputRadial output_radial;
    output_radial.initialize("test_output_radial.dat");
    output_radial.configure(symbol_, ecut_, rcut_, lmax_, l_nchi_.data(), ngrid_, dr_);
    EXPECT_EQ(output_radial.symbol(), symbol_);
    EXPECT_EQ(output_radial.ecut(), ecut_);
    EXPECT_EQ(output_radial.rcut(), rcut_);
    EXPECT_EQ(output_radial.lmax(), lmax_);
    for(int i = 0; i < lmax_ + 1; ++i)
    {
        EXPECT_EQ(output_radial.l_nchi()[i], l_nchi_[i]);
    }
    EXPECT_EQ(output_radial.ngrid(), ngrid_);
    EXPECT_EQ(output_radial.dr(), dr_);
    EXPECT_EQ(output_radial.current_l(), 0);
    EXPECT_EQ(output_radial.current_chi(), 0);
    output_radial.finalize();
}

TEST_F(OutputRadialTest, push)
{
    ModuleIO::OutputRadial output_radial;
    output_radial.initialize("test_output_radial.dat");
    output_radial.configure(symbol_, ecut_, rcut_, lmax_, l_nchi_.data(), ngrid_, dr_);
    output_radial.push(ngrid_, rgrid_, chi_);
    EXPECT_EQ(output_radial.current_l(), 1);
    EXPECT_EQ(output_radial.current_chi(), 0);
    output_radial.finalize();
}
