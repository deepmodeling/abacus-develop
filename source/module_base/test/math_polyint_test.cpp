#include "../math_polyint.h"

#include <math.h>

#include "../realarray.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define doublethreshold 1e-9

/************************************************
*  unit test of class PolyInt
***********************************************/

/**
 * This unit test is to verify the accuracy of 
 * interpolation method on the function sin(x)/x
 * with a interval of 0.01.
 * sin(x)/x is one of the solution of spherical bessel
 * function when l=0.
 * 
 * - Tested function:
 *      - 4 types of Polynomial_Interpolation 
 *      - Polynomial_Interpolation_xy
 */


class bessell0 : public testing::Test
{
    protected:
    
    int TableLength = 400;
    double interval = 0.01;
    ModuleBase::realArray table3,table4;
    ModuleBase::realArray y3;
    double *tablex;
    double *tabley;

    double sinc(double x) {return sin(x)/x;}

    void SetUp()
    {
        tablex = new double[TableLength];
        tabley = new double[TableLength];
        table3.create(1,1,TableLength);
        table4.create(1,1,1,TableLength);
        y3.create(1,1,TableLength);

        for(int i=1;i<TableLength;++i) 
        {
            table3(0,0,i) = sinc(i * interval);
            table4(0,0,0,i) = sinc(i * interval); 
            tablex[i] = i * interval;
            tabley[i] = sinc(i * interval);
        }
    }

    void TearDown()
    {
        delete [] tablex;
        delete [] tabley;
    }
};


TEST_F(bessell0,PolynomialInterpolationThreeDimensionY)
{
    ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,y3,1,TableLength,interval,0.1);
    ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,y3,2,TableLength,interval,1.005);
    ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,y3,3,TableLength,interval,2.005);
    ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,y3,4,TableLength,interval,3.005);
    ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,y3,5,TableLength,interval,3.505);

    EXPECT_NEAR(y3(0,0,1),sinc(0.1),doublethreshold);
    EXPECT_NEAR(y3(0,0,2),sinc(1.005),doublethreshold);
    EXPECT_NEAR(y3(0,0,3),sinc(2.005),doublethreshold);
    EXPECT_NEAR(y3(0,0,4),sinc(3.005),doublethreshold);
    EXPECT_NEAR(y3(0,0,5),sinc(3.505),doublethreshold);
}

TEST_F(bessell0,PolynomialInterpolationThreeDimension)
{
    double y1 = ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,TableLength,interval,0.1);
    double y2 = ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,TableLength,interval,1.005);
    double y3 = ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,TableLength,interval,2.005);
    double y4 = ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,TableLength,interval,3.005);
    double y5 = ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,TableLength,interval,3.505);

    EXPECT_NEAR(y1,sinc(0.1),doublethreshold);
    EXPECT_NEAR(y2,sinc(1.005),doublethreshold);
    EXPECT_NEAR(y3,sinc(2.005),doublethreshold);
    EXPECT_NEAR(y4,sinc(3.005),doublethreshold);
    EXPECT_NEAR(y5,sinc(3.505),doublethreshold);
    testing::internal::CaptureStdout();
    EXPECT_EXIT( double y6 = ModuleBase::PolyInt::Polynomial_Interpolation(table3,0,0,TableLength,interval,3.9997); ,::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Not enough space allocated for radial FFT"));
}

TEST_F(bessell0,PolynomialInterpolationFourDimension)
{
    double y1 = ModuleBase::PolyInt::Polynomial_Interpolation(table4,0,0,0,TableLength,interval,0.1);
    double y2 = ModuleBase::PolyInt::Polynomial_Interpolation(table4,0,0,0,TableLength,interval,1.005);
    double y3 = ModuleBase::PolyInt::Polynomial_Interpolation(table4,0,0,0,TableLength,interval,2.005);
    double y4 = ModuleBase::PolyInt::Polynomial_Interpolation(table4,0,0,0,TableLength,interval,3.005);
    double y5 = ModuleBase::PolyInt::Polynomial_Interpolation(table4,0,0,0,TableLength,interval,3.505);

    EXPECT_NEAR(y1,sinc(0.1),doublethreshold);
    EXPECT_NEAR(y2,sinc(1.005),doublethreshold);
    EXPECT_NEAR(y3,sinc(2.005),doublethreshold);
    EXPECT_NEAR(y4,sinc(3.005),doublethreshold);
    EXPECT_NEAR(y5,sinc(3.505),doublethreshold);
}

TEST_F(bessell0,PolynomialInterpolation)
{
    double y1 = ModuleBase::PolyInt::Polynomial_Interpolation(tabley,TableLength,interval,0.1);
    double y2 = ModuleBase::PolyInt::Polynomial_Interpolation(tabley,TableLength,interval,1.005);
    double y3 = ModuleBase::PolyInt::Polynomial_Interpolation(tabley,TableLength,interval,2.005);
    double y4 = ModuleBase::PolyInt::Polynomial_Interpolation(tabley,TableLength,interval,3.005);
    double y5 = ModuleBase::PolyInt::Polynomial_Interpolation(tabley,TableLength,interval,3.505);

    EXPECT_NEAR(y1,sinc(0.1),doublethreshold);
    EXPECT_NEAR(y2,sinc(1.005),doublethreshold);
    EXPECT_NEAR(y3,sinc(2.005),doublethreshold);
    EXPECT_NEAR(y4,sinc(3.005),doublethreshold);
    EXPECT_NEAR(y5,sinc(3.505),doublethreshold);
}

TEST_F(bessell0,PolynomialInterpolationXY)
{
    double y1 = ModuleBase::PolyInt::Polynomial_Interpolation_xy(tablex,tabley,TableLength,0.1);
    double y2 = ModuleBase::PolyInt::Polynomial_Interpolation_xy(tablex,tabley,TableLength,1.005);
    double y3 = ModuleBase::PolyInt::Polynomial_Interpolation_xy(tablex,tabley,TableLength,2.005);
    double y4 = ModuleBase::PolyInt::Polynomial_Interpolation_xy(tablex,tabley,TableLength,3.005);
    double y5 = ModuleBase::PolyInt::Polynomial_Interpolation_xy(tablex,tabley,TableLength,3.505);

    EXPECT_NEAR(y1,sinc(0.1),doublethreshold);
    EXPECT_NEAR(y2,sinc(1.005),doublethreshold);
    EXPECT_NEAR(y3,sinc(2.005),doublethreshold);
    EXPECT_NEAR(y4,sinc(3.005),doublethreshold);
    EXPECT_NEAR(y5,sinc(3.505),doublethreshold);
}
