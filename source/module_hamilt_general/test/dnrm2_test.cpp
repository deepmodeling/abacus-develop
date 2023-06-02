#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../module_ewald/dnrm2.h"

/************************************************
 *  unit test of dnrm2.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - dnrm2(const int n, const double *x, const int incx): 
 *      - compute the Euclidean length (12 norm) of std::vector x, 
 *         with scaling of input to avoid destructive underflow and overflow
 */

class EwaldTest : public ::testing::Test
{
    
};

TEST_F(EwaldTest,Dnrm2Test)
{
    std::string output;
    double x[3]={1.5,2.5,3.5} ;
    int incx1=1;
    int n1=-1;
    int n2=0;
    int n3=3;
    // case 1
	testing :: internal :: CaptureStderr();
    EXPECT_EQ(dnrm2(n1,x,incx1), 0);
	output = testing::internal::GetCapturedStderr();
	EXPECT_THAT(output,testing::HasSubstr("error in dnrm2, n < 0 or incx <= 0,"));
	testing :: internal :: CaptureStderr();
    EXPECT_EQ(dnrm2(n3,x,0), 0);
	output = testing::internal::GetCapturedStderr();
	EXPECT_THAT(output,testing::HasSubstr("error in dnrm2, n < 0 or incx <= 0,"));
    // case 2
    EXPECT_EQ(dnrm2(n2,x,incx1), 0);
    // case 3
    EXPECT_EQ(dnrm2(n3,x,incx1), sqrt(1.5*1.5+2.5*2.5+3.5*3.5));
}

