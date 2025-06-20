#include "../mathzone.h"
#include "../matrix3.h"
#include "../vector3.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

/************************************************
 *  unit test of class Mathzone
 ***********************************************/

/**
 * - Tested Functions:
 *   - PointwiseProduct
 *     - return a vector, which is the pointwise
 *     - product of another two vectors of the same
 *     - length
 *   - Direct2Cartesian
 *     - change atomic coordinates from direct
 *     - to Cartesian
 *   - Cartesian2Direct
 *     - change atomic coordinates from Cartesian
 *     - to Direct
 */

class MathzoneTest : public testing::Test
{
protected:
	double R11 = 3.68; double R12 = 0.00; double R13 = 0.00;
	double R21 = 0.00; double R22 = 10.1; double R23 = 0.00;
	double R31 = 0.00; double R32 = 0.00; double R33 = 26.7;
	ModuleBase::Matrix3 lattice;
	ModuleBase::Vector3<double> direct, cartesian;
};

TEST_F(MathzoneTest, PointwiseProduct)
{
	std::vector<double> aa, bb, cc;
	for(int i=0;i<10;i++)
	{
		aa.push_back(i*i);
		bb.push_back(i*2);
	}
	cc = ModuleBase::Mathzone::Pointwise_Product(aa,bb);
	for(int i=0;i<10;i++)
	{
		EXPECT_EQ(cc[i],i*i*i*2);
	}
}

TEST_F(MathzoneTest, Direct2Cartesian)
{
	direct.set(0.1,0.2,0.4);
	cartesian.set(0.368,2.02,10.68);
	ModuleBase::Vector3<double> cartnew;
	ModuleBase::Mathzone::Direct_to_Cartesian(direct.x,
			direct.y,
			direct.z,
			R11, R12, R13,
			R21, R22, R23,
			R31, R32, R33,
			cartnew.x,
			cartnew.y,
			cartnew.z);
	EXPECT_NEAR(cartnew.x,cartesian.x, 1e-15);
	EXPECT_NEAR(cartnew.y,cartesian.y, 1e-15);
	EXPECT_NEAR(cartnew.z,cartesian.z, 1e-15);
}

TEST_F(MathzoneTest, Cartesian2Direct)
{
	direct.set(0.1,0.2,0.4);
	cartesian.set(0.368,2.02,10.68);
	ModuleBase::Vector3<double> directnew;
	ModuleBase::Mathzone::Cartesian_to_Direct(cartesian.x,
			cartesian.y,
			cartesian.z,
			R11, R12, R13,
			R21, R22, R23,
			R31, R32, R33,
			directnew.x,
			directnew.y,
			directnew.z);
	EXPECT_NEAR(directnew.x,direct.x, 1e-15);
	EXPECT_NEAR(directnew.y,direct.y, 1e-15);
	EXPECT_NEAR(directnew.z,direct.z, 1e-15);
}
