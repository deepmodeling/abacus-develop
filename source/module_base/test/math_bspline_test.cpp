#include "../math_bspline.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of class Bspline
 ***********************************************/

/**
 * - Tested Functions:
 *   - Init
 *     - norder must be even
 *     - norder mush be positive
 *   - Properties
 *     - \sum_i M_n(u+i) = 1 (i=0,1,2,...n)
 *
 */

class MathBsplineTest : public testing::Test
{
protected:
	ModuleBase::Bspline bp;
	int norder;
};

TEST_F(MathBsplineTest,Init)
{
	EXPECT_DEATH(
		{
			norder = 3; // norder must be even
			bp.init(norder,0.05,0);
		},""
	);
	EXPECT_DEATH(
		{
			norder = 0; // norder must be positive
			bp.init(norder,0.05,0);
		},""
	);
}

// summation over n is unity
TEST_F(MathBsplineTest,Properties)
{
	int by = 2;
	for (norder=2;norder<=20;norder=norder+by)
	{
		bp.init(norder,1.0,0);
		bp.getbspline(0.2);
		double sum=0.0;
		//std::cout << "\n" << "norder : "<< norder<<std::endl;
		for (int i=0;i<=norder;i++)
		{
			//std::cout << i << " " << bp.bezier_ele(i) << std::endl;
			sum += bp.bezier_ele(i);
		}
		//std::cout<<"sum "<< sum<< std::endl;
		EXPECT_NEAR(sum,1.0,1.e-15);
	}
}
