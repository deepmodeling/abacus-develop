#include "module_base/abfs-vector3_order.h"
#include "module_base/vector3.h"
#include "gtest/gtest.h"
#include "module_ri/abfs.h"

/************************************************
 *  unit test of functions in class Abfs::vector3_order
 ***********************************************/

/**
 * - Tested Function
 *   - Two constructor functions of Vector3_Order
 *   - overloaded operator of <
 */

TEST(AbfsVector3Order,Vector3Order)
{
	ModuleBase::Vector3<double> vd31 (10.,10.,10.);
	Abfs::Vector3_Order<double> vdo31(vd31);
	Abfs::Vector3_Order<double> vdo32(10.0,10.0,10.0);
	EXPECT_FALSE(vdo31<vdo32);
	ModuleBase::Vector3<int> vi31 (2,2,2);
	Abfs::Vector3_Order<int> vio31(vi31);
	Abfs::Vector3_Order<int> vio32(2,2,2);
	EXPECT_FALSE(vio31<vio32);
}
