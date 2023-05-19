#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../sltk_grid.h"

/************************************************
 *  unit test of sltk_grid
 ***********************************************/

/**
 * - Tested Functions:
 *   - Grid::get the argu ment
 *     - get the dx, dy, dz, dx_min, dy_min, dz_min
 */

AtomLink::AtomLink(const FAtom &atom, AtomLink* const pointNext)
: fatom(atom), next_p(pointNext) {}
Grid::Grid(const int &test_grid_in):test_grid(test_grid_in)
{
	init_cell_flag = false;
	this->atomlink = new AtomLink[1];
}
Grid::~Grid()
{
	delete[] atomlink;
	this->delete_Cell();
}
FAtom::FAtom()
{
	d_x = 0.0;	
	d_y = 0.0;	
	d_z = 0.0;	
	as = nullptr;
	type = 0;		
	natom = 0;
}
FAtom::~FAtom(){}

class sltkgrid : public testing::Test
{
};

TEST_F(sltkgrid,getTheArgument)
{
    Grid test(1);
    test.dx=1;
    test.dy=2;
    test.dz=3;
    test.d_minX=4;
    test.d_minY=5;
    test.d_minZ=6;
    double testx=test.getCellX();
    double testy=test.getCellY();
    double testz=test.getCellZ();
    double testminx=test.getD_minX();
    double testminy=test.getD_minY();
    double testminz=test.getD_minZ();
    EXPECT_EQ(testx,1);
    EXPECT_EQ(testy,2);
    EXPECT_EQ(testz,3);
    EXPECT_EQ(testminx,4);
    EXPECT_EQ(testminy,5);
    EXPECT_EQ(testminz,6);
}

