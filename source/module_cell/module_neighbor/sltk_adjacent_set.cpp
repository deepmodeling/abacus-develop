#include <stdexcept>
#include <string>
#include "sltk_adjacent_set.h"


long AdjacentSet::call_times = 0;
int AdjacentSet::dx = 3;
int AdjacentSet::dy = 3;
int AdjacentSet::dz = 3;
int AdjacentSet::center = 13;
bool AdjacentSet::expand_flag = false;
int AdjacentSet::trueX = 1;
int AdjacentSet::trueY = 1;
int AdjacentSet::trueZ = 1;


AdjacentSet::AdjacentSet()
{
	this->length = 0;
}

AdjacentSet::~AdjacentSet() {}



void AdjacentSet::set
(
    const int box_x,
    const int box_y,
    const int box_z,
	const int offset_in,
	const int &test_grid_in
)
{
	int index_box=0; // mohan update 2021-06-22

	if (!expand_flag)
	{
		index_box = AdjacentSet::index(box_x, box_y, box_z);
	}
	else
	{
		index_box = AdjacentSet::index_expand(box_x, box_y, box_z);
	}

	if (test_grid_in > 3)
	{
/*
		GlobalV::ofs_running << "\n box(" << box_x
		<< " " << box_y
		<< " " << box_z
		<< " ) offset=" << offset_in
		<< " index = " << index_box;
*/
	}

	this->box.push_back(index_box);

	this->offset.push_back(offset_in);
	this->length++;


	//std::cout << "\n length = " << length << std::endl;
	call_times++;
	return;
}

void AdjacentSet::assertCoordinateIsRight(const int x, const int y, const int z, const char* const s)
{
	/*	if (x > 1 || x < -1 || y > 1 || y < -1 || z > 1 || z < -1)
		{
			const std::string error(static_cast<std::string>("Function ") + s + ": exception std::out_of_range\n"
			                        "\tThe coordinate of grid must be in {-1 0 1}.");
			throw std::out_of_range(error);
		}*/
}

void AdjacentSet::getBox
(
    const int value,
    int &x, // = box_x
    int &y, // = box_y
    int &z  // = box_z
)
{
	if (!expand_flag)
	{
		x = value / 9 - 1;
		int now = value - 9 * (x + 1);
		y = now / 3 - 1;
		now = now -  3 * (y + 1);
		z = now - 1;

//		if (GlobalV::test_grid > 3) GlobalV::ofs_running << "\n value=" << value << " x=" << x << " y=" << y << " z=" << z;
	}
	else
	{
		// mohan remain question??? why fx, fy, fz should be double, not int?
		double fx = static_cast<double>(value / (dy * dz));
		int now = value - dy * dz * (int)fx;
		double fy = static_cast<double>(now / dz);
		double fz = now - fy * dz;
		x = (int)fx - trueX;
		y = (int)fy - trueY;
		z = (int)fz - trueZ;

//		std::cout << "\n value=" << value << " x=" << x << " y=" << y << " z=" << z;
//		std::cout << "\n fx=" << fx << " fy=" << y << " fz=" << z;
//		std::cout << "\n truex=" << trueX << " truey=" << trueY << " truez=" << trueX;		
//		BLOCK_HERE("get box");
	}

	return;
}
//2015-05-07
void AdjacentSet::delete_vector(void)
{
	std::vector<int>().swap(box);
	std::vector<int>().swap(offset);
}
