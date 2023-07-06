#ifndef ADJACENTSET_H
#define ADJACENTSET_H

#include <vector>
#include <stdexcept>
#include <iostream>

class out_of_range;

//==========================================================
// CLASS : Adjacent
// provide the data structure to represent the adjacent Atom
// set of one Atom.
//==========================================================

class AdjacentSet
{
private:
//==========================================================
// MEMBER STATIC FUNCTION :
// NAME : assertCoordinateIsRight
// x, y, z is the relative coordinate of the grid. The value
// of x, y, z must be in {-1, 0, 1}.
// If the value of x, y, z is illegal, this function will
// throw the exception "out_of_range".
//==========================================================
	static
	void assertCoordinateIsRight(const int x, const int y, const int z, const char* const s);	// Peize Lin delete throw(std::out_of_range) 2018-07-14

//==========================================================
// MEMBER VARIABLE :
//==========================================================
	int length;
	static int dx;
	static int dy;
	static int dz;
	static int center;
	static bool expand_flag;
	static int trueX;
	static int trueY;
	static int trueZ;
	//==========================================================
// MEMBER STATIC FUNCTION :
// NAME : index
// NAME : index_expand
// This function transform the coordinate (x, y, z) to the
// index of array,we can also save (x,y,z),but that's a
// big memory cost when AdjacentSet is large.
// RETURN VALUE : the index
//==========================================================
	static short index(const short x, const short y, const short z)
	{ return (x*9 + y*3 + z + 13); }

	static
	short index_expand(const short x, const short y, const short z)
	{ return (x*dy*dz + y*dz + z + center); }

public:
//==========================================================
// Constructors and destructor
//==========================================================
	AdjacentSet();
	~AdjacentSet();
	//2015-05-07
	void delete_vector(void);

	int getLength(void) const {return length;}

	static long call_times;
	static void setDx(const int dx_in) { dx = dx_in; }
	static void setDy(const int dy_in) { dy = dy_in; }
	static void setDz(const int dz_in) { dz = dz_in; }
	static void setCenter(const int center_in) { center = center_in; }
	static void setExpandFlag(const int flag) { expand_flag = flag; }
	static void setTrueX(const int trueX_in) { trueX = trueX_in; }
	static void setTrueY(const int trueY_in) { trueY = trueY_in; }
	static void setTrueZ(const int trueZ_in) { trueZ = trueZ_in; }

	void set
	(
	    const int box_x,
	    const int box_y,
	    const int box_z,
		const int offset,
		const int& test_grid_in
	);

	static void getBox
	(
	    const int value,
	    int &box_x,
	    int &box_y,
	    int &box_z
	);

//==========================================================
// MEMBER VARIABLES :
// NAME : offset
// record position of adjacent atoms in the workspace
// NAME : box
// record the box the adjacent atom belongs to
//==========================================================
//	std::vector<short> offset;
//	std::vector<short> box;

	// mohan modify 2010-07-01
	// short is available about 32000.
	// but when I use Si containing 768 atoms,
	// the periodic condition will expand
	// the atoms to more than 40000.
	// so I change short to int.
	std::vector<int> offset;
	std::vector<int> box;

};

#endif
