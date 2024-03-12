#ifndef INCLUDE_FATOM
#define INCLUDE_FATOM

#include <memory>
#include "sltk_util.h"
#include "sltk_adjacent_set.h"

class AdjacentSet;

// a class contains the atom position, 
// the type and the index,
class FAtom
{
private:
	double d_x;
	double d_y;
	double d_z;
	std::shared_ptr<AdjacentSet> as;

	int type;
	int natom;

public:
//==========================================================
// Default Constructor and deconstructor
//==========================================================

	FAtom();
	~FAtom();
//2015-05-07
	void delete_vector(void);

//	static int count1;
//	static int count2;

//==========================================================
// MEMBER FUNCTION :
// NAME : setAdjacent
// Dangerous but high performance interface function!
// no exception test.
//
// NAME : getAdjacentSet
//
// NAME : setAdjacentSet
//==========================================================

	std::shared_ptr<AdjacentSet> getAdjacentSet() const
	{ return this->as; }

	void allocate_AdjacentSet(void)
	{ this->as = std::make_shared<AdjacentSet>(); }

//==========================================================
// MEMBER FUNCTION :
// EXPLAIN : get value
//==========================================================
	const double& x() const { return d_x; }
	const double& y() const { return d_y; }
	const double& z() const { return d_z; }
	const int& getType() const { return type;}
	const int& getNatom() const { return natom;}

//==========================================================
// MEMBER FUNCTION :
// EXPLAIN : set value
//==========================================================
	void setX(const double& r) { d_x = r; }
	void setY(const double& r) { d_y = r; }
	void setZ(const double& r) { d_z = r; }
	void setType(const int ntype) {type = ntype;}
	void setNatom(const int atom) {natom = atom;}
};

#endif
