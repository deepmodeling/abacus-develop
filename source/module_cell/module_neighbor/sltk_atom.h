#ifndef INCLUDE_FATOM
#define INCLUDE_FATOM

#include <memory>
#include "sltk_util.h"
#include "module_base/timer.h"
#include <vector>

// a class contains the atom position, 
// the type and the index,
class FAtom
{
private:
	double d_x;
	double d_y;
	double d_z;
	std::vector<FAtom *> adjacent;

	int type;
	int natom;

	int cell_x;
	int cell_y;
	int cell_z;
public:
//==========================================================
// Default Constructor and deconstructor
//==========================================================

	FAtom();
	FAtom(const double& x_in, const double& y_in, const double& z_in, 
			const int& type_in, const int& natom_in, 
			const int& cell_x_in, const int& cell_y_in, const int& cell_z_in)
	{
		d_x = x_in;
		d_y = y_in;
		d_z = z_in;
		type = type_in;
		natom = natom_in;
		cell_x = cell_x_in;
		cell_y = cell_y_in;
		cell_z = cell_z_in;
	}
	~FAtom()
	{
		adjacent.clear();
	}

	void addAdjacent(FAtom& atom_in)
	{
		adjacent.push_back( &atom_in);
	}
	const std::vector<FAtom *>& getAdjacent() const { return adjacent; }
	void clearAdjacent() { adjacent.clear(); }
//==========================================================
// MEMBER FUNCTION :
// EXPLAIN : get value
//==========================================================
	const double& x() const { return d_x; }
	const double& y() const { return d_y; }
	const double& z() const { return d_z; }
	const int& getType() const { return type;}
	const int& getNatom() const { return natom;}
	const int& getCellX() const { return cell_x; }
	const int& getCellY() const { return cell_y; }
	const int& getCellZ() const { return cell_z; }
};

#endif
