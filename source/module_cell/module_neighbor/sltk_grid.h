#ifndef GRID_H
#define GRID_H

#include <stdexcept>
#include <functional>
#include "sltk_util.h"
#include "sltk_atom.h"
#include "sltk_atom_input.h"

#include "module_cell/unitcell.h"
//extern Structure_Factor sf;

//==========================================================
// STRUCT :
// NAME : AtomLink
// NAME : CellSet
//==========================================================

struct AtomLink
{
	FAtom fatom;
	AtomLink* next_p;

	// Constructors and destructor
	AtomLink
	(
	    const FAtom& atom = FAtom(),
	    AtomLink* const pointNext = NULL //mohan fix bug 2011/09/26, from NullPtr->NULL
	);

};

struct CellSet
{
	AtomLink* address;
	int length;
	int in_grid[3];
	CellSet();
};

//==========================================================
// CLASS NAME :
// Atom_input : defined elsewhere
//==========================================================

class Atom_input;

//==========================================================
// CLASS NAME :
// Grid :
//==========================================================

class Grid
{
public:

	// Constructors and destructor
	Grid(const int &test_grid_in);
	virtual ~Grid();

	void init(
		std::ofstream &ofs,
		const UnitCell &ucell, 
		const Atom_input &input);

	//2015-05-07
	void delete_vector(const Atom_input &input);


	//Static data
	static const double TOLERATE_ERROR;
	static const std::hash<int> INT_HASHER;
	static const char* const ERROR[3];

	//Data
	int natom;// Total atoms.
	bool pbc; // periodic boundary condition
	bool expand_flag;
	double sradius;// searching radius
	double d_minX;// origin of all cells
	double d_minY;
	double d_minZ;
	int dx;
	int dy;
	int dz;
	int layer;
	double cell_x_length;
	double cell_y_length;
	double cell_z_length;
	CellSet ***Cell; //dx , dy ,dz is cell number in each direction,respectly.
	void delete_Cell() //it will replace by container soon!
	{
		if (this->init_cell_flag)
		{
			for (int i = 0;i < this->dx;i++)
			{
				for (int j = 0;j < this->dy;j++)
				{
					delete[] this->Cell[i][j];
				}
			}

			for (int i = 0;i < this->dx;i++)
			{
				delete[] this->Cell[i];
			}

			delete[] this->Cell;
			this->init_cell_flag = false;
		}
	}

	double grid_length[3];
	double vec1[3];
	double vec2[3];
	double vec3[3];
	double lat_now;
	bool init_cell_flag;
    //LiuXh add 2019-07-15
    const double& getD_minX(void) const {return d_minX;}
    const double& getD_minY(void) const {return d_minY;}
    const double& getD_minZ(void) const {return d_minZ;}

    const int& getCellX(void) const {return dx;}
    const int& getCellY(void) const {return dy;}
    const int& getCellZ(void) const {return dz;}

	// Inner Function
protected:
	AtomLink* getHashCode(const UnitCell &ucell, const FAtom &atom)const;		// Peize Lin delete const 2018-07-14
//	AtomLink* const getHashCode(const FAtom &atom)const;
	AtomLink* atomlink;
	AtomLink* cordon_p;// Warning! A guard! Don't delete it!

private:

	const int test_grid;
//==========================================================
// MEMBER FUNCTIONS :
// Three Main Steps:
// NAME : setMemberVariables (read in datas from Atom_input,
// 			init cells.)
// NAME : setAtomLinkArray( set the AtomLinkArray twice,
// 			first use Build_Hash,second use Fold_Hash)
// NAME : setBoundaryAdjacent( Consider different situations,
// 			if not_expand case : nature/periodic boundary
// 			condition , if expand_case)
//==========================================================
	void setMemberVariables(
		std::ofstream &ofs_in, 
		const Atom_input &input);

	void setAtomLinkArray(
		const UnitCell &ucell, 
		const Atom_input &input);

	void setBoundaryAdjacent(
		std::ofstream &ofs_in, 
		const Atom_input &input);

//==========================================================
//
//==========================================================
	AtomLink* Build_Cache(const UnitCell &ucell, const Atom_input &input);		
	// Peize Lin delete const and throw(std::out_of_range, std::logic_error) 2018-07-14

	//	AtomLink* const Build_Cache(const Atom_input &input) throw(std::out_of_range, std::logic_error);
	bool Push(const UnitCell &ucell, const FAtom& atom);
	void In_Which_Cell(const UnitCell &ucell, int &a, int &b, int &c, const FAtom &atom)const;
	void Build_Cell(void);
	void Build_Hash_Table(const UnitCell &ucell, AtomLink* const pointCache);
	void Fold_Hash_Table(void);		// Peize Lin delete const and throw(std::logic_error) 2018-07-14
	static int Hash_one_hit;

//==========================================================
//
//==========================================================
	void Construct_Adjacent_expand(const int i, const int j, const int k);
	void Construct_Adjacent_expand_periodic(
	    const int i, const int j, const int k, const int ia);

	void Construct_Adjacent_begin(void);
	void Construct_Adjacent_nature(
	    const int i, const int j, const int k, const int ia);
	void Construct_Adjacent_periodic(
	    const int i, const int j, const int k, const int ia);
	void Construct_Adjacent_final(
	    const int i, const int j, const int k, const int ia,
	    const int i2, const int j2, const int k2, const int ia2);
};

#endif
