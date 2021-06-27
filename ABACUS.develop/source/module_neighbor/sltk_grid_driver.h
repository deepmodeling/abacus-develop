#ifndef GRID_DRIVER_H
#define GRID_DRIVER_H

#include <stdexcept>
#include "sltk_atom.h"
#include "sltk_atom_input.h"
#include "sltk_grid.h"
#include "../src_pw/tools.h"
#include "../src_pw/pw_basis.h"

class Grid_Driver : public Grid
{
public:

	//==========================================================
	// THE INTERFACE WITH USER :
	// MEMBRE FUNCTIONS :
	// NAME : Find_atom (input cartesian position,find the
	//		adjacent of this atom,and store the information
	//		in 'adj_num','ntype','natom'
	//==========================================================
	Grid_Driver(
		const int &test_d_in, 
		const int &test_gd_in, 
		const int &test_grid_in);

	~Grid_Driver();

	void Find_atom(
		const UnitCell &ucell, 
		const Vector3<double> &cartesian_posi, 
		const int &ntype, 
		const int &nnumber)const;

	//==========================================================
	// EXPLAIN : The adjacent information for the input 
	// cartesian_pos
	// MEMBER VARIABLES : 
	// NAME : getAdjacentNum
	// NAME : getNtype
	// NAME : getNatom
	// NAME : getAdjaentTau
	//==========================================================
	const int& getAdjacentNum(void)const { return adj_num; }
	const int& getType(const int i) const { return ntype[i]; }
	const int& getNatom(const int i) const { return natom[i]; }
	const Vector3<double>& getAdjacentTau(const int i) const { return adjacent_tau[i]; } 
	const Vector3<int>& getBox(const int i) const {return box[i];}

private:

	mutable int adj_num;
	mutable int* ntype;
	mutable int* natom;
	mutable Vector3<double>* adjacent_tau;
	mutable Vector3<int>* box;

	const int test_deconstructor;//caoyu reconst 2021-05-24
	const int test_grid_driver;	//caoyu reconst 2021-05-24

	//==========================================================
	// MEMBER FUNCTIONS :
	// NAME : Locate_offset (find the atom index according to pos)
	// NAME : Find_adjacent_atom ( find adjacent atmos for offset)
	// NAME : Distance ( between a1 and a2)
	//==========================================================
	int Locate_offset(
		const UnitCell &ucell, 
		const Vector3<double> &cartesian_pos, 
		const int &ntype, 
		const int &nnumber)const;

	void Find_adjacent_atom(
		const int offset, 
		AdjacentSet* as)const;

	double Distance(const AtomLink& a1, const Vector3<double> &adjacent_site)const;
	double Distance(const AtomLink& a1, const AtomLink& a2)const;

//==========================================================
// MEMBER FUNCTIONS :
// NAME : Calculate_adjacent_site
//==========================================================
	Vector3<double> Calculate_adjacent_site
	(
	    const short offset, // use offset cartesian coordiante
	    const double &box11, const double &box12, const double &box13,
	    const double &box21, const double &box22, const double &box23,
	    const double &box31, const double &box32, const double &box33,
	    const short box_x, // three dimensions of the target box
	    const short box_y,
	    const short box_z
	)const;
};

#endif
