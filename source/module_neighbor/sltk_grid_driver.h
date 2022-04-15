#ifndef GRID_DRIVER_H
#define GRID_DRIVER_H

#include <memory>
#include <stdexcept>
#include "sltk_atom.h"
#include "sltk_atom_input.h"
#include "sltk_grid.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/vector3.h"
#include "src_pw/pw_basis.h"
#include "module_cell/unitcell_pseudo.h"
#include <tuple>

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
		const ModuleBase::Vector3<double> &cartesian_posi, 
		const int &ntype, 
		const int &nnumber);

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
	const ModuleBase::Vector3<double>& getAdjacentTau(const int i) const { return adjacent_tau[i]; } 
	const ModuleBase::Vector3<int>& getBox(const int i) const {return box[i];}

    std::vector<std::tuple<int, int, ModuleBase::Vector3<int>, ModuleBase::Vector3<double>>> get_adjs(const UnitCell_pseudo& ucell_in, const size_t &iat);
    std::vector<std::vector<std::tuple<int, int, ModuleBase::Vector3<int>, ModuleBase::Vector3<double>>>> get_adjs(const UnitCell_pseudo& ucell_in);

private:

	mutable int adj_num;
	std::vector<int> ntype;
	std::vector<int> natom;
	std::vector<ModuleBase::Vector3<double>> adjacent_tau;
	std::vector<ModuleBase::Vector3<int>> box;

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
		const ModuleBase::Vector3<double> &cartesian_pos, 
		const int &ntype, 
		const int &nnumber)const;

	void Find_adjacent_atom(
		const int offset, 
		std::shared_ptr<AdjacentSet> as);

	double Distance(const AtomLink& a1, const ModuleBase::Vector3<double> &adjacent_site)const;
	double Distance(const AtomLink& a1, const AtomLink& a2)const;

//==========================================================
// MEMBER FUNCTIONS :
// NAME : Calculate_adjacent_site
//==========================================================
	ModuleBase::Vector3<double> Calculate_adjacent_site
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
