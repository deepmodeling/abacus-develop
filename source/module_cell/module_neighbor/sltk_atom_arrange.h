#ifndef ATOM_ARRANGE_H
#define ATOM_ARRANGE_H

#include "sltk_grid.h"
#include "sltk_grid_driver.h"
#include "sltk_atom_input.h"


class atom_arrange
{
public:

	atom_arrange();
	~atom_arrange();
	
	static void search(
		const bool flag,
		std::ofstream &ofs,
		Grid_Driver &grid_d, 
		const UnitCell &ucell, 
		const double& search_radius_bohr, 
		const int &test_atom_in,
		const bool test_only = false);

	//caoyu modify 2021-05-24
	static double set_sr_NL(
		std::ofstream &ofs_in,
		std::string &output_level,
		const double& rcutmax_Phi, 
		const double& rcutmax_Beta, 
		const bool gamma_only_local);

	//2015-05-07
	static void delete_vector(
		std::ofstream &ofs_in,
		const bool pbc_flag,
		Grid_Driver &grid_d, 
		const UnitCell &ucell, 
		const double &search_radius_bohr, 
		const int &test_atom_in);

private:

};

#endif
