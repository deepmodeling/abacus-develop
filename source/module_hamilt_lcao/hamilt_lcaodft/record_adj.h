#ifndef RECORD_ADJ_H
#define RECORD_ADJ_H

#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_basis/module_ao/parallel_orbitals.h"


//---------------------------------------------------
// FUNCTION: record the adjacent atoms for each atom
//---------------------------------------------------
class Record_adj
{
	private:
		bool info_modified=false;
	public:

	Record_adj();
	~Record_adj();
	
	//--------------------------------------------
	// This will record the orbitals according to
	// HPSEPS's 2D block division.
	//--------------------------------------------
	void for_2d(Parallel_Orbitals &pv, bool gamma_only);

	//--------------------------------------------
	// This will record the orbitals according to
	// grid division (cut along z direction) 
	//--------------------------------------------
	void for_grid(const Grid_Technique &gt);

	void delete_grid();

	int na_proc;
	int* na_each;

	//--------------------------------------------
	// record sparse atom index in for_grid(const Grid_Technique &gt);
	// Map iat(dense atom index) to sparse atom index
	// Mainly removing the index dependency for OpenMP parallel loop
	// 
	// Meaning:
	// 1. if iat2ca[iat] > 0, it contains the sparse atom index
	// 2. if iat2ca[iat] < 0, the sparse atom index of iat does not exist 
	//
	// Usage: 
	// 1. iat2ca[iat] > 0 ? na_each[iat2ca[iat]] : 0
	// 2. iat2ca[iat] > 0 ? info[iat2ca[iat]] : nullptr
	//--------------------------------------------
	int* iat2ca;

	//------------------------------------------------
	// info will identify each atom in each unitcell.
	//------------------------------------------------
	int*** info;

	private:

};

#endif
