#ifndef RECORD_ADJ_H
#define RECORD_ADJ_H

#include "grid_technique.h"
#include "module_orbital/parallel_orbitals.h"


//---------------------------------------------------
// FUNCTION: record the adjacent atoms for each atom
//---------------------------------------------------
class Record_adj
{
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

	//------------------------------------------------
	// info will identify each atom in each unitcell.
	//------------------------------------------------
	int*** info;

	private:

};

#endif
