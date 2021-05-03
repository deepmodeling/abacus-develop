#ifndef GRID_TECHNIQUE_H
#define GRID_TECHNIQUE_H

#include "grid_meshball.h"

// Author: mohan
// Date: 2009-10-17
class Grid_Technique : public Grid_MeshBall
{
	// public variables.
	public:
	
	//------------------------------------
	// 1: Info about atom number on grid.
	//------------------------------------
	// record how many atoms on each grid.
	int* how_many_atoms;
	// max atom on grid
	int max_atom;
	// sum of how_many_atoms
	int total_atoms_on_grid;

	//------------------------------------
	// 2: Info about which atom on grid.
	//------------------------------------
	// save the start position of each big cell's adjacent
	// atoms in 1D grid.
	int* bcell_start;
	// save the 'iat' atom. 
	// dim: total_atoms_on_grid.
	int* which_atom;

	//--------------------------------------
	// save the bigcell index in meshball.
	// dim: total_atoms_on_grid.
	//--------------------------------------
	int* which_bigcell;
	int* which_unitcell;
	
	//------------------------------------
	// 3: which atom on local grid.
	//------------------------------------
	bool* in_this_processor;
	int lnat; // local nat.
	int lgd; // local grid dimension.  lgd * lgd symmetry matrix. 
	int* trace_lo; // trace local orbital.
	int lgbeta; // mohan add
	int* trace_beta; // sunzhiyuan add, trace to nonlocal projector beta.

	int* atomip; // atom index in this processor
	
	// public functions
	public:

	Grid_Technique();
	~Grid_Technique();

	void set_pbc_grid(
			const int &ncx_in,
			const int &ncy_in,
			const int &ncz_in,
			const int &bx_in,
			const int &by_in,
			const int &bz_in,
			const int &nbx_in,
			const int &nby_in,
			const int &nbz_in,
			const int &nbxx_in,
			const int &nbzp_start_in,
			const int &nbzp_in);

	private:
		
	// atoms on meshball
	void init_atoms_on_grid(void);
	void init_atoms_on_grid2(const int* index2normal);
	void cal_grid_integration_index(void);
	void cal_trace_lo(void);
	void cal_trace_beta(void);//mohan add 2012-04-13
};

extern Grid_Technique GridT;

#endif
