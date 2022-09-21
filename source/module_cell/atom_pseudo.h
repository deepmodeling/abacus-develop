#ifndef ATOM_PSEUDO_H
#define ATOM_PSEUDO_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/vector3.h"
#include "../src_io/output.h"
#include "../module_base/complexarray.h"
#include "../module_base/complexmatrix.h"
#include "pseudo_nc.h"
using namespace std;

class Atom_pseudo : public pseudo_nc
{
public:

	Atom_pseudo();
	~Atom_pseudo();

	ModuleBase::Vector3<int> *mbl; // whether the atoms can move or not
	std::string pseudo_fn; // File name of pseudopotentials
	std::string pseudo_type; // Type of pseudopotential. It is added to be on the safe side, thought it seems useless. sunliang 2022-09-15
	double mass; // the mass of atom
	bool flag_empty_element = false; // whether is the empty element for bsse.	Peize Lin add 2021.04.07

	// mohan add 2021-05-07
	ModuleBase::ComplexArray d_so; //(:,:,:), spin-orbit case
	int nproj;
	int nproj_soc; // dimension of D_ij^so
	int non_zero_count_soc[4];
	int *index1_soc[4];
	int *index2_soc[4];

	void set_d_so( // mohan add 2021-05-07
		ModuleBase::ComplexMatrix &d_so_in,
		const int &nproj_in,
		const int &nproj_in_so,
		const bool has_so);
	
protected:

	void print_atom(std::ofstream &ofs);

#ifdef __MPI
	void bcast_atom_pseudo(const int &na);
	void bcast_atom_pseudo2(void); // for upf201
#endif

};

#endif
