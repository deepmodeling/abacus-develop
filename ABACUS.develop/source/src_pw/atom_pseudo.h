#ifndef ATOM_PSEUDO_H
#define ATOM_PSEUDO_H

#include "tools.h"
#include "pseudo_nc.h"
using namespace std;

class Atom_pseudo : public pseudo_nc
{
public:

	Atom_pseudo();
	~Atom_pseudo();

	Vector3<int> *mbl; // whether the atoms can move or not
	string pseudo_fn; // File name of pseudopotentials
	double mass; // the mass of atom
	bool flag_empty_element = false; // whether is the empty element for bsse.	Peize Lin add 2021.04.07

protected:

	void print_atom(ofstream &ofs);

#ifdef __MPI
	void bcast_atom_pseudo(const int &na);
	void bcast_atom_pseudo2(void); // for upf201
#endif

};

#endif
