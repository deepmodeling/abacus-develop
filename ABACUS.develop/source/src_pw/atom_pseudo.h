//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-08
//==========================================================
#ifndef ATOM_PSEUDO_H
#define ATOM_PSEUDO_H

#include "tools.h"
#include "pseudo_us.h"
using namespace std;

//==========================================================
// CLASS :
// NAME : Atom_pseudo
//==========================================================
class Atom_pseudo : public pseudo_us
{
public:

	Atom_pseudo();
	~Atom_pseudo();

	Vector3<int> *mbl; //If this atom can move
	string pseudo_fn;// File name of pseudopotentia
	double mass; // the mass of atom

protected:

	void print_atom(ofstream &ofs);

#ifdef __MPI
	void bcast_atom_pseudo(const int &na);
	void bcast_atom_pseudo2(void);
#endif

};

#endif
