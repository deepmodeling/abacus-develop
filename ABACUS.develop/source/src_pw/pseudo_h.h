#ifndef PSEUDOH_H
#define PSEUDOH_H

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#include "pseudopot_upf.h"

class pseudo_h
{
public:

	pseudo_h();
	~pseudo_h();

	//<PP_HEADER>
	bool has_so;		// if .true. includes spin-orbit
	int  nv;			// UPF file version number
	string psd;			// Element label
	string pp_type;		// Pseudo type ( NC or US )
	bool tvanp;			// .true. if Ultrasoft
	bool nlcc;			// Non linear core corrections(bool)
	string dft[4];		// Exch-Corr type
	int  zv;			// z valence
	double etotps;		// total energy
	double ecutwfc;		// suggested cut-off for wfc
	double ecutrho;		// suggested cut-off for rho
	int lmax;			// maximum angular momentum component
	int mesh;			// number of point in the radial mesh
	int nchi;			// nwfc,number of wavefunctions
	int nbeta;			// number of projectors
	string *els;		// els[nchi]
	int *lchi;			// lchi[nchi]
	double *oc;			// oc[nchi]

	double *jjj;	  	// total angual momentum, jjj[nbeta]
	double *jchi;		//jchi(nwfc), added by zhengdy-soc
	int *nn;

	// member functions
	void set_pseudo_h(const Pseudopot_upf &upf);
	void print_pseudo_h(ofstream &ofs);
};

#endif // PSEUDOH_H
