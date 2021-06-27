#ifndef UNITCELL_PSEUDO_H
#define UNITCELL_PSEUDO_H

#include "atom_spec.h"
#include "src_pw/tools.h"
#include "unitcell.h"

class UnitCell_pseudo : public UnitCell
{
public: // member variables
	//============================================================
	// meshx : max number of mesh point in pseudopotential file
	// natomwfc : number of starting wavefunctions
	// lmax  : Max L used for localized orbital.
	// nmax  : Max N used for localized orbital.
	// lmax_ppwf : Max L of pseudo wave functinos
	// nelec : total number of electrons
	//============================================================
	int meshx;
	int natomwfc;
	int lmax;
	int nmax;
	int nmax_total;//mohan add 2009-09-10
	int lmax_ppwf;
	double nelec;

public: // member functions
	UnitCell_pseudo();
	~UnitCell_pseudo();
	void setup_cell(const string &s_pseudopot_dir, const string &fn ,ofstream &log);
	void read_atom_species(ifstream &ifa); // read in the atom information for each type of atom
	bool read_atom_positions(ifstream &ifpos); // read in atomic positions
	int find_type(const string &label);
	void print_tau(void)const;
	void print_stru_file(const string &fn, const int &type=1)const; // mohan add 2011-03-22
	void check_dtau(void);
    void setup_cell_after_vc(const string &s_pseudopot_dir, const string &fn, ofstream &log); //LiuXh add 20180515

	bool set_atom_flag;//added on 2009-3-8 by mohan

	// read in pseudopotential from files for each type of atom
	void read_cell_pseudopots(const string &fn);

	//================================================================
	// cal_natomwfc : calculate total number of atomic wavefunctions
	// cal_nwfc     : calculate total number of local basis and lmax
	// cal_nelec    : calculate total number of electrons
	// cal_meshx	: calculate max number of mesh points in pp file
	//================================================================
	void cal_nwfc();
	void cal_nelec();
	void cal_meshx();
	void cal_natomwfc(); 
	void print_unitcell_pseudo(const string &fn);
	bool check_tau(void)const; //mohan add 2011-03-03

#ifdef __MPI
	void bcast_unitcell_pseudo(void);
	void bcast_unitcell_pseudo2(void);
#endif
};

#endif
