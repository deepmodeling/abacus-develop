#ifndef ATOM_H
#define ATOM_H

#include "atom_pseudo.h"
#include "../src_io/output.h"
class Atom: public Atom_pseudo
{
public:

    // constructor and destructor
    Atom();
    ~Atom();

    int *iw2m; // use iw to find m
    int *iw2n; // use iw to find n
    int *iw2l; // use iw to find L
	int *iw2_ylm;
	bool *iw2_new;
    int nw; // number of local orbitals (l,n,m) of this type

    void set_index(void);

    int type; // Index of atom type
    int na; // Number of atoms in this type.

    int nwl; // max L(Angular momentum) (for local basis)
    double Rcut; //pengfei Li 16-2-29
    int *l_nchi; // number of chi for each L

    int stapos_wf; // start position of wave functions

    std::string label; // atomic symbol
    ModuleBase::Vector3<double> *tau;// Cartesian coordinates of each atom in this type.
    ModuleBase::Vector3<double> *tau_original;// Cartesian coordinates of each atom in this type, but without periodic adjustment.
	ModuleBase::Vector3<double> *taud;// Direct coordinates of each atom in this type.
    ModuleBase::Vector3<double> *vel;// velocities of each atom in this type.
    ModuleBase::Vector3<double> *force; // force acting on each atom in this type.

	double* mag;
	double* angle1;//spin angle, added by zhengdy-soc
	double* angle2;
    ModuleBase::Vector3<double> *m_loc_;

    void print_Atom(std::ofstream &ofs);
    void update_force(ModuleBase::matrix &fcs);
#ifdef __MPI
    void bcast_atom(void);
    void bcast_atom2(void);
#endif

};

#endif //Atomspec

