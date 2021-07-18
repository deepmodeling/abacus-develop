#ifndef SYMMETRY_H
#define SYMMETRY_H

//#include "../src_pw/tools.h"
#include "../module_cell/unitcell_pseudo.h"
#include "symmetry_basic.h"

class Symmetry : public Symmetry_Basic
{
public:

	 Symmetry();
	~Symmetry();

	// mohan add 2021-02-21
	static bool symm_flag;

	void analy_sys(const UnitCell_pseudo &ucell, const output &out);
	bool available;

	Vector3<double> s1, s2, s3;
	Vector3<double> a1, a2, a3;	//primitive cell vectors(might be changed during the process of the program)
	Vector3<double>	p1, p2, p3;	//primitive cell vectors
	
	int ntype;	//the number of atomic species
	int nat; 	//the number of all atoms
 	int *na;	//number of atoms for each species
	int *istart; //start number of atom.
	int itmin_type; //the type has smallest number of atoms
	int itmin_start;

	// direct coordinates of atoms.
	double *dirpos;
	// cartesian coordinates of atoms.
	double *newpos;
	// positions of atoms after rotation.
	double *rotpos;
	
	
	double *ptrans;
    double ncell;	//the number of primitive cells within one supercell
	int *index;
	
	double cel_const[6];
	double pcel_const[6];
	int change; //whether the lattice vectors have been changed

	bool symflag_fft[48];
	int sym_test;
	int pbrav;
	int ibrav;
	string ilattname;	//the bravais lattice type of the supercell
	string plattname;	//the bravais lattice type of the primitive cell

	Matrix3 gmatrix[48];	//the rotation matrices for all space group operations
	Vector3<double> gtrans[48];
	
	Matrix3 symop[48];	//the rotation matrices for the pure bravais lattice
	int nop;	//the number of point group operations of the pure bravais lattice without basis
	int s_flag;	//whether the current matrix is one of all space group operations
	int nrot;	//the number of pure point group rotations
	int nrotk; 	//the number of all space group operations
	int pgnumber;	//the serial number of point group
	string pgname;	//the Schoenflies name of the point group

	int tab;

	int standard_lat(Vector3<double> &a,Vector3<double> &b,Vector3<double> &c,double *celconst );

	void lattice_type(Vector3<double> &v1,Vector3<double> &v2,Vector3<double> &v3, 
			int &ibrav,double *cel_const,string &bravname, const UnitCell_pseudo &ucell);

	void recip(
			const double a,
			const Vector3<double> &a1,
			const Vector3<double> &a2,
			const Vector3<double> &a3,
			Vector3<double> &b1,
			Vector3<double> &b2,
			Vector3<double> &b3
			);
	
	void change_lattice(void);

	// check if the input cell is a primitive cell.
	//void pricell(const UnitCell_pseudo &ucell);
	void getgroup(int &nrot, int &nrotk);
	void checksym(Matrix3 &s, Vector3<double> &gtrans, double *pos);
	void rho_symmetry(double *rho, const int &nr1, const int &nr2, const int &nr3);
	void force_symmetry(matrix &force, double* pos, const UnitCell_pseudo &ucell);
	void stress_symmetry(matrix &sigma, const UnitCell_pseudo &ucell);
	void write();

	void print_pos(const double* pos, const int &nat);


	private:

	// (s)tart (p)osition of atom (t)ype which
	// has (min)inal number.
	Vector3<double> sptmin;

};

#endif
