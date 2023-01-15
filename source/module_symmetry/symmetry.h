#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "../module_cell/unitcell.h"
#include "symmetry_basic.h"

namespace ModuleSymmetry
{
class Symmetry : public Symmetry_Basic
{
public:

	 Symmetry();
	~Symmetry();

	//symmetry flag for levels
	//-1 : no symmetry at all, k points would be total nks in KPT
	//0 : only basic time-reversal symmetry is considered, point k and -k would fold to k
	//1 : point group symmetry is considered
	static int symm_flag;

	void analy_sys(const UnitCell &ucell, std::ofstream &ofs_running);
	bool available;

	ModuleBase::Vector3<double> s1, s2, s3;
	ModuleBase::Vector3<double> a1, a2, a3;	//primitive cell vectors(might be changed during the process of the program)
	ModuleBase::Vector3<double>	p1, p2, p3;	//primitive cell vectors
	
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
	
	
	std::vector<ModuleBase::Vector3<double>> ptrans;
    int ncell;	//the number of primitive cells within one supercell
	int *index;
	
	double cel_const[6];
	double pcel_const[6];	//cel_const of primitive cell
	double pre_const[6];	//cel_const of input configuration
	int change; //whether the lattice vectors have been changed

	bool symflag_fft[48];
	int sym_test;
	int ibrav;
	int pbrav;		//ibrav of primitive cell
	int real_brav;    // the real ibrav for the cell     pengfei Li 3-15-2022
	std::string ilattname;	//the bravais lattice type of the supercell
	std::string plattname;	//the bravais lattice type of the primitive cell

	ModuleBase::Matrix3 gmatrix[48];	//the rotation matrices for all space group operations
	ModuleBase::Matrix3 kgmatrix[48];	//the rotation matrices in reciprocal space
	ModuleBase::Vector3<double> gtrans[48];
	
	ModuleBase::Matrix3 symop[48];	//the rotation matrices for the pure bravais lattice
	int nop;	//the number of point group operations of the pure bravais lattice without basis
	int s_flag;	//whether the current matrix is one of all space group operations
	int nrot;	//the number of pure point group rotations
	int nrotk; 	//the number of all space group operations
	int pgnumber;	//the serial number of point group
	int spgnumber;	//the serial number of point group in space group
	std::string pgname;	//the Schoenflies name of the point group R in {R|0}
	std::string spgname;	//the Schoenflies name of the point group R in the space group {R|t}

	ModuleBase::Matrix3 optlat;		//the optimized-symmetry lattice
	ModuleBase::Matrix3 plat;		//the primitive lattice

	int tab;

	int standard_lat(ModuleBase::Vector3<double> &a,ModuleBase::Vector3<double> &b,ModuleBase::Vector3<double> &c,double *celconst );

	void lattice_type(ModuleBase::Vector3<double> &v1,ModuleBase::Vector3<double> &v2,ModuleBase::Vector3<double> &v3, 
			double *cel_const,std::string &bravname, const UnitCell &ucell);

	void recip(
			const double a,
			const ModuleBase::Vector3<double> &a1,
			const ModuleBase::Vector3<double> &a2,
			const ModuleBase::Vector3<double> &a3,
			ModuleBase::Vector3<double> &b1,
			ModuleBase::Vector3<double> &b2,
			ModuleBase::Vector3<double> &b3
			);
	
	void change_lattice(void);

	// check if the input cell is a primitive cell.
	//void pricell(const UnitCell &ucell);
	void getgroup(int &nrot, int &nrotk, std::ofstream &ofs_running);
	void checksym(ModuleBase::Matrix3 &s, ModuleBase::Vector3<double> &gtrans, double *pos);
	void pricell(double* pos);
	void rho_symmetry(double *rho, const int &nr1, const int &nr2, const int &nr3);
	void rhog_symmetry(std::complex<double> *rhogtot, int* ixyz2ipw, const int &nx, 
			const int &ny, const int &nz, const int & fftnx, const int &fftny, const int &fftnz);
	void force_symmetry(ModuleBase::matrix &force, double* pos, const UnitCell &ucell);
	void stress_symmetry(ModuleBase::matrix &sigma, const UnitCell &ucell);
	void write();

	void print_pos(const double* pos, const int &nat);

	//convert n rotation-matrices from sa on basis {a1, a2, a3} to sb on basis {b1, b2, b3}
	void gmatrix_convert(const ModuleBase::Matrix3* sa, ModuleBase::Matrix3* sb, 
			const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b);
	void gmatrix_convert_int(const ModuleBase::Matrix3* sa, ModuleBase::Matrix3* sb, 
			const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b);
	//convert n translation-vectors from va on basis {a1, a2, a3} to vb on basis {b1, b2, b3}
	void gtrans_convert(const ModuleBase::Vector3<double>* va, ModuleBase::Vector3<double>* vb, 
			const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b);
	void gmatrix_invmap(const ModuleBase::Matrix3* s, const int n, int* invmap);
	private:

	// (s)tart (p)osition of atom (t)ype which
	// has (min)inal number.
	ModuleBase::Vector3<double> sptmin;

	// to be called in lattice_type and plat_type
	void get_shortest_latvec(ModuleBase::Vector3<double> &a1, 
			ModuleBase::Vector3<double> &a2, ModuleBase::Vector3<double> &a3);
	void get_optlat(ModuleBase::Vector3<double> &v1, ModuleBase::Vector3<double> &v2, 
			ModuleBase::Vector3<double> &v3, ModuleBase::Vector3<double> &w1, 
			ModuleBase::Vector3<double> &w2, ModuleBase::Vector3<double> &w3, 
			int &ibrav, int& real_brav, double* cel_const, double* tmp_const);
};
}

#endif
