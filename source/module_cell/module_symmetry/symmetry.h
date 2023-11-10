#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "module_cell/unitcell.h"
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
    static bool symm_autoclose;
    static bool pricell_loop;   ///< whether to loop primitive cell in rhog_symmetry

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
	double *newpos;
	// positions of atoms after rotation.
	double *rotpos;
	
	
	std::vector<ModuleBase::Vector3<double>> ptrans;
    int ncell=1;	//the number of primitive cells within one supercell
	int *index;
	
	double cel_const[6];
	double pcel_const[6];	//cel_const of primitive cell
	double pre_const[6];	//cel_const of input configuration

	bool symflag_fft[48];
	int sym_test;
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
    int nrotk = -1; 	//the number of all space group operations
    int max_nrotk = -1;  ///< record the maximum number of symmetry operations during cell-relax
    int pgnumber;	//the serial number of point group
	int spgnumber;	//the serial number of point group in space group
	std::string pgname;	//the Schoenflies name of the point group R in {R|0}
	std::string spgname;	//the Schoenflies name of the point group R in the space group {R|t}

	ModuleBase::Matrix3 optlat;		//the optimized-symmetry lattice
	ModuleBase::Matrix3 plat;		//the primitive lattice

	int tab;

	int standard_lat(ModuleBase::Vector3<double> &a,ModuleBase::Vector3<double> &b,ModuleBase::Vector3<double> &c,double *celconst )const;

	void lattice_type(ModuleBase::Vector3<double> &v1,ModuleBase::Vector3<double> &v2,ModuleBase::Vector3<double> &v3, 
	    	ModuleBase::Vector3<double> &v01, ModuleBase::Vector3<double> &v02, ModuleBase::Vector3<double> &v03,
			double *cel_const, double* pre_const, int& real_brav, std::string &bravname, const UnitCell &ucell, 
			bool convert_atoms, double* newpos=nullptr)const;

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

    /// symmetrize a vector3 with nat elements, which can be forces or variation of atom positions in relax
    void symmetrize_vec3_nat(double* v)const;
    /// symmetrize a 3*3 tensor, which can be stress or variation of unitcell in cell-relax
    void symmetrize_mat3(ModuleBase::matrix& sigma, const UnitCell& ucell)const;

    void write();

	void print_pos(const double* pos, const int &nat);

	//convert n rotation-matrices from sa on basis {a1, a2, a3} to sb on basis {b1, b2, b3}
	void gmatrix_convert(const ModuleBase::Matrix3* sa, ModuleBase::Matrix3* sb, 
			const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b)const;
	void gmatrix_convert_int(const ModuleBase::Matrix3* sa, ModuleBase::Matrix3* sb, 
			const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b)const;
	//convert n translation-vectors from va on basis {a1, a2, a3} to vb on basis {b1, b2, b3}
	void gtrans_convert(const ModuleBase::Vector3<double>* va, ModuleBase::Vector3<double>* vb, 
			const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b)const;
	void gmatrix_invmap(const ModuleBase::Matrix3* s, const int n, int* invmap);
	void hermite_normal_form(const ModuleBase::Matrix3 &s, ModuleBase::Matrix3 &H, ModuleBase::Matrix3 &b) const;
	private:

	// (s)tart (p)osition of atom (t)ype which
	// has (min)inal number.
	ModuleBase::Vector3<double> sptmin;

    /// atom-map for each symmetry operation: isym_rotiat[isym][iat]=rotiat
    std::vector<std::vector<int>> isym_rotiat_;


    /// @brief  set atom map for each symmetry operation
    void set_atom_map(const UnitCell& ucell);
    // to be called in lattice_type
	void get_shortest_latvec(ModuleBase::Vector3<double> &a1, 
			ModuleBase::Vector3<double> &a2, ModuleBase::Vector3<double> &a3)const;
	void get_optlat(ModuleBase::Vector3<double> &v1, ModuleBase::Vector3<double> &v2, 
			ModuleBase::Vector3<double> &v3, ModuleBase::Vector3<double> &w1, 
			ModuleBase::Vector3<double> &w2, ModuleBase::Vector3<double> &w3, 
        int& real_brav, double* cel_const, double* tmp_const)const;

    /// Loop the magmom of each atoms in its type when NSPIN>1. If not all the same, primitive cells should not be looped in rhog_symmetry.
    bool magmom_same_check(const UnitCell& ucell)const;
};
}

#endif
