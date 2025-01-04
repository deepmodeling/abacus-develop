#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "module_cell/unitcell_data.h"
#include "module_cell/atom_spec.h"
#include "symmetry_basic.h"

namespace ModuleSymmetry
{
class Symmetry : public Symmetry_Basic
{
public:
    Symmetry() {
        this->epsilon = 1e-6;
    };
    ~Symmetry() {};

	//symmetry flag for levels
	//-1 : no symmetry at all, k points would be total nks in KPT
	//0 : only basic time-reversal symmetry is considered, point k and -k would fold to k
	//1 : point group symmetry is considered
    static int symm_flag;
    static bool symm_autoclose; // controled by INPUT
    static bool pricell_loop;   ///< whether to loop primitive cell in rhog_symmetry, Only for AFM

    /// @brief analyze the symmetry of the system
    /// @param lat structure of lattice
    /// @param st 
    /// @param atoms all atoms
    /// @param ofs_running 
	/// get the symmetry information of the system, gmatries (rotation 3*3 matrixs), gtrans (transfer a collections vector3), etc.
    void analy_sys(const Lattice& lat, const Statistics& st, Atom* atoms, std::ofstream& ofs_running);

	ModuleBase::Vector3<double> s1, s2, s3;
	ModuleBase::Vector3<double> a1, a2, a3;	//primitive cell vectors(might be changed during the process of the program)
	ModuleBase::Vector3<double>	p1, p2, p3;	//primitive cell vectors
	
	int ntype=0;	  //the number of atomic species
	int nat  =0; 	  //the number of all atoms
 	int *na  =nullptr;//number of atoms for each species
	int *istart=nullptr; //start number of atom.
	int itmin_type=0; //the type has smallest number of atoms
	int itmin_start=0;

	// direct coordinates of atoms.
	double *newpos=nullptr;
	// positions of atoms after rotation.
	double *rotpos=nullptr;
	
	
	std::vector<ModuleBase::Vector3<double>> ptrans; // the translation vectors of the primitive cell in the input structure
    int ncell=1;	//the number of primitive cells within one supercell
	int *index=nullptr;
	
	double cel_const[6]={0.0};
	double pcel_const[6]={0.0};	//cel_const of primitive cell
	double pre_const[6]={0.0};	//cel_const of input configuration, first 3 is moduli of a1, a2, a3, last 3 is eular angle

	bool symflag_fft[48]={false};
	int sym_test=0;
	int pbrav=0;		//ibrav of primitive cell
	int real_brav=0;    // the real ibrav for the cell     pengfei Li 3-15-2022
	std::string ilattname;	//the bravais lattice type of the supercell
	std::string plattname;	//the bravais lattice type of the primitive cell

	ModuleBase::Matrix3 gmatrix[48];	//the rotation matrices for all space group operations
	ModuleBase::Matrix3 kgmatrix[48];	//the rotation matrices in reciprocal space
	ModuleBase::Vector3<double> gtrans[48];
	
	ModuleBase::Matrix3 symop[48];	//the rotation matrices for the pure bravais lattice
    int nop=0;	//the number of point group operations of the pure bravais lattice without basis
	int nrot=0;	//the number of pure point group rotations
    int nrotk = -1; 	//the number of all space group operations, >0 means the nrotk has been analyzed
    int max_nrotk = -1;  ///< record the maximum number of symmetry operations during cell-relax
    int pgnumber=0;	//the serial number of point group
	int spgnumber=0;	//the serial number of point group in space group
	std::string pgname;	//the Schoenflies name of the point group R in {R|0}
	std::string spgname;	//the Schoenflies name of the point group R in the space group {R|t}

	ModuleBase::Matrix3 optlat;		//the optimized-symmetry lattice
	ModuleBase::Matrix3 plat;		//the primitive lattice

    bool all_mbl = true;    ///< whether all the atoms are movable in all the directions

    int standard_lat(ModuleBase::Vector3<double>& a, ModuleBase::Vector3<double>& b, ModuleBase::Vector3<double>& c, double* celconst)const;

	void lattice_type(ModuleBase::Vector3<double> &v1,ModuleBase::Vector3<double> &v2,ModuleBase::Vector3<double> &v3, 
        ModuleBase::Vector3<double>& v01, ModuleBase::Vector3<double>& v02, ModuleBase::Vector3<double>& v03,
        double* cel_const, double* pre_const, int& real_brav, std::string& bravname, const Atom* atoms,
        bool convert_atoms, double* newpos = nullptr)const;

    void getgroup(int& nrot, int& nrotk, std::ofstream& ofs_running, const int& nop,
        const ModuleBase::Matrix3* symop, ModuleBase::Matrix3* gmatrix, ModuleBase::Vector3<double>* gtrans,
        double* pos, double* rotpos, int* index, const int itmin_type, const int itmin_start, int* istart, int* na)const;
    bool checksym(const ModuleBase::Matrix3& s, ModuleBase::Vector3<double>& gtrans,
        double* pos, double* rotpos, int* index, const int itmin_type, const int itmin_start, int* istart, int* na)const;
    /// @brief  primitive cell analysis
    void pricell(double* pos, const Atom* atoms);

	/// -----------------------
	/// Symmetrize the charge density, the forces, and the stress
	/// -----------------------
	void rho_symmetry(double *rho, const int &nr1, const int &nr2, const int &nr3);
	void rhog_symmetry(std::complex<double> *rhogtot, int* ixyz2ipw, const int &nx, 
			const int &ny, const int &nz, const int & fftnx, const int &fftny, const int &fftnz);
    /// symmetrize a vector3 with nat elements, which can be forces or variation of atom positions in relax
    void symmetrize_vec3_nat(double* v)const;   // force
    /// symmetrize a 3*3 tensor, which can be stress or variation of unitcell in cell-relax
    void symmetrize_mat3(ModuleBase::matrix& sigma, const Lattice& lat)const; // stress

	//convert n rotation-matrices from sa on basis {a1, a2, a3} to sb on basis {b1, b2, b3}
	void gmatrix_convert(const ModuleBase::Matrix3* sa, ModuleBase::Matrix3* sb, 
			const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b)const;
	void gmatrix_convert_int(const ModuleBase::Matrix3* sa, ModuleBase::Matrix3* sb, 
			const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b)const;
	//convert n translation-vectors from va on basis {a1, a2, a3} to vb on basis {b1, b2, b3}
	void gtrans_convert(const ModuleBase::Vector3<double>* va, ModuleBase::Vector3<double>* vb, 
			const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b)const;
	void gmatrix_invmap(const ModuleBase::Matrix3* s, const int n, int* invmap) const;
	void hermite_normal_form(const ModuleBase::Matrix3 &s, ModuleBase::Matrix3 &H, ModuleBase::Matrix3 &b) const;
    int get_rotated_atom(int isym, int iat)const
    {
        if (!this->isym_rotiat_.empty()) { return this->isym_rotiat_[isym][iat]; }
        else { return -1; }
    }
	private:

    /// atom-map for each symmetry operation: isym_rotiat[isym][iat]=rotiat
    std::vector<std::vector<int>> isym_rotiat_;

    /// @brief  set atom map for each symmetry operation
    void set_atom_map(const Atom* atoms);
    /// @brief check if all the atoms are movable
    ///  delta_pos symmetrization in relax is only meaningful when all the atoms are movable in all the directions.
    bool is_all_movable(const Atom* atoms, const Statistics& st)const;

    // to be called in lattice_type
	void get_shortest_latvec(ModuleBase::Vector3<double> &a1, 
			ModuleBase::Vector3<double> &a2, ModuleBase::Vector3<double> &a3)const;
	void get_optlat(ModuleBase::Vector3<double> &v1, ModuleBase::Vector3<double> &v2, 
			ModuleBase::Vector3<double> &v3, ModuleBase::Vector3<double> &w1, 
			ModuleBase::Vector3<double> &w2, ModuleBase::Vector3<double> &w3, 
        int& real_brav, double* cel_const, double* tmp_const)const;

    /// Loop the magmom of each atoms in its type when NSPIN>1. If not all the same, primitive cells should not be looped in rhog_symmetry.
    bool magmom_same_check(const Atom* atoms)const;

};
}

#endif
