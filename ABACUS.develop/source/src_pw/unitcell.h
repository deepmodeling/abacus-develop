#ifndef UNITCELL_H
#define UNITCELL_H

#include "tools.h"
#include "atom_spec.h"

// the "base class" of UnitCell_pseudo.
// provide the basic information about unitcell.
class UnitCell
{
public:
    Atom *atoms;

    int ntype;// number of atom species in UnitCell
    int nat; // total number of atoms of all species in unitcell
    string Coordinate; // "Direct" or "Cartesian" or "Cartesian_angstrom"
    string latName; // Lattice name
    double lat0; // Lattice constant(bohr)(a.u.)
    double lat0_angstrom;// Lattice constant(angstrom)
    double tpiba;// 2*pi / lat0;
    double tpiba2; // tpiba ^ 2
    double omega;// the volume of the unit cell

    Matrix3 latvec; // Unitcell lattice vectors
	int *lc;  // Change the lattice vectors or not
	Vector3<double> a1,a2,a3; // Same as latvec, just at another form.
	Vector3<double> latcenter; // (a1+a2+a3)/2 the center of vector
    Matrix3 latvec_supercell; // Supercell lattice vectors
    Matrix3 G; // reciprocal lattice vector (2pi*inv(R) )
    Matrix3 GT; // traspose of G
    Matrix3 GGT; // GGT = G*GT
    Matrix3 invGGT; // inverse G

    //========================================================
    // relationship between:
    // ntype, it
    // nat, iat
    // atoms[it].na, ia,
    // atoms[it].nw, iw
    //
    // if know it ==> atoms[it].na; atoms[it].nw
    // if know iat ==> it; ia;
    // if know ia, mush have known it ==> iat
    // if know iwt, must have known it, ia ==> iwt
    //========================================================
    int namax;// the max na among all atom species
    int nwmax;// the max nw among all atom species

    int *iat2it; //iat==>it, distinguish a atom belong to which type
    int *iat2ia; //iat==>ia
	int *iwt2iat; // iwt ==> iat.
	int *iwt2iw; // iwt ==> iw, Peize Lin add 2018-07-02
    IntArray itia2iat;//(it, ia)==>iat, the index in nat, add 2009-3-2 by mohan
    IntArray itiaiw2iwt;//(it, ia, iw)==>iwt, the index in nwfc, add 2009-3-2 by mohan
    //LiuXh add 20180515
    Matrix3 G0;
    Matrix3 GT0;
    Matrix3 GGT0;
    Matrix3 invGGT0;
	
public:
    UnitCell();
    ~UnitCell();
    void print_cell(ofstream &ofs)const;
    void print_cell_xyz(const string &fn)const;
    void print_cell_cif(const string &fn)const;
    const double& getNelec(void)const {return electrons_number;}

    void update_pos_tau(const double* pos);
    void update_pos_taud(const Vector3<double>* posd_in);
    void periodic_boundary_adjustment();
    void bcast_atoms_tau();
    void save_cartesian_position(double* pos)const;

protected:

    double electrons_number;

    double *atom_mass;
    string *atom_label;
    string *pseudo_fn;

#ifdef __MPI
    void bcast_unitcell(void);
    void bcast_unitcell2(void);
#endif

	void set_iat2it(void);

};

#endif //unitcell class
