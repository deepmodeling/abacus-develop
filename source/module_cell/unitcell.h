#ifndef UNITCELL_H
#define UNITCELL_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix3.h"
#include "../module_base/intarray.h"
#include "../src_io/output.h"
#ifndef __CMD
#include "../src_pw/magnetism.h"
#endif

#include "atom_spec.h"

// the "base class" of UnitCell_pseudo.
// provide the basic information about unitcell.
class UnitCell
{
public:
    Atom *atoms;

#ifndef __CMD
    Magnetism magnet;  // magnetism Yu Liu 2021-07-03
    void cal_ux();
#endif
    bool judge_parallel(double a[3],ModuleBase::Vector3<double> b);
	double *atom_mag;
	int n_mag_at;

    int ntype;// number of atom species in UnitCell
    int nat; // total number of atoms of all species in unitcell
    std::string Coordinate; // "Direct" or "Cartesian" or "Cartesian_angstrom"
    std::string latName; // Lattice name
    double lat0; // Lattice constant(bohr)(a.u.)
    double lat0_angstrom;// Lattice constant(angstrom)
    double tpiba;// 2*pi / lat0;
    double tpiba2; // tpiba ^ 2
    double omega;// the volume of the unit cell

    ModuleBase::Matrix3 latvec; // Unitcell lattice vectors
	int *lc;  // Change the lattice vectors or not
	ModuleBase::Vector3<double> a1,a2,a3; // Same as latvec, just at another form.
	ModuleBase::Vector3<double> latcenter; // (a1+a2+a3)/2 the center of vector
    ModuleBase::Matrix3 latvec_supercell; // Supercell lattice vectors
    ModuleBase::Matrix3 G; // reciprocal lattice vector (2pi*inv(R) )
    ModuleBase::Matrix3 GT; // traspose of G
    ModuleBase::Matrix3 GGT; // GGT = G*GT
    ModuleBase::Matrix3 invGGT; // inverse G

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
    ModuleBase::IntArray itia2iat;//(it, ia)==>iat, the index in nat, add 2009-3-2 by mohan
    ModuleBase::IntArray itiaiw2iwt;//(it, ia, iw)==>iwt, the index in nwfc, add 2009-3-2 by mohan
    //LiuXh add 20180515
    ModuleBase::Matrix3 G0;
    ModuleBase::Matrix3 GT0;
    ModuleBase::Matrix3 GGT0;
    ModuleBase::Matrix3 invGGT0;
	
    //I'm doing a bad thing here! Will change later
    bool ionic_position_updated = false; //whether the ionic position has been updated
    bool cell_parameter_updated = false; //whether the cell parameters are updated

private:
    ModuleBase::Matrix3 stress; //calculate stress on the cell

public:
    UnitCell();
    ~UnitCell();
    void print_cell(std::ofstream &ofs)const;
    void print_cell_xyz(const std::string &fn)const;
    void print_cell_cif(const std::string &fn)const;

    void update_pos_tau(const double* pos);
    void update_pos_tau(const ModuleBase::Vector3<double>* posd_in);
    void update_pos_taud(const ModuleBase::Vector3<double>* posd_in);
    void update_vel(const ModuleBase::Vector3<double>* vel_in);
    void periodic_boundary_adjustment();
    void bcast_atoms_tau();
    void save_cartesian_position(double* pos)const;
    void save_cartesian_position(ModuleBase::Vector3<double>* pos)const;
    void save_cartesian_position_original(double* pos)const;
    void save_cartesian_position_original(ModuleBase::Vector3<double>* pos)const;
    bool judge_big_cell(void)const;

    void update_stress(ModuleBase::matrix &scs); //updates stress
    void update_force(ModuleBase::matrix &fcs); //updates force in Atom

    double *atom_mass;
    std::string *atom_label;
    std::string *pseudo_fn;

#ifdef __MPI
    void bcast_unitcell(void);
    void bcast_unitcell2(void);
#endif

	void set_iat2itia(void);

};

#endif //unitcell class
