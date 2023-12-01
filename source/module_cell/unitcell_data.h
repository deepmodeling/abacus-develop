#include "module_base/matrix3.h"
#include "module_base/intarray.h"
/// @brief info of lattice 
struct Lattice
{
    std::string Coordinate; // "Direct" or "Cartesian" or "Cartesian_angstrom"
    std::string latName; // Lattice name
    double lat0; // Lattice constant(bohr)(a.u.)
    double lat0_angstrom;// Lattice constant(angstrom)
    double tpiba;// 2*pi / lat0;
    double tpiba2; // tpiba ^ 2
    double omega;// the volume of the unit cell
    int* lc;  // Change the lattice vectors or not

    ModuleBase::Matrix3 latvec; // Unitcell lattice vectors
    ModuleBase::Vector3<double> a1, a2, a3; // Same as latvec, just at another form.
    ModuleBase::Vector3<double> latcenter; // (a1+a2+a3)/2 the center of vector
    ModuleBase::Matrix3 latvec_supercell; // Supercell lattice vectors
    ModuleBase::Matrix3 G; // reciprocal lattice vector (2pi*inv(R) )
    ModuleBase::Matrix3 GT; // traspose of G
    ModuleBase::Matrix3 GGT; // GGT = G*GT
    ModuleBase::Matrix3 invGGT; // inverse G
};

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
/// @brief usefull data and index maps
struct Statistics
{
    int ntype;// number of atom species in UnitCell
    int nat; // total number of atoms of all species in unitcell
    int* iat2it; //iat==>it, distinguish a atom belong to which type
    int* iat2ia; //iat==>ia
    int* iwt2iat; // iwt ==> iat.
    int* iwt2iw; // iwt ==> iw, Peize Lin add 2018-07-02
    ModuleBase::IntArray itia2iat;//(it, ia)==>iat, the index in nat, add 2009-3-2 by mohan
    int namax;// the max na among all atom species
    int nwmax;// the max nw among all atom species
};