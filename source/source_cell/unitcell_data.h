/**
 * @file unitcell_data.h
 * @brief Data structures for unit cell information.
 */
#ifndef UNITCELL_DATA_H
#define UNITCELL_DATA_H

#include <vector>

#include "source_base/intarray.h"
#include "source_base/matrix3.h"

/**
 * @brief Lattice information.
 */
struct Lattice
{
    std::string Coordinate = "Direct"; ///< "Direct" or "Cartesian" or "Cartesian_angstrom"
    std::string latName = "user_defined_lattice";  ///< Lattice name
    double lat0 = 0.0;                 ///< Lattice constant(bohr)(a.u.)
    double lat0_angstrom = 0.0;        ///< Lattice constant(angstrom)
    double tpiba = 0.0;                ///< 2*pi / lat0
    double tpiba2 = 0.0;               ///< tpiba ^ 2
    double omega = 0.0;                ///< the volume of the unit cell
    std::vector<int> lat_axis_free{0, 0, 0};      ///< whether each lattice axis (a,b,c) is allowed to relax (0=fixed, 1=free)

    ModuleBase::Matrix3 latvec = ModuleBase::Matrix3();           ///< Unitcell lattice vectors, in unit of lat0
    ModuleBase::Vector3<double> a1, a2, a3;                       ///< Same as latvec, just at another form
    ModuleBase::Vector3<double> latcenter;                        ///< (a1+a2+a3)/2 the center of vector
    ModuleBase::Matrix3 latvec_supercell = ModuleBase::Matrix3(); ///< Supercell lattice vectors
    ModuleBase::Matrix3 G = ModuleBase::Matrix3();                ///< reciprocal lattice vector (2pi*inv(R))
    ModuleBase::Matrix3 GT = ModuleBase::Matrix3();               ///< transpose of G
    ModuleBase::Matrix3 GGT = ModuleBase::Matrix3();              ///< GGT = G*GT
    ModuleBase::Matrix3 invGGT = ModuleBase::Matrix3();           ///< inverse G
};

#endif
