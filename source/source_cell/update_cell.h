/**
 * @file update_cell.h
 * @brief Functions for updating cell information.
 *
 * This file contains functions for:
 * 1. remake_cell: for constrained vc-relaxation where type of lattice
 *    is fixed, adjust the lattice vectors
 * 2. setup_cell_after_vc: setup cell after vc-relaxation
 * 3. periodic_boundary_adjustment: adjust the boundary of the cell
 * 4. update_pos_tau: update the Cartesian coordinate position of the atoms
 */
#ifndef UPDATE_CELL_H
#define UPDATE_CELL_H

#include "unitcell_data.h"
#include "unitcell.h"  

namespace unitcell
{
    /**
     * @brief Adjust lattice vectors for constrained vc-relaxation.
     *
     * For constrained vc-relaxation where type of lattice is fixed,
     * adjust the lattice vectors.
     *
     * @param lat lattice object [in/out]
     */
    void remake_cell(Lattice& lat);

    /**
     * @brief Setup cell after vc-relaxation.
     *
     * @param ucell unit cell [in/out]
     * @param log output file stream [in]
     * @param nspin number of spin components [in]
     */
    void setup_cell_after_vc(UnitCell& ucell, std::ofstream& log, const int nspin);
    
    /**
     * @brief Check the boundary of the cell.
     *
     * For each atom, the taud in three directions should be in the range of [-1,1).
     *
     * @param atoms the atoms to be adjusted [in/out]
     * @param latvec the lattice of the atoms [in]
     * @param ntype the number of types of the atoms [in]
    */
    void periodic_boundary_adjustment(Atom* atoms,
                                      const ModuleBase::Matrix3& latvec,
                                      const int ntype);

    /** 
    * @brief Update the position and tau of the atoms.
    * 
    * @param lat the lattice of the atoms [in]
    * @param pos the position of the atoms [in]
    * @param ntype the number of types of the atoms [in]
    * @param nat the number of atoms [in]
    * @param atoms the atoms to be updated [out]
    */
    void update_pos_tau(const Lattice& lat,
                        const double* pos,
                        const int ntype,
                        const int nat,
                        Atom* atoms);
    
    /**
     * @brief Update the position and taud of the atoms.
     * 
     * @param lat the lattice of the atoms [in]
     * @param posd_in the position of the atoms in direct coordinate system [in]
     * @param ntype the number of types of the atoms [in]
     * @param nat the number of atoms [in]
     * @param atoms the atoms to be updated [out]
     */
    void update_pos_taud(const Lattice& lat,
                         const double* posd_in,
                         const int ntype,
                         const int nat,
                         Atom* atoms);

    /**
     * @brief Update the position and taud of the atoms (Vector3 version).
     * 
     * @param lat the lattice of the atoms [in]
     * @param posd_in the position of the atoms in direct coordinate system [in]
     * @param ntype the number of types of the atoms [in]
     * @param nat the number of atoms [in]
     * @param atoms the atoms to be updated [out]
     */
    void update_pos_taud(const Lattice& lat,
                         const ModuleBase::Vector3<double>* posd_in,
                         const int ntype,
                         const int nat,
                         Atom* atoms);

    /**
     * @brief Update the velocity of the atoms.
     * 
     * @param vel_in the velocity of the atoms [in]
     * @param ntype the number of types of the atoms [in]
     * @param nat the number of atoms [in]
     * @param atoms the atoms to be updated [out]
    */
    void update_vel(const ModuleBase::Vector3<double>* vel_in,
                    const int ntype,
                    const int nat,
                    Atom* atoms);
}

#endif // UPDATE_CELL_H