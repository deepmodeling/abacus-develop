/**
 * @file print_cell.h
 * @brief Functions for printing cell information.
 */
#ifndef PRINT_CELL_H
#define PRINT_CELL_H

#include "atom_spec.h"
#include "source_cell/unitcell.h"

namespace unitcell
{
    /**
     * @brief Print atom positions (tau).
     *
     * @param atoms atom pointer [in]
     * @param Coordinate coordinate system type [in]
     * @param ntype number of atom types [in]
     * @param lat0 lattice constant [in]
     * @param ofs output file stream [in]
     */
    void print_tau(Atom* atoms,
                   const std::string& Coordinate,
                   const int ntype,
                   const double lat0,
                   std::ofstream &ofs);
    
    /**
     * @brief Print STRU file according to given settings.
     *
     * @note UnitCell class is too heavy, this function would be moved elsewhere.
     *
     * @param ucell reference of unitcell [in]
     * @param atoms Atom list [in]
     * @param latvec lattice parameter vector [in]
     * @param fn STRU file name [in]
     * @param nspin number of spin channels [in]
     * @param direct true for direct coords, false for cartesian coords [in]
     * @param vel true for printing velocities [in]
     * @param magmom true for printing Mulliken population analysis produced magmom [in]
     * @param orb true for printing NUMERICAL_ORBITAL section [in]
     * @param dpks_desc true for printing NUMERICAL_DESCRIPTOR section [in]
     * @param iproc GlobalV::MY_RANK [in]
     */
    void print_stru_file(const UnitCell& ucell,
                         const Atom*     atoms,
                         const ModuleBase::Matrix3& latvec,
                         const std::string& fn,
                         const int& nspin = 1,
                         const bool& direct = false,
                         const bool& vel = false,
                         const bool& magmom = false,
                         const bool& orb = false,
                         const bool& dpks_desc = false,
                         const int& iproc = 0);

    /**
     * @brief Print basic unitcell information to output stream.
     *
     * @param ucell reference of unitcell [in]
     * @param ofs output file stream [in]
     */
    void print_cell(const UnitCell& ucell, std::ofstream& ofs);
}

#endif
