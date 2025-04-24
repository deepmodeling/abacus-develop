#ifndef PRINT_CELL_H
#define PRINT_CELL_H

#include "atom_spec.h"
#include "module_cell/unitcell.h"
namespace unitcell
{
    void print_tau(Atom* atoms,
                   const std::string& Coordinate,
                   const int ntype,
                   const double lat0,
                   std::ofstream &ofs);
    
      /**
     * @brief UnitCell class is too heavy, this function would be moved
     * elsewhere. Print STRU file respect to given setting
     *
     * @param ucell reference of unitcell
     * @param atoms Atom list 
     * @param latvec lattice const parmater vector 
     * @param fn STRU file name
     * @param nspin PARAM.inp.nspin feed in
     * @param direct true for direct coords, false for cartesian coords
     * @param vol true for printing velocities
     * @param magmom true for printing Mulliken population analysis produced
     * magmom
     * @param orb true for printing NUMERICAL_ORBITAL section
     * @param dpks_desc true for printing NUMERICAL_DESCRIPTOR section
     * @param iproc GlobalV::MY_RANK feed in
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
}

#endif
