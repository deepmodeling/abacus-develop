#ifndef READ_STRU_H
#define READ_STRU_H

#include "atom_spec.h"
#include "module_cell/unitcell.h"
namespace unitcell
{
    bool check_tau(const Atom* atoms,
                   const int& ntype,
                   const double& lat0);
    void check_dtau(Atom* atoms,
                    const int& ntype,
                    const double& lat0,
                    ModuleBase::Matrix3& latvec);
    
    // read in the atom information for each type of atom
    bool read_atom_species(std::ifstream& ifa,
                          std::ofstream& ofs_running,
                          UnitCell& ucell); 
    
    bool read_lattice_constant(std::ifstream& ifa,
                               std::ofstream& ofs_running,
                               Lattice& lat);
                               
    // Read atomic positions
    // return 1: no problem.
    // return 0: some problems.
    bool read_atom_positions(UnitCell& ucell,
                            std::ifstream &ifpos, 
                            std::ofstream &ofs_running, 
                            std::ofstream &ofs_warning);
}
#endif // READ_STRU_H