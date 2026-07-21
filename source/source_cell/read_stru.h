#ifndef READ_STRU_H
#define READ_STRU_H

#include "atom_spec.h"
#include "source_cell/unitcell.h"
namespace unitcell
{
    bool check_tau(const Atom* atoms,
                   const int& ntype,
                   const double& lat0);
    void check_dtau(Atom* atoms,
                    const int& ntype,
                    const double& lat0,
                    ModuleBase::Matrix3& latvec);
    
    bool read_atom_species(std::ifstream& ifa,
                          std::ofstream& ofs_running,
                          UnitCell& ucell,
                          const std::string& basis_type,
                          const std::string& orbital_dir,
                          const std::string& init_wfc,
                          const double onsite_radius,
                          const bool deepks_setorb,
                          const bool rpa); 
    
    bool read_lattice_constant(std::ifstream& ifa,
                               std::ofstream& ofs_running,
                               Lattice& lat);
                               
    // Read atomic positions
    // return 1: no problem.
    // return 0: some problems.
    bool read_atom_positions(UnitCell& ucell,
                            std::ifstream &ifpos, 
                            std::ofstream &ofs_running, 
                            std::ofstream &ofs_warning,
                            const int nspin,
                            const std::string& basis_type,
                            const std::string& orbital_dir,
                            const std::string& init_wfc,
                            const double onsite_radius,
                            const bool fixed_atoms,
                            const bool noncolin,
                            const std::string& calculation,
                            const std::string& esolver_type);
}
#endif // READ_STRU_H