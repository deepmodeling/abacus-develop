/**
 * @file read_stru.h
 * @brief Functions for reading STRU file.
 */
#ifndef READ_STRU_H
#define READ_STRU_H

#include "atom_spec.h"
#include "source_cell/unitcell.h"

namespace unitcell
{
    /**
     * @brief Check atom positions (tau).
     *
     * @param atoms atom pointer [in]
     * @param ntype number of atom types [in]
     * @param lat0 lattice constant [in]
     * @return true if check passes
     */
    bool check_tau(const Atom* atoms,
                   const int& ntype,
                   const double& lat0);

    /**
     * @brief Check atom displacements (dtau).
     *
     * @param atoms atom pointer [in/out]
     * @param ntype number of atom types [in]
     * @param lat0 lattice constant [in]
     * @param latvec lattice vectors [in]
     */
    void check_dtau(Atom* atoms,
                    const int& ntype,
                    const double& lat0,
                    ModuleBase::Matrix3& latvec);
    
    /**
     * @brief Read atom species information.
     *
     * @param ifa input file stream [in]
     * @param ofs_running output file stream [in]
     * @param ucell unit cell [in/out]
     * @param basis_type basis type [in]
     * @param orbital_dir orbital directory [in]
     * @param init_wfc initial wavefunction type [in]
     * @param onsite_radius onsite radius [in]
     * @param deepks_setorb whether to set orb for deepks [in]
     * @param rpa RPA flag [in]
     * @return true if reading succeeds
     */
    bool read_atom_species(std::ifstream& ifa,
                          std::ofstream& ofs_running,
                          UnitCell& ucell,
                          const std::string& basis_type,
                          const std::string& orbital_dir,
                          const std::string& init_wfc,
                          const double onsite_radius,
                          const bool deepks_setorb,
                          const bool rpa); 
    
    /**
     * @brief Read lattice constant.
     *
     * @param ifa input file stream [in]
     * @param ofs_running output file stream [in]
     * @param lat lattice object [in/out]
     * @return true if reading succeeds
     */
    bool read_lattice_constant(std::ifstream& ifa,
                               std::ofstream& ofs_running,
                               Lattice& lat);
                               
    /**
     * @brief Read atomic positions.
     *
     * @return true (1): no problem.
     * @return false (0): some problems.
     *
     * @param ucell unit cell [in/out]
     * @param ifpos input file stream [in]
     * @param ofs_running output file stream [in]
     * @param ofs_warning warning output file stream [in]
     * @param nspin number of spin components [in]
     * @param basis_type basis type [in]
     * @param orbital_dir orbital directory [in]
     * @param init_wfc initial wavefunction type [in]
     * @param onsite_radius onsite radius [in]
     * @param fixed_atoms whether atoms are fixed [in]
     * @param noncolin non-collinear flag [in]
     * @param calculation calculation type [in]
     * @param esolver_type solver type [in]
     */
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
                            const std::string& esolver_type,
                            const int symmetry);
}

#endif // READ_STRU_H