#include <vector>
#include <fstream>
#include <string>

namespace ModuleIO
{
    /**
     * @brief static version of read_abacus_orb. A delete-new operation may cause the memory leak, 
     * it is better to use std::vector to replace the raw pointer.
     * 
     * @param ifs [in] ifstream from the orbital file, via `std::ifstream ifs(forb);`
     * @param elem [out] element symbol
     * @param ecut [out] planewave energy cutoff
     * @param nr [out] number of radial grid points
     * @param dr [out] radial grid spacing
     * @param nzeta [out] number of zeta functions for each angular momentum
     * @param radials [out] radial orbitals
     * @param rank [in] MPI rank
     */
    void read_abacus_orb(std::ifstream& ifs,
                         std::string& elem,
                         double& ecut,
                         int& nr,
                         double& dr,
                         std::vector<int>& nzeta,
                         std::vector<std::vector<double>>& radials,
                         const int rank = 0);
    
    void write_abacus_orb(std::ofstream& ofs,
                          const std::string& elem,
                          const double& ecut,
                          const int nr,
                          const double dr,
                          const std::vector<int>& nzeta,
                          const std::vector<std::vector<double>>& radials,
                          const int rank = 0);
}