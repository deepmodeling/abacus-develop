#ifndef RepIn_NAO_H
#define RepIn_NAO_H

#include "module_basis/module_representation/repin.h"

template<typename T, typename Device>
class RepIn_NAO : public RepIn<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        RepIn_NAO() = delete;
        #ifdef __MPI
        /// @brief constructor of RepIn
        /// @param sf_in link to Structure_Factor
        /// @param pw_wfc_in link to ModulePW::PW_Basis_K
        /// @param p_ucell_in link to UnitCell
        /// @param p_parakpts_in link to Parallel_Kpoints
        /// @param p_pspot_nl_in link to pseudopot_cell_vnl
        RepIn_NAO(Structure_Factor* sf_in, 
                  ModulePW::PW_Basis_K* pw_wfc_in, 
                  UnitCell* p_ucell_in, 
                  Parallel_Kpoints* p_parakpts_in, 
                  pseudopot_cell_vnl* p_pspot_nl_in);
        #else
        /// @brief constructor of RepIn
        /// @param sf_in link to Structure_Factor
        /// @param pw_wfc_in link to ModulePW::PW_Basis_K
        /// @param p_ucell_in link to UnitCell
        /// @param p_pspot_nl_in link to pseudopot_cell_vnl
        RepIn_NAO(Structure_Factor* sf_in, 
                  ModulePW::PW_Basis_K* pw_wfc_in, 
                  UnitCell* p_ucell_in, 
                  pseudopot_cell_vnl* p_pspot_nl_in);
        #endif
        ~RepIn_NAO();
        /*
            central function
        */
        /// @brief calculate pw representation of basis function of input wavefunction
        /// @param psig basis function in pw representation
        void cal_psig(psi::Psi<T, Device>* psig) override;
        /*
            initialization and subroutines
        */
        /// @brief call set_orbital_files, create_ovlp_flzjlq, read_orbital_files, cal_ovlp_flzjlq
        /// @param orbital_files array storing numerical orbital files
        void initialize(std::string* orbital_files) override;
        /// @brief setter of numerical orbital files
        /// @param orbital_files array storing numerical orbital files
        void set_orbital_files(std::string* orbital_files);
        /// @brief allocate memory for ovlp_flzjlq and initialize all elements to 0
        /// @attention warning! p_ucell must be set in advance!
        void create_ovlp_flzjlq();
        /// @brief before refactor and reorganization of UnitCell class, it is temporary to write this function here.
        /// In future version, it will be moved into UnitCell class.
        void read_orbital_files();
        /// @brief calculate overlap integral between f_{l\\zeta} the radial numerical orbital and spherical Bessel function
        void cal_ovlp_flzjlq();
        /*
            getters
        */
        /// @brief getter of orbital filenames
        /// @return orbital filenames in array
        std::vector<std::string> get_orbital_files() const { return orbital_files; }
        /// @brief getter of matrix of overlap between numerical orbital and Spherical Bessel function
        /// @return ovlp_flzjlq
        ModuleBase::realArray get_ovlp_flzjlq() const { return ovlp_flzjlq; }
    private:
        std::vector<std::string> orbital_files;
        ModuleBase::realArray ovlp_flzjlq;
        /// @brief number of realspace grids per type per chi, [itype][ichi]
        std::vector<std::vector<int>> n_rgrid;
        /// @brief data of numerical atomic orbital per type per chi per position, [itype][ichi][ir]
        std::vector<std::vector<std::vector<double>>> flz;
        /// @brief r of numerical atomic orbital per type per chi per position, [itype][ichi][ir]
        std::vector<std::vector<std::vector<double>>> rgrid;
};

#endif // RepIn_NAO_H