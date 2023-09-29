#ifndef PSI_INITIALIZER_NAO_H
#define PSI_INITIALIZER_NAO_H
#include "psi_initializer.h"
#include "module_base/realarray.h"

/*
Psi (planewave based wavefunction) initializer: numerical atomic orbital method
*/
class psi_initializer_nao : public psi_initializer
{
    public:
        #ifdef __MPI
        /// @brief parameterized constructor of psi initializer (with MPI support)
        /// @param sf_in interface, link with Structure_Factor ESolver_FP::sf
        /// @param pw_wfc_in interface, link with ModulePW::PW_Basis_K* ESolver_FP::pw_wfc
        /// @param p_ucell_in interface, link with UnitCell GlobalC::ucell
        /// @param p_parakpts_in interface, link with Parallel_Kpoints GlobalC::Pkpoints
        /// @param random_seed_in random seed
        psi_initializer_nao(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, Parallel_Kpoints* p_parakpts_in, int random_seed_in = 1);
        #else
        /// @brief parameterized constructor of psi initializer (without MPI support)
        /// @param sf_in interface, link with Structure_Factor ESolver_FP::sf
        /// @param pw_wfc_in interface, link with ModulePW::PW_Basis_K* ESolver_FP::pw_wfc
        /// @param p_ucell_in interface, link with UnitCell GlobalC::ucell
        /// @param random_seed_in random seed
        psi_initializer_nao(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, int random_seed_in = 1);
        #endif
        /// @brief default destructor
        ~psi_initializer_nao();

        // methods

        /// @brief calculate and output planewave wavefunction
        /// @param ik kpoint index
        /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
        psi::Psi<std::complex<double>>* cal_psig(int ik) override;

        /// @brief initialize only once, for nao, it should be, read numerical orbitals, create ovlp_Xjlq(, calculate ovlp_flzjlq)
        /// @param p_pspot_nl_in (for atomic) interfaces to pseudopot_cell_vnl object, in GlobalC now
        /// @attention if one variable is necessary for all methods, initialize it in constructor, not here.
        void initialize_only_once(pseudopot_cell_vnl* p_pspot_nl_in = nullptr) override;
        // setters

        /// @brief setter of numerical orbital files
        /// @param orbital_files array storing numerical orbital files
        void set_orbital_files(std::string* orbital_files);
        // I wont write a function to set ovlp_flzjlq, it is totally useless
        
        /// @brief allocate memory for ovlp_flzjlq and initialize all elements to 0
        /// @attention warning! p_ucell must be set in advance!
        void create_ovlp_Xjlq() override;
        /// @brief before refactor and reorganization of UnitCell class, it is temporary to write this function here.
        /// In future version, it will be moved into UnitCell class.
        void read_orbital_files();
        /// @brief calculate overlap integral between f_{l\\zeta} the radial numerical orbital and spherical Bessel function
        void cal_ovlp_flzjlq();
        
        // getters

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
#endif