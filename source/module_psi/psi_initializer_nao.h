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
        /// @brief constructor of psi initializer (nao)
        /// @param sf_in Structure factor interface, link ESolver
        /// @param pw_wfc_in ModulePW::PW_Basis_K interface, link ESolver
        psi_initializer_nao(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        /// @brief default destructor
        ~psi_initializer_nao();

        // methods

        /// @brief calculate and output planewave wavefunction
        /// @param ik kpoint index
        /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
        psi::Psi<std::complex<double>>* cal_psig(int ik) override;

        // setters

        /// @brief setter of numerical orbital files
        /// @param orbital_files array storing numerical orbital files
        void set_orbital_files(std::string* orbital_files);
        // I wont write a function to set ovlp_flzjlq, it is totally useless
        
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
};
#endif