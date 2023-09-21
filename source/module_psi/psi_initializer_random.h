#ifndef PSI_INITIALIZER_RANDOM_H
#define PSI_INITIALIZER_RANDOM_H

#include "psi_initializer.h"
/*
Psi (planewave based wavefunction) initializer: random method
*/
class psi_initializer_random : public psi_initializer
{
    public:
        /// @brief constructor of psi initializer (random)
        /// @param sf_in Structure factor interface, link ESolver
        /// @param pw_wfc_in ModulePW::PW_Basis_K interface, link ESolver
        psi_initializer_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        /// @brief default destructor
        ~psi_initializer_random();

        // methods
        
        /// @brief kernel function, generating psi value in certain band range
        /// @param psi psi data carrier
        /// @param iw_start the starting band index to fill random value
        /// @param iw_end the end band index
        /// @param ik kpoint index
        /// @param wfc_basis ModulePW::PW_Basis_K interface, will be automatically linked
        void random(std::complex<double>* psi,
                    const int iw_start,
                    const int iw_end,
                    const int ik,
                    const ModulePW::PW_Basis_K* wfc_basis);
        /// @brief calculate and output planewave wavefunction
        /// @param ik kpoint index
        /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
        psi::Psi<std::complex<double>>* cal_psig(int ik) override;
};
#endif