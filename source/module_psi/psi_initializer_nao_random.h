#ifndef PSI_INITIALIZER_NAO_RANDOM_H
#define PSI_INITIALIZER_NAO_RANDOM_H
#include "psi_initializer_nao.h"

/*
Psi (planewave based wavefunction) initializer: numerical atomic orbital + random method
*/
class psi_initializer_nao_random : public psi_initializer_nao
{
    public:
        
        /// @brief constructor of psi initializer (nao+random)
        /// @param sf_in Structure factor interface, link ESolver
        /// @param pw_wfc_in ModulePW::PW_Basis_K interface, link ESolver
        psi_initializer_nao_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        /// @brief default destructor
        ~psi_initializer_nao_random();

        /// @brief calculate and output planewave wavefunction
        /// @param ik kpoint index
        /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
        psi::Psi<std::complex<double>>* cal_psig(int ik) override;
    private:
};
#endif