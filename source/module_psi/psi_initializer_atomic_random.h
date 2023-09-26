#ifndef PSI_INITIALIZER_ATOMIC_RANDOM_H
#define PSI_INITIALIZER_ATOMIC_RANDOM_H
#include "psi_initializer_atomic.h"

/*
Psi (planewave based wavefunction) initializer: atomic+random
*/
class psi_initializer_atomic_random : public psi_initializer_atomic
{
    public:

        /// @brief constructor of psi initializer (atomic+random)
        /// @param sf_in Structure factor interface, link ESolver
        /// @param pw_wfc_in ModulePW::PW_Basis_K interface, link ESolver
        psi_initializer_atomic_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        
        /// @brief default destructor
        ~psi_initializer_atomic_random();

        // methods
        /// @brief calculate and output planewave wavefunction
        /// @param ik kpoint index
        /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
        psi::Psi<std::complex<double>>* cal_psig(int ik) override;
    private:
};
#endif