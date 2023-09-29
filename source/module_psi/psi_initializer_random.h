#ifndef PSI_INITIALIZER_RANDOM_H
#define PSI_INITIALIZER_RANDOM_H

#include "psi_initializer.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
/*
Psi (planewave based wavefunction) initializer: random method
*/
class psi_initializer_random : public psi_initializer
{
    public:
        #ifdef __MPI
        /// @brief parameterized constructor of psi initializer (with MPI support)
        /// @param sf_in interface, link with Structure_Factor ESolver_FP::sf
        /// @param pw_wfc_in interface, link with ModulePW::PW_Basis_K* ESolver_FP::pw_wfc
        /// @param p_ucell_in interface, link with UnitCell GlobalC::ucell
        /// @param p_parakpts_in interface, link with Parallel_Kpoints GlobalC::Pkpoints
        /// @param random_seed_in random seed
        psi_initializer_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, Parallel_Kpoints* p_parakpts_in, int random_seed_in = 1);
        #else
        /// @brief parameterized constructor of psi initializer (without MPI support)
        /// @param sf_in interface, link with Structure_Factor ESolver_FP::sf
        /// @param pw_wfc_in interface, link with ModulePW::PW_Basis_K* ESolver_FP::pw_wfc
        /// @param p_ucell_in interface, link with UnitCell GlobalC::ucell
        /// @param random_seed_in random seed
        psi_initializer_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, int random_seed_in = 1);
        #endif
        /// @brief default destructor
        ~psi_initializer_random();

        // methods
        
        /// @brief kernel function, generating psi value in certain band range
        /// @param psi psi data carrier
        /// @param iw_start the starting band index to fill random value
        /// @param iw_end the end band index
        /// @param ik kpoint index
        void random(std::complex<double>* psi,
                    const int iw_start,
                    const int iw_end,
                    const int ik);
        /// @brief calculate and output planewave wavefunction
        /// @param ik kpoint index
        /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
        psi::Psi<std::complex<double>>* cal_psig(int ik) override;
        /// @brief for variables can be only initialized for once.
        /// @param p_pspot_nl_in (for atomic) interfaces to pseudopot_cell_vnl object, in GlobalC now
        /// @attention if one variable is necessary for all methods, initialize it in constructor, not here.
        void initialize_only_once(pseudopot_cell_vnl* p_pspot_nl_in = nullptr) override {};
};
#endif