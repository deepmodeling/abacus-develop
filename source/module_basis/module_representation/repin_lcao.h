#ifndef RepIn_LCAO_H
#define RepIn_LCAO_H

/*
    This class is for converting the actual psi in ABACUS to the so called ctot
*/

#include "module_basis/module_representation/repin.h"
template<typename T, typename Device>
class RepIn_LCAO : public RepIn<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        RepIn_LCAO() = delete;
        #ifdef __MPI
        /// @brief constructor of RepIn_LCAO
        /// @param sf_in link to Structure_Factor
        /// @param pw_wfc_in link to ModulePW::PW_Basis_K
        /// @param p_ucell_in link to UnitCell
        /// @param p_parakpts_in link to Parallel_Kpoints
        /// @param p_pspot_nl_in link to pseudopot_cell_vnl
        RepIn_LCAO(Structure_Factor* sf_in, 
                   ModulePW::PW_Basis_K* pw_wfc_in, 
                   UnitCell* p_ucell_in, 
                   Parallel_Kpoints* p_parakpts_in, 
                   pseudopot_cell_vnl* p_pspot_nl_in);
        #else
        /// @brief constructor of RepIn_LCAO
        /// @param sf_in link to Structure_Factor
        /// @param pw_wfc_in link to ModulePW::PW_Basis_K
        /// @param p_ucell_in link to UnitCell
        /// @param p_pspot_nl_in link to pseudopot_cell_vnl
        RepIn_LCAO(Structure_Factor* sf_in, 
                   ModulePW::PW_Basis_K* pw_wfc_in, 
                   UnitCell* p_ucell_in, 
                   pseudopot_cell_vnl* p_pspot_nl_in);
        #endif
        ~RepIn_LCAO();
        /*
            central function
        */
        /// @brief calculate pw representation of basis function of input wavefunction
        /// @param psig basis function in pw representation
        void cal_psig(psi::Psi<T, Device>* psig) override;
        /*
            subroutines
        */
        

};

#endif // RepIn_LCAO_H