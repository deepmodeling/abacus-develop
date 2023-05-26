#pragma once
#include "Exx_LRI.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include <memory>

template<typename Tdata>
class Exx_LRI_Interface
{
public:
    /// @brief  Constructor for Exx_LRI_Interface
    /// @param exx_lri
    Exx_LRI_Interface(Exx_LRI<Tdata>& exx_lri) : exx_lri(&exx_lri) {}
    Exx_LRI_Interface() = delete;

    // Processes in ESolver_KS_LCAO
    /// @brief in beforescf: set xc type, opt_orb, do DM mixing
    void exx_beforescf(const K_Vectors& kv);

        /// @brief in eachiterinit:  do DM mixing and calculate Hexx when entering 2nd SCF
    void exx_eachiterinit(const Local_Orbital_Charge& loc, const Charge_Mixing& chgmix, const int& iter);

    /// @brief in hamilt2density: calcate Hexx and Eexx
    void exx_hamilt2density(elecstate::ElecState& elec, const Parallel_Orbitals& pv);

    /// @brief: in do_after_converge: add exx operators; do DM mixing if seperate loop
    bool exx_after_converge(
        hamilt::Hamilt<double>& hamilt,
        LCAO_Matrix& lm,
        const Local_Orbital_Charge& loc,
        const K_Vectors& kv,
        int& iter);
    
private:
    Exx_LRI<Tdata>* exx_lri;
};
#include "Exx_LRI_interface.hpp"