#ifndef CHGMIXING_H
#define CHGMIXING_H

#include "source_estate/elecstate.h" // use pelec
#include "source_estate/module_charge/charge.h" // use chr
#include "source_estate/module_charge/charge_mixing.h" // use p_chgmix
#include "source_io/module_parameter/input_parameter.h" // use Input_para
#include "source_cell/unitcell.h"
#include "source_lcao/module_dftu/dftu.h" // mohan add 2025-11-06

namespace module_charge
{

void chgmixing_ks(const int iter, // scf iteration number
		UnitCell& ucell,
        elecstate::ElecState* pelec, 
        Charge &chr, // charge density
        Charge_Mixing* p_chgmix, // charge mixing class
        const int nrxx, // charge density
        double &drho, // charge density deviation
        bool &oscillate_esolver, // whether the esolver has oscillation of charge density
        bool &conv_esolver,
        const double &hsolver_error,
        const double &scf_thr,
        const double &scf_ene_thr,
        const bool converged_u, // mohan add 2025-11-06
		const Input_para& inp); // input parameters

void chgmixing_ks_pw(const int iter,
        Charge_Mixing* p_chgmix,
        Plus_U &dftu, // mohan add 2025-11-06
		const Input_para& inp); // input parameters

void chgmixing_ks_lcao(const int iter, // scf iteration number
        Charge_Mixing* p_chgmix, // charge mixing class
        Plus_U &dftu, // mohan add 2025-11-06
        const int nnr, // dimension of density matrix
		const Input_para& inp); // input parameters

}


#endif
