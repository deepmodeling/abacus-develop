#ifndef SETUP_POT_H
#define SETUP_POT_H

#include "source_base/module_device/device.h" // use Device
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_estate/elecstate.h"
#include "source_pw/module_pwdft/VL_in_pw.h"
#include "source_hamilt/hamilt.h"

namespace pw
{

template <typename T, typename Device>
void setup_pot(const int istep, 
		UnitCell& ucell, // unitcell 
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<T, Device>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<T, Device>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
		const Input_para& inp); // input parameters

}



#endif
