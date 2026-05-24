#ifndef WRITE_INIT_H
#define WRITE_INIT_H

#include "source_io/module_parameter/input_parameter.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/fp_energy.h"
#include "source_estate/elecstate.h"

namespace ModuleIO
{

// Write initial charge density to cube file in real space.
// Triggered when inp.out_chg[0] == 2.
// Output frequency is controlled by out_freq_ion (output at step 0 or every out_freq_ion steps).
// Output file naming convention:
//   nspin=1: chgg{geom_step}_ini.cube (e.g., chgg1_ini.cube)
//   nspin=2/4: chgs{spin}g{geom_step}_ini.cube (e.g., chgs1g1_ini.cube, chgs2g1_ini.cube)
// Note: geom_step starts from 1 (geom_step = istep + 1).
void write_chg_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    const Charge &chr,
    const elecstate::Efermi &efermi,
    const int istep,
    const Input_para& inp);

// Write initial effective potential to cube file in real space.
// Triggered when inp.out_pot[0] == 3.
// Output: pot_ini.cube (nspin=1) or potsN_ini.cube (nspin=2/4)
void write_pot_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    elecstate::ElecState *pelec,
    const int istep,
    const Input_para& inp);

}

#endif
