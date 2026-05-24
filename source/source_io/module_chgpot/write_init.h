// =====================================================================
// This module handles the output of initial charge density and potential
// to cube files in real space. It is part of the module_io package.
//
// Output files are named according to the following convention:
//   - Initial charge density (out_chg = 2): chgg{#}_ini.cube or chgs{#}g{#}_ini.cube
//   - Initial potential (out_pot = 3): potg{#}_ini.cube or pots{#}g{#}_ini.cube
//   - Geometry step index starts from 1 (geom_step = istep + 1)
//
// Usage:
//   ModuleIO::write_chg_init(ucell, para_grid, chr, efermi, istep, out_dir, inp);
//   ModuleIO::write_pot_init(ucell, para_grid, pelec, istep, out_dir, inp);
//
// Module: module_io/module_chgpot
// =====================================================================

#ifndef WRITE_INIT_H
#define WRITE_INIT_H

#include <string>
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
    const std::string& out_dir,
    const Input_para& inp);

// Write initial effective potential to cube file in real space.
// Triggered when inp.out_pot[0] == 3.
// Output frequency is controlled by out_freq_ion (output at step 0 or every out_freq_ion steps).
// Output file naming convention:
//   nspin=1: potg{geom_step}_ini.cube (e.g., potg1_ini.cube)
//   nspin=2/4: pots{spin}g{geom_step}_ini.cube (e.g., pots1g1_ini.cube, pots2g1_ini.cube)
// Note: geom_step starts from 1 (geom_step = istep + 1).
void write_pot_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    elecstate::ElecState *pelec,
    const int istep,
    const std::string& out_dir,
    const Input_para& inp);

}

#endif
