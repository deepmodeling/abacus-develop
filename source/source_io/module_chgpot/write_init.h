// =====================================================================
// This module handles the output of initial charge density and potential
// to cube files in real space. It is part of the module_io package.
//
// Output files are named according to the following convention:
//   - out_freq_ion = 0: chg_ini.cube/pot_ini.cube, chgs1_ini.cube/chgs2_ini.cube
//   - out_freq_ion > 0: chgg{#}_ini.cube/potg{#}_ini.cube, chgs{#}g{#}_ini.cube
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

// Generate initial data file name.
// prefix: "chg" or "pot"
std::string gen_ini_filename(
    const std::string& prefix,
    const std::string& out_dir,
    const int nspin,
    const int is,
    const int istep,
    const bool include_geom_step);

// Write initial charge density to cube file in real space.
// Triggered when inp.out_chg[0] == 2.
// Output frequency is controlled by out_freq_ion:
//   out_freq_ion = 0: every step output (overwrite same file)
//   out_freq_ion > 0: output every out_freq_ion steps
// Output file naming convention:
//   out_freq_ion = 0: chg_ini.cube (nspin=1), chgs1_ini.cube/chgs2_ini.cube (nspin=2/4)
//   out_freq_ion > 0: chgg{geom_step}_ini.cube (nspin=1), chgs{spin}g{geom_step}_ini.cube (nspin=2/4)
// Note: geom_step starts from 1 (geom_step = istep + 1).
void write_chg_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    const Charge &chr,
    const elecstate::Efermi &efermi,
    const int istep,
    const std::string& out_dir,
    const Input_para& inp,
    const bool two_fermi);

// Write initial effective potential to cube file in real space.
// Triggered when inp.out_pot[0] == 3.
// Output frequency is controlled by out_freq_ion:
//   out_freq_ion = 0: every step output (overwrite same file)
//   out_freq_ion > 0: output every out_freq_ion steps
// Output file naming convention:
//   out_freq_ion = 0: pot_ini.cube (nspin=1), pots1_ini.cube/pots2_ini.cube (nspin=2/4)
//   out_freq_ion > 0: potg{geom_step}_ini.cube (nspin=1), pots{spin}g{geom_step}_ini.cube (nspin=2/4)
// Note: geom_step starts from 1 (geom_step = istep + 1).
void write_pot_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    elecstate::ElecState *pelec,
    const int istep,
    const std::string& out_dir,
    const Input_para& inp,
    const bool two_fermi);

}

#endif
