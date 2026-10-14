// =====================================================================
// This module handles the output of initial charge density and potential
// to cube files in real space. It is part of the module_io package.
//
// Output files are named according to the following convention:
//   - Frequency is out_freq_td for RT-TDDFT and out_freq_ion otherwise.
//   - Frequency = 0: chg_ini.cube/pot_ini.cube, chgs1_ini.cube/chgs2_ini.cube
//   - Frequency > 0: chgg{#}_ini.cube/potg{#}_ini.cube, chgs{#}g{#}_ini.cube
//   - File indices start from 1 (istep + 1); RT-TDDFT uses electronic steps.
//
// Usage:
//   ModuleIO::write_chg_init(ucell, para_grid, chr, efermi, istep, out_dir, inp);
//   ModuleIO::write_pot_init(ucell, para_grid, pelec, istep, out_dir, inp);
//
// Module: module_io/module_chgpot
// =====================================================================

#ifndef WRITE_INIT_H
#define WRITE_INIT_H

#include "source_estate/elecstate.h"
#include "source_estate/fp_energy.h"
#include "source_estate/module_charge/charge.h"
#include "source_io/module_parameter/input_parameter.h"

#include <fstream>
#include <string>

namespace ModuleIO
{

// Generate initial data file name.
// prefix: "chg" or "pot"
std::string gen_ini_filename(const std::string& prefix,
                             const std::string& out_dir,
                             const int nspin,
                             const int is,
                             const int istep,
                             const bool include_geom_step);

// Write initial charge density to cube file in real space.
// Triggered when inp.out_chg[0] == 2.
// Use out_freq_td and an electronic istep for RT-TDDFT, out_freq_ion and an ionic istep otherwise.
// Zero frequency overwrites the same file; positive frequency selects istep % frequency == 0.
// File naming and one-based indices follow the convention above.
void write_chg_init(const UnitCell& ucell,
                    const Parallel_Grid& para_grid,
                    const Charge& chr,
                    const elecstate::Efermi& efermi,
                    const int istep,
                    const std::string& out_dir,
                    const Input_para& inp,
                    const bool two_fermi,
                    std::ofstream& ofs_running);

// Write initial effective potential to cube file in real space.
// Triggered when inp.out_pot[0] == 3.
// Use out_freq_td and an electronic istep for RT-TDDFT, out_freq_ion and an ionic istep otherwise.
// Zero frequency overwrites the same file; positive frequency selects istep % frequency == 0.
// File naming and one-based indices follow the convention above.
void write_pot_init(const UnitCell& ucell,
                    const Parallel_Grid& para_grid,
                    elecstate::ElecState* pelec,
                    const int istep,
                    const std::string& out_dir,
                    const Input_para& inp,
                    const bool two_fermi,
                    std::ofstream& ofs_running);

} // namespace ModuleIO

#endif
