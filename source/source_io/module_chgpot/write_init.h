#ifndef WRITE_INIT_H
#define WRITE_INIT_H

#include "source_io/module_parameter/input_parameter.h" // use inp
#include "source_cell/parallel_kpoints.h" // use para_grid
#include "source_estate/module_charge/charge.h" // use chg
#include "source_estate/fp_energy.h" // use efermi
#include "source_estate/elecstate.h" // use pelec

namespace ModuleIO
{

void write_chg_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    const Charge &chr,
    const elecstate::Efermi &efermi,
    const int istep,
    const Input_para& inp);

void write_pot_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    elecstate::ElecState *pelec,
    const int istep,
    const Input_para& inp);

}
    

#endif
