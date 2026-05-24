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

#include "source_io/module_chgpot/write_init.h"
#include "source_io/module_output/cube_io.h"
#include "source_base/tool_quit.h"

#include <sstream>
#include <cassert>

void ModuleIO::write_chg_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    const Charge &chr,
    const elecstate::Efermi &efermi,
    const int istep,
    const std::string& out_dir,
    const Input_para& inp)
{
    const int nspin = inp.nspin;
    assert(nspin == 1 || nspin == 2 || nspin == 4);

    if (istep < 0)
    {
        ModuleBase::WARNING_QUIT("write_chg_init", "istep must be >= 0");
    }

    if (inp.out_chg[0] == 2)
    {
        bool should_output = (istep == 0) ||
                             (inp.out_freq_ion > 0 && istep % inp.out_freq_ion == 0);

        if (should_output)
        {
            int geom_step = istep + 1;

            for (int is = 0; is < nspin; is++)
            {
                std::stringstream ss;
                ss << out_dir << "chg";

                if (nspin == 1)
                {
                    ss << "g" << geom_step << "_ini.cube";
                }
                else if (nspin == 2 || nspin == 4)
                {
                    ss << "s" << is + 1 << "g" << geom_step << "_ini.cube";
                }

                double fermi_energy = 0.0;
                if (nspin == 1 || nspin == 4)
                {
                    fermi_energy = efermi.ef;
                }
                else if (nspin == 2)
                {
                    if (is == 0)
                    {
                        fermi_energy = efermi.ef_up;
                    }
                    else if (is == 1)
                    {
                        fermi_energy = efermi.ef_dw;
                    }
                }

                ModuleIO::write_vdata_palgrid(para_grid,
                                              chr.rho[is],
                                              is,
                                              nspin,
                                              istep,
                                              ss.str(),
                                              fermi_energy,
                                              &(ucell),
                                              inp.out_chg[1]);
            }
        }
    }
    return;
}


void ModuleIO::write_pot_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    elecstate::ElecState *pelec,
    const int istep,
    const std::string& out_dir,
    const Input_para& inp)
{
    const int nspin = inp.nspin;
    assert(nspin == 1 || nspin == 2 || nspin == 4);

    if (istep < 0)
    {
        ModuleBase::WARNING_QUIT("write_pot_init", "istep must be >= 0");
    }

    if (inp.out_pot[0] == 3)
    {
        bool should_output = (istep == 0) ||
                             (inp.out_freq_ion > 0 && istep % inp.out_freq_ion == 0);

        if (should_output)
        {
            int geom_step = istep + 1;

            for (int is = 0; is < nspin; is++)
            {
                std::stringstream ss;
                ss << out_dir << "pot";

                if (nspin == 1)
                {
                    ss << "g" << geom_step << "_ini.cube";
                }
                else if (nspin == 2 || nspin == 4)
                {
                    ss << "s" << is + 1 << "g" << geom_step << "_ini.cube";
                }

                ModuleIO::write_vdata_palgrid(para_grid,
                                              pelec->pot->get_eff_v(is),
                                              is,
                                              nspin,
                                              istep,
                                              ss.str(),
                                              0.0,
                                              &(ucell),
                                              11,
                                              0);
            }
        }
    }
    return;
}
