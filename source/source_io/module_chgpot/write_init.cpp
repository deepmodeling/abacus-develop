#include "source_io/module_chgpot/write_init.h"
#include "source_io/module_output/cube_io.h"

#include <sstream>
#include <cassert>

void ModuleIO::write_chg_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    const Charge &chr,
    const elecstate::Efermi &efermi,
    const int istep,
    const Input_para& inp)
{
    const int nspin = inp.nspin;
    assert(nspin == 1 || nspin == 2 || nspin == 4);

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
                ss << PARAM.globalv.global_out_dir << "chg";

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
    const Input_para& inp)
{
    const int nspin = inp.nspin;
    assert(nspin == 1 || nspin == 2 || nspin == 4);

    if (inp.out_pot[0] == 3)
    {
        for (int is = 0; is < nspin; is++)
        {
            std::stringstream ss;
            ss << PARAM.globalv.global_out_dir << "pot";

            if (nspin == 1)
            {
                ss << "ini.cube";
            }
            else if (nspin == 2 || nspin == 4)
            {
                ss << "s" << is + 1 << "ini.cube";
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
