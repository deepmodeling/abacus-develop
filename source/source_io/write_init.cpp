#include "source_io/write_init.h"
#include "source_io/cube_io.h"

void ModuleIO::write_chg_init(
    const UnitCell& ucell,
    const Parallel_Grid &para_grid,
    const Charge &chr,
    const elecstate::Efermi &efermi,
    const int istep,
    const Input_para& inp)
{
    const int nspin = inp.nspin;
    assert(nspin == 1 || nspin ==2 || nspin == 4);

    if (inp.out_chg[0] == 2)
    {
        for (int is = 0; is < nspin; is++)
        {
            std::stringstream ss;
            ss << PARAM.globalv.global_out_dir << "chg";

            if(nspin==1)
            {
                ss << "ini.cube";
            }
            else if(nspin==2 || nspin==4)
            {
                ss << "s" << is + 1 << "ini.cube";
            }

            // mohan add 2025-10-18
            double fermi_energy = 0.0;
            if(nspin == 1 || nspin ==4)
			{
				fermi_energy = efermi.ef;
			}
			else if(nspin == 2)
			{
				if(is==0) 
				{
					fermi_energy = efermi.ef_up;
				}
				else if(is==1)
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
                                          &(ucell));
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
    //! output total local potential of the initial charge density
    const int nspin = inp.nspin;
    assert(nspin == 1 || nspin ==2 || nspin == 4);

    if (inp.out_pot == 3)
    {
        for (int is = 0; is < nspin; is++)
        {
            std::stringstream ss;
            ss << PARAM.globalv.global_out_dir << "pot";

            if(nspin==1)
            {
                ss << "ini.cube";
            }
            else if(nspin==2 || nspin==4)
            {
                ss << "s" << is + 1 << "ini.cube";
            }

            ModuleIO::write_vdata_palgrid(para_grid,
                                          pelec->pot->get_effective_v(is),
                                          is,
                                          nspin,
                                          istep,
                                          ss.str(),
                                          0.0, // efermi
                                          &(ucell),
                                          11, // precsion
                                          0); // out_fermi
        }
    }

}
