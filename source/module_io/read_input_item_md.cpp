#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_md()
{
    // 9. Molecular dynamics
    {
        Input_Item item("md_type");
        item.annotation = "choose ensemble";
        read_sync_string(mdp.md_type);
        this->add_item(item);
    }
    {
        Input_Item item("md_thermostat");
        item.annotation = "choose thermostat";
        read_sync_string(mdp.md_thermostat);
        this->add_item(item);
    }
    {
        Input_Item item("md_nstep");
        item.annotation = "md steps";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.mdp.md_nstep == 0)
            {
                GlobalV::ofs_running << "md_nstep should be set. Autoset md_nstep to 50!" << std::endl;
                para.input.mdp.md_nstep = 50;
            }
        };
        read_sync_int(mdp.md_nstep);
        this->add_item(item);
    }
    {
        Input_Item item("md_dt");
        item.annotation = "time step";
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.mdp.md_dt < 0)
                ModuleBase::WARNING_QUIT("ReadInput", "time interval of MD calculation should be positive");
        };
        read_sync_double(mdp.md_dt);
        this->add_item(item);
    }
    {
        Input_Item item("md_tchain");
        item.annotation = "number of Nose-Hoover chains";
        read_sync_int(mdp.md_tchain);
        this->add_item(item);
    }
    {
        Input_Item item("md_tfirst");
        item.annotation = "temperature first";
        read_sync_double(mdp.md_tfirst);
        this->add_item(item);
    }
    {
        Input_Item item("md_tlast");
        item.annotation = "temperature last";
        read_sync_double(mdp.md_tlast);
        this->add_item(item);
    }
    {
        Input_Item item("md_dumpfreq");
        item.annotation = "The period to dump MD information";
        read_sync_int(mdp.md_dumpfreq);
        this->add_item(item);
    }
    {
        Input_Item item("md_restartfreq");
        item.annotation = "The period to output MD restart information";
        read_sync_int(mdp.md_restartfreq);
        this->add_item(item);
    }
    {
        Input_Item item("md_seed");
        item.annotation = "random seed for MD";
        read_sync_int(mdp.md_seed);
        this->add_item(item);
    }
    {
        Input_Item item("md_prec_level");
        item.annotation = "precision level for vc-md";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation != "md")
            {
                para.input.mdp.md_prec_level = 0;
            }
            // md_prec_level only used in vc-md  liuyu 2023-03-27
            else if (para.input.mdp.md_type != "msst" && para.input.mdp.md_type != "npt")
            {
                para.input.mdp.md_prec_level = 0;
            }
        };
        read_sync_int(mdp.md_prec_level);
        this->add_item(item);
    }
    {
        Input_Item item("ref_cell_factor");
        item.annotation = "construct a reference cell bigger than the initial cell";
        read_sync_double(ref_cell_factor);
        this->add_item(item);
    }
    {
        Input_Item item("md_restart");
        item.annotation = "whether restart";
        read_sync_bool(mdp.md_restart);
        this->add_item(item);
    }
    {
        Input_Item item("lj_rcut");
        item.annotation = "cutoff radius of LJ potential";
        read_sync_double(mdp.lj_rcut);
        this->add_item(item);
    }
    {
        Input_Item item("lj_epsilon");
        item.annotation = "the value of epsilon for LJ potential";
        read_sync_double(mdp.lj_epsilon);
        this->add_item(item);
    }
    {
        Input_Item item("lj_sigma");
        item.annotation = "the value of sigma for LJ potential";
        read_sync_double(mdp.lj_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("pot_file");
        item.annotation = "the filename of potential files for CMD such as DP";
        read_sync_string(mdp.pot_file);
        this->add_item(item);
    }
    {
        Input_Item item("msst_direction");
        item.annotation = "the direction of shock wave";
        read_sync_int(mdp.msst_direction);
        this->add_item(item);
    }
    {
        Input_Item item("msst_vel");
        item.annotation = "the velocity of shock wave";
        read_sync_double(mdp.msst_vel);
        this->add_item(item);
    }
    {
        Input_Item item("msst_vis");
        item.annotation = "artificial viscosity";
        read_sync_double(mdp.msst_vis);
        this->add_item(item);
    }
    {
        Input_Item item("msst_tscale");
        item.annotation = "reduction in initial temperature";
        read_sync_double(mdp.msst_tscale);
        this->add_item(item);
    }
    {
        Input_Item item("msst_qmass");
        item.annotation = "mass of thermostat";
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.mdp.msst_qmass <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "msst_qmass must be greater than 0!");
            }
        };
        read_sync_double(mdp.msst_qmass);
        this->add_item(item);
    }
    {
        Input_Item item("md_tfreq");
        item.annotation = "oscillation frequency, used to determine qmass of NHC";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.mdp.md_tfreq == 0 && para.input.calculation == "md")
            {
                para.input.mdp.md_tfreq = 1.0 / 40 / para.input.mdp.md_dt;
            }
        };
        read_sync_double(mdp.md_tfreq);
        this->add_item(item);
    }
    {
        Input_Item item("md_damp");
        item.annotation = "damping parameter (time units) used to add force in "
                          "Langevin method";
        read_sync_double(mdp.md_damp);
        this->add_item(item);
    }
    {
        Input_Item item("md_nraise");
        item.annotation = "parameters used when md_type=nvt";
        read_sync_int(mdp.md_nraise);
        this->add_item(item);
    }
    {
        Input_Item item("cal_syns");
        item.annotation = "calculate asynchronous overlap matrix to output for Hefei-NAMD";
        read_sync_bool(cal_syns);
        this->add_item(item);
    }
    {
        Input_Item item("dmax");
        item.annotation = "maximum displacement of all atoms in one step (bohr)";
        read_sync_double(dmax);
        this->add_item(item);
    }
    {
        Input_Item item("md_tolerance");
        item.annotation = "tolerance for velocity rescaling (K)";
        read_sync_double(mdp.md_tolerance);
        this->add_item(item);
    }
    {
        Input_Item item("md_pmode");
        item.annotation = "NPT ensemble mode: iso, aniso, tri";
        read_sync_string(mdp.md_pmode);
        this->add_item(item);
    }
    {
        Input_Item item("md_pcouple");
        item.annotation = "whether couple different components: xyz, xy, yz, xz, none";
        read_sync_string(mdp.md_pcouple);
        this->add_item(item);
    }
    {
        Input_Item item("md_pchain");
        item.annotation = "num of thermostats coupled with barostat";
        read_sync_int(mdp.md_pchain);
        this->add_item(item);
    }
    {
        Input_Item item("md_pfirst");
        item.annotation = "initial target pressure";
        read_sync_double(mdp.md_pfirst);
        this->add_item(item);
    }
    {
        Input_Item item("md_plast");
        item.annotation = "final target pressure";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (!item.is_read()) // no md_plast in INPUT
                para.input.mdp.md_plast = para.input.mdp.md_pfirst;
        };
        read_sync_double(mdp.md_plast);
        this->add_item(item);
    }
    {
        Input_Item item("md_pfreq");
        item.annotation = "oscillation frequency, used to determine qmass of "
                          "thermostats coupled with barostat";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.mdp.md_pfreq == 0 && para.input.calculation == "md")
            {
                para.input.mdp.md_pfreq = 1.0 / 400 / para.input.mdp.md_dt;
            }
        };
        read_sync_double(mdp.md_pfreq);
        this->add_item(item);
    }
    {
        Input_Item item("dump_force");
        item.annotation = "output atomic forces into the file MD_dump or not";
        read_sync_bool(mdp.dump_force);
        this->add_item(item);
    }
    {
        Input_Item item("dump_vel");
        item.annotation = "output atomic velocities into the file MD_dump or not";
        read_sync_bool(mdp.dump_vel);
        this->add_item(item);
    }
    {
        Input_Item item("dump_virial");
        item.annotation = "output lattice virial into the file MD_dump or not";
        read_sync_bool(mdp.dump_virial);
        this->add_item(item);
    }
}
} // namespace ModuleIO