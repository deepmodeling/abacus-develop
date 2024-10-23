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
        read_sync_string(input.mdp.md_type);
        this->add_item(item);
    }
    {
        Input_Item item("md_thermostat");
        item.annotation = "choose thermostat";
        read_sync_string(input.mdp.md_thermostat);
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
        read_sync_int(input.mdp.md_nstep);
        this->add_item(item);
    }
    {
        Input_Item item("md_dt");
        item.annotation = "time step";
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.mdp.md_dt < 0) {
                ModuleBase::WARNING_QUIT("ReadInput", "time interval of MD calculation should be positive");
}
        };
        read_sync_double(input.mdp.md_dt);
        this->add_item(item);
    }
    {
        Input_Item item("md_tchain");
        item.annotation = "number of Nose-Hoover chains";
        read_sync_int(input.mdp.md_tchain);
        this->add_item(item);
    }
    {
        Input_Item item("md_tfirst");
        item.annotation = "temperature first";
        read_sync_double(input.mdp.md_tfirst);
        this->add_item(item);
    }
    {
        Input_Item item("md_tlast");
        item.annotation = "temperature last";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.mdp.md_tlast < 0)
            {
                para.input.mdp.md_tlast = para.mdp.md_tfirst;
            }
        };
        read_sync_double(input.mdp.md_tlast);
        this->add_item(item);
    }
    {
        Input_Item item("md_dumpfreq");
        item.annotation = "The period to dump MD information";
        read_sync_int(input.mdp.md_dumpfreq);
        this->add_item(item);
    }
    {
        Input_Item item("md_restartfreq");
        item.annotation = "The period to output MD restart information";
        read_sync_int(input.mdp.md_restartfreq);
        this->add_item(item);
    }
    {
        Input_Item item("md_seed");
        item.annotation = "random seed for MD";
        read_sync_int(input.mdp.md_seed);
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
        read_sync_int(input.mdp.md_prec_level);
        this->add_item(item);
    }
    {
        Input_Item item("ref_cell_factor");
        item.annotation = "construct a reference cell bigger than the initial cell";
        read_sync_double(input.ref_cell_factor);
        this->add_item(item);
    }
    {
        Input_Item item("md_restart");
        item.annotation = "whether restart";
        read_sync_bool(input.mdp.md_restart);
        this->add_item(item);
    }
    {
        Input_Item item("lj_rule");
        item.annotation = "combination rules used to construct the parameter matrix for LJ potential";
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.esolver_type == "lj" && para.input.mdp.lj_rule != 1 && para.input.mdp.lj_rule != 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "lj_rule must be 1 or 2");
            }
        };
        read_sync_int(input.mdp.lj_rule);
        this->add_item(item);
    }
    {
        Input_Item item("lj_eshift");
        item.annotation = "whether to use energy shift for LJ potential";
        read_sync_bool(input.mdp.lj_eshift);
        this->add_item(item);
    }
    {
        Input_Item item("lj_rcut");
        item.annotation = "cutoff radius of LJ potential";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            para.input.mdp.lj_rcut.resize(count);
            std::transform(begin(item.str_values),
                           end(item.str_values),
                           begin(para.input.mdp.lj_rcut),
                           [](std::string str) { return std::stod(str); });
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (!item.is_read()) {
                return;
}
            size_t n_ljrcut = para.input.mdp.lj_rcut.size();
            if (n_ljrcut != 1 && n_ljrcut != para.input.ntype * (para.input.ntype + 1) / 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", " the number of lj_rcut should be 1 or ntype(ntype+1)/2 ");
            }
            for (auto rcut: para.input.mdp.lj_rcut)
            {
                if (rcut <= 0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "lj_rcut must > 0");
                }
            }
        };
        sync_doublevec(input.mdp.lj_rcut, para.input.mdp.lj_rcut.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("lj_epsilon");
        item.annotation = "the value of epsilon for LJ potential";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            para.input.mdp.lj_epsilon.resize(count);
            std::transform(begin(item.str_values),
                           end(item.str_values),
                           begin(para.input.mdp.lj_epsilon),
                           [](std::string str) { return std::stod(str); });
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (!item.is_read()) {
                return;
}
            size_t n_ljepsilon = para.input.mdp.lj_epsilon.size();
            if (n_ljepsilon != para.input.ntype && n_ljepsilon != para.input.ntype * (para.input.ntype + 1) / 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", " the number of lj_epsilon should be ntype or ntype(ntype+1)/2 ");
            }
        };
        sync_doublevec(input.mdp.lj_epsilon, para.input.mdp.lj_epsilon.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("lj_sigma");
        item.annotation = "the value of sigma for LJ potential";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            para.input.mdp.lj_sigma.resize(count);
            std::transform(begin(item.str_values),
                           end(item.str_values),
                           begin(para.input.mdp.lj_sigma),
                           [](std::string str) { return std::stod(str); });
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (!item.is_read()) {
                return;
}
            size_t n_ljsigma = para.input.mdp.lj_sigma.size();
            if (n_ljsigma != para.input.ntype && n_ljsigma != para.input.ntype * (para.input.ntype + 1) / 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", " the number of lj_sigma should be ntype or ntype(ntype+1)/2 ");
            }
        };
        sync_doublevec(input.mdp.lj_sigma, para.input.mdp.lj_sigma.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("pot_file");
        item.annotation = "the filename of potential files for CMD such as DP";
        read_sync_string(input.mdp.pot_file);
        this->add_item(item);
    }
    {
        Input_Item item("dp_rescaling");
        item.annotation = "rescaling factor for dp potential";
        read_sync_double(input.mdp.dp_rescaling);
        this->add_item(item);
    }
    {
        Input_Item item("dp_fparam");
        item.annotation = "the frame parameter for dp potential";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            para.input.mdp.dp_fparam.resize(count);
            std::transform(begin(item.str_values),
                           end(item.str_values),
                           begin(para.input.mdp.dp_fparam),
                           [](std::string str) { return std::stod(str); });
        };
        sync_doublevec(input.mdp.dp_fparam, para.input.mdp.dp_fparam.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("dp_aparam");
        item.annotation = "the atomic parameter for dp potential";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            para.input.mdp.dp_aparam.resize(count);
            std::transform(begin(item.str_values),
                           end(item.str_values),
                           begin(para.input.mdp.dp_aparam),
                           [](std::string str) { return std::stod(str); });
        };
        sync_doublevec(input.mdp.dp_aparam, para.input.mdp.dp_aparam.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("msst_direction");
        item.annotation = "the direction of shock wave";
        read_sync_int(input.mdp.msst_direction);
        this->add_item(item);
    }
    {
        Input_Item item("msst_vel");
        item.annotation = "the velocity of shock wave";
        read_sync_double(input.mdp.msst_vel);
        this->add_item(item);
    }
    {
        Input_Item item("msst_vis");
        item.annotation = "artificial viscosity";
        read_sync_double(input.mdp.msst_vis);
        this->add_item(item);
    }
    {
        Input_Item item("msst_tscale");
        item.annotation = "reduction in initial temperature";
        read_sync_double(input.mdp.msst_tscale);
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
        read_sync_double(input.mdp.msst_qmass);
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
        read_sync_double(input.mdp.md_tfreq);
        this->add_item(item);
    }
    {
        Input_Item item("md_damp");
        item.annotation = "damping parameter (time units) used to add force in "
                          "Langevin method";
        read_sync_double(input.mdp.md_damp);
        this->add_item(item);
    }
    {
        Input_Item item("md_nraise");
        item.annotation = "parameters used when md_type=nvt";
        read_sync_int(input.mdp.md_nraise);
        this->add_item(item);
    }
    {
        Input_Item item("cal_syns");
        item.annotation = "calculate asynchronous overlap matrix to output for Hefei-NAMD";
        read_sync_bool(input.cal_syns);
        this->add_item(item);
    }
    {
        Input_Item item("dmax");
        item.annotation = "maximum displacement of all atoms in one step (bohr)";
        read_sync_double(input.dmax);
        this->add_item(item);
    }
    {
        Input_Item item("md_tolerance");
        item.annotation = "tolerance for velocity rescaling (K)";
        read_sync_double(input.mdp.md_tolerance);
        this->add_item(item);
    }
    {
        Input_Item item("md_pmode");
        item.annotation = "NPT ensemble mode: iso, aniso, tri";
        read_sync_string(input.mdp.md_pmode);
        this->add_item(item);
    }
    {
        Input_Item item("md_pcouple");
        item.annotation = "whether couple different components: xyz, xy, yz, xz, none";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.mdp.md_pmode == "iso")
            {
                para.input.mdp.md_pcouple = "xyz";
            }
        };
        read_sync_string(input.mdp.md_pcouple);
        this->add_item(item);
    }
    {
        Input_Item item("md_pchain");
        item.annotation = "num of thermostats coupled with barostat";
        read_sync_int(input.mdp.md_pchain);
        this->add_item(item);
    }
    {
        Input_Item item("md_pfirst");
        item.annotation = "initial target pressure";
        read_sync_double(input.mdp.md_pfirst);
        this->add_item(item);
    }
    {
        Input_Item item("md_plast");
        item.annotation = "final target pressure";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (!item.is_read()) { // no md_plast in INPUT
                para.input.mdp.md_plast = para.input.mdp.md_pfirst;
}
        };
        read_sync_double(input.mdp.md_plast);
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
        read_sync_double(input.mdp.md_pfreq);
        this->add_item(item);
    }
    {
        Input_Item item("dump_force");
        item.annotation = "output atomic forces into the file MD_dump or not";
        read_sync_bool(input.mdp.dump_force);
        this->add_item(item);
    }
    {
        Input_Item item("dump_vel");
        item.annotation = "output atomic velocities into the file MD_dump or not";
        read_sync_bool(input.mdp.dump_vel);
        this->add_item(item);
    }
    {
        Input_Item item("dump_virial");
        item.annotation = "output lattice virial into the file MD_dump or not";
        read_sync_bool(input.mdp.dump_virial);
        this->add_item(item);
    }
}
} // namespace ModuleIO