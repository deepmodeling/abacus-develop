#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_md()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    // 9. Molecular dynamics
    {
        Input_Item item("md_type");
        item.annotation = "choose ensemble";
        item.category = "Molecular dynamics";
        item.type = "String";
        item.description = R"(Control the algorithm to integrate the equation of motion for molecular dynamics (MD), see md.md in detail.

* fire: a MD-based relaxation algorithm, named fast inertial relaxation engine.
* nve: NVE ensemble with velocity Verlet algorithm.
* nvt: NVT ensemble, see md_thermostat in detail.
* npt: Nose-Hoover style NPT ensemble, see md_pmode in detail.
* langevin: NVT ensemble with Langevin thermostat, see md_damp in detail.
* msst: MSST method, see msst_direction, msst_vel, msst_qmass, msst_vis, msst_tscale in detail.)";
        item.default_value = "nvt";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.mdp.md_type);
        this->add_item(item);
    }
    {
        Input_Item item("md_nstep");
        item.annotation = "md steps";
        item.category = "Molecular dynamics";
        item.type = "Integer";
        item.description = "The total number of molecular dynamics steps.";
        item.default_value = "10";
        item.unit = "";
        item.availability = "";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.mdp.md_nstep == 0 && para.input.esolver_type != "tddft")
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
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The time step used in molecular dynamics calculations.";
        item.default_value = "1.0";
        item.unit = "fs";
        item.availability = "";
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.mdp.md_dt < 0) {
                ModuleBase::WARNING_QUIT("ReadInput", "time interval of MD calculation should be positive");
}
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.td_dt != -1.0)
            {
                GlobalV::ofs_running << "td_dt exist, set md_dt with td_dt" << std::endl;
                para.input.mdp.md_dt = para.input.td_dt * para.input.estep_per_md;
            }
        };
        read_sync_double(input.mdp.md_dt);
        this->add_item(item);
    }
    {
        Input_Item item("md_thermostat");
        item.annotation = "choose thermostat";
        item.category = "Molecular dynamics";
        item.type = "String";
        item.description = R"(Specify the temperature control method used in NVT ensemble.

* nhc: Nose-Hoover chain, see md_tfreq and md_tchain in detail.
* anderson: Anderson thermostat, see md_nraise in detail.
* berendsen: Berendsen thermostat, see md_nraise in detail.
* rescaling: velocity Rescaling method 1, see md_tolerance in detail.
* rescale_v: velocity Rescaling method 2, see md_nraise in detail.)";
        item.default_value = "nhc";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.mdp.md_thermostat);
        this->add_item(item);
    }
    {
        Input_Item item("md_tfirst");
        item.annotation = "temperature first";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = R"(The temperature used in molecular dynamics calculations.

If md_tfirst is unset or less than zero, init_vel is autoset to be true. If init_vel is true, the initial temperature will be determined by the velocities read from STRU. In this case, if velocities are unspecified in STRU, the initial temperature is set to zero.

If md_tfirst is set to a positive value and init_vel is true simultaneously, please make sure they are consistent, otherwise abacus will exit immediately.

Note that md_tlast is only used in NVT/NPT simulations. If md_tlast is unset or less than zero, md_tlast is set to md_tfirst. If md_tlast is set to be different from md_tfirst, ABACUS will automatically change the temperature from md_tfirst to md_tlast.)";
        item.default_value = "No default";
        item.unit = "K";
        item.availability = "";
        read_sync_double(input.mdp.md_tfirst);
        this->add_item(item);
    }
    {
        Input_Item item("md_tlast");
        item.annotation = "temperature last";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = R"(The temperature used in molecular dynamics calculations.

If md_tfirst is unset or less than zero, init_vel is autoset to be true. If init_vel is true, the initial temperature will be determined by the velocities read from STRU. In this case, if velocities are unspecified in STRU, the initial temperature is set to zero.

If md_tfirst is set to a positive value and init_vel is true simultaneously, please make sure they are consistent, otherwise abacus will exit immediately.

Note that md_tlast is only used in NVT/NPT simulations. If md_tlast is unset or less than zero, md_tlast is set to md_tfirst. If md_tlast is set to be different from md_tfirst, ABACUS will automatically change the temperature from md_tfirst to md_tlast.)";
        item.default_value = "No default";
        item.unit = "K";
        item.availability = "";
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
        Input_Item item("md_prec_level");
        item.annotation = "precision level for vc-md";
        item.category = "Molecular dynamics";
        item.type = "Integer";
        item.description = R"(Determine the precision level of variable-cell molecular dynamics calculations.
* 0: FFT grids do not change, only G vectors and K vectors are changed due to the change of lattice vector. This level is suitable for cases where the variation of the volume and shape is not large, and the efficiency is relatively higher.
* 2: FFT grids change per step. This level is suitable for cases where the variation of the volume and shape is large, such as the MSST method. However, accuracy comes at the cost of efficiency.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("md_restart");
        item.annotation = "whether restart";
        item.category = "Molecular dynamics";
        item.type = "Boolean";
        item.description = R"(Control whether to restart molecular dynamics calculations and time-dependent density functional theory calculations.
* True: ABACUS will read in {md_step}, then read in the corresponding STRU_MD_suffix/STRU/ automatically. For tddft, ABACUS will also read in WFC_NAO_K${kpoint} of the last step (You need to set out_wfc_lcao=1 and out_app_flag=0 to obtain this file).
* False: ABACUS will start molecular dynamics calculations normally from the first step.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.mdp.md_restart);
        this->add_item(item);
    }
    {
        Input_Item item("md_restartfreq");
        item.annotation = "The period to output MD restart information";
        item.category = "Molecular dynamics";
        item.type = "Integer";
        item.description = "The output frequency of OUT.{suffix}/STRIU/, which are used to restart molecular dynamics calculations, see md_restart in detail.";
        item.default_value = "5";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.mdp.md_restartfreq);
        this->add_item(item);
    }
    {
        Input_Item item("md_dumpfreq");
        item.annotation = "The period to dump MD information";
        item.category = "Molecular dynamics";
        item.type = "Integer";
        item.description = "The output frequency of OUT.${suffix}/MD_dump in molecular dynamics calculations, which including the information of lattices and atoms.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.mdp.md_dumpfreq);
        this->add_item(item);
    }
    {
        Input_Item item("dump_force");
        item.annotation = "output atomic forces into the file MD_dump or not";
        item.category = "Molecular dynamics";
        item.type = "Boolean";
        item.description = "Whether to output atomic forces into the file OUT.${suffix}/MD_dump.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.mdp.dump_force);
        this->add_item(item);
    }
    {
        Input_Item item("dump_vel");
        item.annotation = "output atomic velocities into the file MD_dump or not";
        item.category = "Molecular dynamics";
        item.type = "Boolean";
        item.description = "Whether to output atomic velocities into the file OUT.${suffix}/MD_dump.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.mdp.dump_vel);
        this->add_item(item);
    }
    {
        Input_Item item("dump_virial");
        item.annotation = "output lattice virial into the file MD_dump or not";
        item.category = "Molecular dynamics";
        item.type = "Boolean";
        item.description = "Whether to output lattice virials into the file OUT.${suffix}/MD_dump.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.mdp.dump_virial);
        this->add_item(item);
    }
    {
        Input_Item item("md_seed");
        item.annotation = "random seed for MD";
        item.category = "Molecular dynamics";
        item.type = "Integer";
        item.description = R"(The random seed to initialize random numbers used in molecular dynamics calculations.
* < 0: No srand() function is called.
* >= 0: The function srand(md_seed) is called.)";
        item.default_value = "-1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.mdp.md_seed);
        this->add_item(item);
    }
    {
        Input_Item item("md_tfreq");
        item.annotation = "oscillation frequency, used to determine qmass of NHC";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = R"(Control the frequency of temperature oscillations during the simulation. If it is too large, the temperature will fluctuate violently; if it is too small, the temperature will take a very long time to equilibrate with the atomic system.

Note: It is a system-dependent empirical parameter, ranging from 1/(40*md_dt) to 1/(100*md_dt). An improper choice might lead to the failure of jobs.)";
        item.default_value = "1/40/md_dt";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("md_tchain");
        item.annotation = "number of Nose-Hoover chains";
        item.category = "Molecular dynamics";
        item.type = "Integer";
        item.description = "Number of thermostats coupled with the particles in the NVT/NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.mdp.md_tchain);
        this->add_item(item);
    }
    {
        Input_Item item("md_pmode");
        item.annotation = "NPT ensemble mode: iso, aniso, tri";
        item.category = "Molecular dynamics";
        item.type = "String";
        item.description = R"(Determine the precision level of variable-cell molecular dynamics calculations.
* 0: FFT grids do not change, only G vectors and K vectors are changed due to the change of lattice vector. This level is suitable for cases where the variation of the volume and shape is not large, and the efficiency is relatively higher.
* 2: FFT grids change per step. This level is suitable for cases where the variation of the volume and shape is large, such as the MSST method. However, accuracy comes at the cost of efficiency.)";
        item.default_value = "iso";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.mdp.md_pmode);
        this->add_item(item);
    }
    {
        Input_Item item("ref_cell_factor");
        item.annotation = "construct a reference cell bigger than the initial cell";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "Construct a reference cell bigger than the initial cell. The reference cell has to be large enough so that the lattice vectors of the fluctuating cell do not exceed the reference lattice vectors during MD. Typically, 1.02 ~ 1.10 is sufficient. However, the cell fluctuations depend on the specific system and thermodynamic conditions. So users must test for a proper choice. This parameters should be used in conjunction with erf_ecut, erf_height, and erf_sigma.";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.ref_cell_factor);
        this->add_item(item);
    }
    {
        Input_Item item("md_pcouple");
        item.annotation = "whether couple different components: xyz, xy, yz, xz, none";
        item.category = "Molecular dynamics";
        item.type = "String";
        item.description = R"(The coupled lattice vectors will scale proportionally in NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.
* none: Three lattice vectors scale independently.
* xyz: Lattice vectors x, y, and z scale proportionally.
* xy: Lattice vectors x and y scale proportionally.
* xz: Lattice vectors x and z scale proportionally.
* yz: Lattice vectors y and z scale proportionally.)";
        item.default_value = "none";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("md_pfirst");
        item.annotation = "initial target pressure";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The target pressure used in NPT ensemble simulations, the default value of md_plast is md_pfirst. If md_plast is set to be different from md_pfirst, ABACUS will automatically change the target pressure from md_pfirst to md_plast.";
        item.default_value = "-1.0";
        item.unit = "kbar";
        item.availability = "";
        read_sync_double(input.mdp.md_pfirst);
        this->add_item(item);
    }
    {
        Input_Item item("md_plast");
        item.annotation = "final target pressure";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The target pressure used in NPT ensemble simulations, the default value of md_plast is md_pfirst. If md_plast is set to be different from md_pfirst, ABACUS will automatically change the target pressure from md_pfirst to md_plast.";
        item.default_value = "-1.0";
        item.unit = "kbar";
        item.availability = "";
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
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = R"(The frequency of pressure oscillations during the NPT ensemble simulation. If it is too large, the pressure will fluctuate violently; if it is too small, the pressure will take a very long time to equilibrate with the atomic system.

Note: It is a system-dependent empirical parameter. An improper choice might lead to the failure of jobs.)";
        item.default_value = "1/400/md_dt";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("md_pchain");
        item.annotation = "num of thermostats coupled with barostat";
        item.category = "Molecular dynamics";
        item.type = "Integer";
        item.description = "The number of thermostats coupled with the barostat in the NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.mdp.md_pchain);
        this->add_item(item);
    }
    {
        Input_Item item("lj_rule");
        item.annotation = "combination rules used to construct the parameter matrix for LJ potential";
        item.category = "Molecular dynamics";
        item.type = "Integer";
        item.description = "The Lennard-Jones potential between two atoms equals: $\\sigma_k\\sigma(i,j)$";
        item.default_value = "2";
        item.unit = "";
        item.availability = "";
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
        item.category = "Molecular dynamics";
        item.type = "Boolean";
        item.description = "It True, the LJ potential is shifted by a constant such that it is zero at the cut-off distance.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.mdp.lj_eshift);
        this->add_item(item);
    }
    {
        Input_Item item("lj_rcut");
        item.annotation = "cutoff radius of LJ potential";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "Cut-off radius for Leonard Jones potential, beyond which the interaction will be neglected. It can be a single value, which means that all pairs of atoms types share the same cut-off radius. Otherwise, it should be a multiple-component vector, containing values, see details in lj_rule.";
        item.default_value = "No default";
        item.unit = "Angstrom";
        item.availability = "";
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
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The vector representing the matrix for Leonard Jones potential. See details in lj_rule.";
        item.default_value = "No default";
        item.unit = "eV";
        item.availability = "";
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
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The vector representing the matrix for Leonard Jones potential. See details in lj_rule.";
        item.default_value = "No default";
        item.unit = "Angstrom";
        item.availability = "";
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
        item.category = "Molecular dynamics";
        item.type = "String";
        item.description = "The filename of DP/NEP potential files, see md.md in detail.";
        item.default_value = "graph.pb";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.mdp.pot_file);
        this->add_item(item);
    }
    {
        Input_Item item("dp_rescaling");
        item.annotation = "rescaling factor for dp potential";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "Rescaling factor to use a temperature-dependent DP. Energy, stress and force calculated by DP will be multiplied by this factor.";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "esolver_type = dp.";
        read_sync_double(input.mdp.dp_rescaling);
        this->add_item(item);
    }
    {
        Input_Item item("dp_fparam");
        item.annotation = "the frame parameter for dp potential";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The frame parameter for dp potential. The array size is dim_fparam, then all frames are assumed to be provided with the same fparam.";
        item.default_value = "{}";
        item.unit = "";
        item.availability = "esolver_type = dp.";
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
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The atomic parameter for dp potential. The array size can be (1) natoms x dim_aparam, then all frames are assumed to be provided with the same aparam; (2) dim_aparam, then all frames and atoms are assumed to be provided with the same aparam.";
        item.default_value = "{}";
        item.unit = "";
        item.availability = "esolver_type = dp.";
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
        item.category = "Molecular dynamics";
        item.type = "Integer";
        item.description = R"(The direction of the shock wave in the MSST method.
* 0: x direction
* 1: y direction
* 2: z direction)";
        item.default_value = "2";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.mdp.msst_direction);
        this->add_item(item);
    }
    {
        Input_Item item("msst_vel");
        item.annotation = "the velocity of shock wave";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The velocity of the shock wave in the MSST method.";
        item.default_value = "0.0";
        item.unit = "Angstrom/fs";
        item.availability = "";
        read_sync_double(input.mdp.msst_vel);
        this->add_item(item);
    }
    {
        Input_Item item("msst_vis");
        item.annotation = "artificial viscosity";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "Artificial viscosity in the MSST method.";
        item.default_value = "0.0";
        item.unit = "g/(mol*Angstrom*fs)";
        item.availability = "";
        read_sync_double(input.mdp.msst_vis);
        this->add_item(item);
    }
    {
        Input_Item item("msst_tscale");
        item.annotation = "reduction in initial temperature";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The reduction percentage of the initial temperature used to compress volume in the MSST method.";
        item.default_value = "0.01";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.mdp.msst_tscale);
        this->add_item(item);
    }
    {
        Input_Item item("msst_qmass");
        item.annotation = "mass of thermostat";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "Inertia of the extended system variable. You should set a number larger than 0.";
        item.default_value = "No default";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("md_damp");
        item.annotation = "damping parameter (time units) used to add force in "
                          "Langevin method";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The damping parameter used to add fictitious force in the Langevin method.";
        item.default_value = "1.0";
        item.unit = "fs";
        item.availability = "";
        read_sync_double(input.mdp.md_damp);
        this->add_item(item);
    }
    {
        Input_Item item("md_tolerance");
        item.annotation = "tolerance for velocity rescaling (K)";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The temperature tolerance for velocity rescaling. Velocities are rescaled if the current and target temperature differ more than md_tolerance.";
        item.default_value = "100.0";
        item.unit = "K";
        item.availability = "";
        read_sync_double(input.mdp.md_tolerance);
        this->add_item(item);
    }
    {
        Input_Item item("md_nraise");
        item.annotation = "parameters used when md_type=nvt";
        item.category = "Molecular dynamics";
        item.type = "Integer";
        item.description = R"(* Anderson: The "collision frequency" parameter is given as 1/md_nraise.
* Berendsen: The "rise time" parameter is given in units of the time step: tau = md_nraise*md_dt, so md_dt/tau = 1/md_nraise.
* Rescale_v: Every md_nraise steps the current temperature is rescaled to the target temperature.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.mdp.md_nraise);
        this->add_item(item);
    }
    {
        Input_Item item("cal_syns");
        item.annotation = "calculate asynchronous overlap matrix to output for Hefei-NAMD";
        item.category = "Molecular dynamics";
        item.type = R"(Boolean [Integer](optional))";
        item.description = R"(Whether to calculate and output asynchronous overlap matrix for Hefei-NAMD interface. When enabled, calculates <phi(t-1)|phi(t)> by computing overlap between basis functions at atomic positions from previous time step and current time step. The overlap is calculated by shifting atom positions backward by velocity x md_dt. Output file: OUT.*/syns_nao.csr in CSR format.

* 0 or false: disable
* 1 or true: enable with default precision (8 digits)
* 1 5: enable with custom precision (5 digits)

[NOTE] Only works with LCAO basis and molecular dynamics calculations. Requires atomic velocities. Output starts from the second MD step (istep > 0).)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count < 1) ModuleBase::WARNING_QUIT("ReadInput", "cal_syns needs at least 1 value");
            para.input.cal_syns[0] = assume_as_boolean(item.str_values[0]);
            para.input.cal_syns[1] = 8;
            if (count >= 2) try { para.input.cal_syns[1] = std::stoi(item.str_values[1]); }
            catch (const std::invalid_argument&) { /* do nothing */ }
            catch (const std::out_of_range&) {/* do nothing */}
        };

        sync_intvec(input.cal_syns, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("dmax");
        item.annotation = "maximum displacement of all atoms in one step (bohr)";
        item.category = "Molecular dynamics";
        item.type = "Real";
        item.description = "The maximum displacement of all atoms in one step. This parameter is useful when cal_syns = True.";
        item.default_value = "0.01";
        item.unit = "bohr";
        item.availability = "";
        read_sync_double(input.dmax);
        this->add_item(item);
    }
}
} // namespace ModuleIO
