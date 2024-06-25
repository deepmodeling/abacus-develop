#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

#include <fstream>
#include <unistd.h>

namespace ModuleIO
{
// There are some examples:
// Generallly:
// {
//      Input_Item item("suffix");
//      item.annotation = "the name of main output directory";
//      read_sync_string(suffix);
//      this->add_item(item);
// }
//
// Specially:
// {
//      Input_Item item("kspacing");
//      item.annotation = "unit in 1/bohr, should be > 0, default is 0 which means read KPT file";
//
//      item.readvalue = [](const Input_Item& item, Parameter& para) {
//          para.input.kspacing[0] = std::stod(item.str_values[0]);
//          para.input.kspacing[1] = std::stod(item.str_values[1]);
//          para.input.kspacing[2] = std::stod(item.str_values[2]);
//      };
//
//      item.resetvalue = [](const Input_Item& item, Parameter& para) {
//          if(para.input.kspacing[0] <= 0) para.input.kspacing[0] = 1;
//      };
//
//      item.checkvalue = [](const Input_Item& item, const Parameter& para) {assert(para.input.kspacing[0]>0);};
//
//      item.getfinalvalue = [](Input_Item& item, const Parameter& para) {
//          item.final_value << para.input.kspacing[0] << " " << para.input.kspacing[1] << " " <<
//          para.input.kspacing[2];
//      };
//
//      add_doublevec_bcast(&Parameter::PARAMETER, N);
//      this->add_item(item);
//  }
void ReadInput::item_general()
{
    {
        Input_Item item("suffix");
        item.annotation = "the name of main output directory";
        read_sync_string(suffix);
        this->add_item(item);
    }
    {
        Input_Item item("latname");
        item.annotation = "the name of lattice name";
        read_sync_string(latname);
        this->add_item(item);
    }
    {
        Input_Item item("stru_file");
        item.annotation = "the filename of file containing atom positions";
        read_sync_string(stru_file);
        this->add_item(item);
    }
    {
        Input_Item item("kpoint_file");
        item.annotation = "the name of file containing k points";
        read_sync_string(kpoint_file);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_dir");
        item.annotation = "the directory containing pseudo files";
        read_sync_string(pseudo_dir);
        this->add_item(item);
    }
    {
        Input_Item item("orbital_dir");
        item.annotation = "the directory containing orbital files";
        read_sync_string(orbital_dir);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_rcut");
        item.annotation = "default #exchange correlation functional";
        read_sync_double(pseudo_rcut);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_mesh");
        item.annotation = "0: use our own mesh to do radial renormalization; 1: use mesh as in QE";
        read_sync_bool(pseudo_mesh);
        this->add_item(item);
    }
    {
        Input_Item item("lmaxmax");
        item.annotation = "maximum of l channels used";
        read_sync_int(lmaxmax);
        this->add_item(item);
    }
    {
        Input_Item item("dft_functional");
        item.annotation = "exchange correlation functional";
        read_sync_string(dft_functional);
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
#ifndef USE_PAW
            if (para.input.use_paw)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "to use PAW, compile with USE_PAW");
                if (para.input.basis_type != "pw")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "PAW is for pw basis only");
                }
            }
#endif
            if (para.input.dft_functional == "default" && para.input.use_paw)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dft_functional must be set when use_paw is true");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("xc_temperature");
        item.annotation = "temperature for finite temperature functionals";
        read_sync_double(xc_temperature);
        this->add_item(item);
    }
    {
        Input_Item item("calculation");
        item.annotation = "test; scf; relax; nscf; get_wf; get_pchg";
        read_sync_string(calculation);
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            std::string& calculation = para.input.calculation;
            if (calculation != "md")
            {
                para.input.mdp.md_prec_level = 0;
            }
            if (calculation == "scf")
            {
                if (para.input.mem_saver == 1)
                {
                    para.input.mem_saver = 0;
                    ModuleBase::GlobalFunc::AUTO_SET("mem_saver", "0");
                }
                para.input.relax_nmax = 1;
            }
            else if (calculation == "relax") // pengfei 2014-10-13
            {
                if (para.input.mem_saver == 1)
                {
                    para.input.mem_saver = 0;
                    ModuleBase::GlobalFunc::AUTO_SET("mem_saver", "0");
                }
                para.input.cal_force = 1;
                if (!para.input.relax_nmax)
                    para.input.relax_nmax = 50;
            }
            else if (calculation == "nscf" || calculation == "get_S")
            {
                GlobalV::CALCULATION = "nscf";
                para.input.relax_nmax = 1;
                para.input.out_stru = 0;

                if (para.input.basis_type == "pw" && calculation == "get_S")
                {
                    if (para.input.pw_diag_thr > 1.0e-3)
                    {
                        para.input.pw_diag_thr = 1.0e-5;
                    }
                }
                if (para.input.cal_force) // mohan add 2010-09-07
                {
                    para.input.cal_force = false;
                    ModuleBase::GlobalFunc::AUTO_SET("cal_force", "false");
                }
                if (para.input.init_chg != "file")
                {
                    para.input.init_chg = "file";
                    ModuleBase::GlobalFunc::AUTO_SET("init_chg", para.input.init_chg);
                }
            }
            else if (calculation == "get_pchg")
            {
                GlobalV::CALCULATION = "get_pchg";
                para.input.relax_nmax = 1;
                para.input.out_stru = 0;
                para.input.out_dos = 0;
                para.input.out_band[0] = 0;
                para.input.out_proj_band = 0;
                para.input.cal_force = 0;
                para.input.init_wfc = "file";
                para.input.init_chg = "atomic";   // useless,
                para.input.chg_extrap = "atomic"; // xiaohui modify 2015-02-01
                para.input.out_chg = 1;           // this leads to the calculation of state charge.
                para.input.out_dm = 0;
                para.input.out_dm1 = 0;
                para.input.out_pot = 0;
            }
            else if (calculation == "get_wf")
            {
                GlobalV::CALCULATION = "get_wf"; // mohan fix 2011-11-04
                para.input.relax_nmax = 1;
                para.input.out_stru = 0;
                para.input.out_dos = 0;
                para.input.out_band[0] = 0;
                para.input.out_proj_band = 0;
                para.input.cal_force = 0;
                para.input.init_wfc = "file";
                para.input.init_chg = "atomic";
                para.input.chg_extrap = "atomic"; // xiaohui modify 2015-02-01
                para.input.out_chg = 1;
                para.input.out_dm = 0;
                para.input.out_dm1 = 0;
                para.input.out_pot = 0;
            }
            else if (calculation == "md") // mohan add 2011-11-04
            {
                GlobalV::CALCULATION = "md";
                para.input.symmetry = "0";
                para.input.cal_force = 1;
                if (para.input.mdp.md_nstep == 0)
                {
                    GlobalV::ofs_running << "md_nstep should be set. Autoset md_nstep to 50!" << std::endl;
                    para.input.mdp.md_nstep = 50;
                }
                if (!para.input.sup.out_md_control)
                    para.input.out_level = "m"; // zhengdy add 2019-04-07

                if (para.input.mdp.md_tfreq == 0)
                {
                    para.input.mdp.md_tfreq = 1.0 / 40 / para.input.mdp.md_dt;
                }
                if (para.input.mdp.md_pfreq == 0)
                {
                    para.input.mdp.md_pfreq = 1.0 / 400 / para.input.mdp.md_dt;
                }
                if (para.input.mdp.md_tfirst < 0 || para.input.mdp.md_restart)
                {
                    para.input.init_vel = 1;
                }
                if (para.input.esolver_type == "lj" || para.input.esolver_type == "dp"
                    || para.input.mdp.md_type == "msst" || para.input.mdp.md_type == "npt")
                {
                    para.input.cal_stress = 1;
                }

                // md_prec_level only used in vc-md  liuyu 2023-03-27
                if (para.input.mdp.md_type != "msst" && para.input.mdp.md_type != "npt")
                {
                    para.input.mdp.md_prec_level = 0;
                }
            }
            else if (calculation == "cell-relax") // mohan add 2011-11-04
            {
                para.input.cal_force = 1;
                para.input.cal_stress = 1;
                if (!para.input.relax_nmax)
                    para.input.relax_nmax = 50;
            }
            else if (calculation == "test_memory" || calculation == "test_neighbour" || calculation == "gen_bessel")
            {
                para.input.relax_nmax = 1;
            }
        };
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            const std::string& calculation = para.input.calculation;
            if (calculation == "nscf" || calculation == "get_S")
            {
                if (para.input.out_dos == 3 && para.input.symmetry == "1")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "symmetry can't be used for out_dos==3(Fermi Surface Plotting) by now.");
                }
            }
            else if (calculation == "get_pchg")
            {
                if (para.input.basis_type == "pw") // xiaohui add 2013-09-01
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "calculate = get_pchg is only availble for LCAO.");
                }
            }
            else if (calculation == "get_wf")
            {
                if (para.input.basis_type == "pw") // xiaohui add 2013-09-01
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "calculate = get_wf is only availble for LCAO.");
                }
            }
            else if (calculation == "md") // mohan add 2011-11-04
            {
                // deal with input parameters , 2019-04-30
                if (para.input.mdp.md_dt < 0)
                    ModuleBase::WARNING_QUIT("ReadInput", "time interval of MD calculation should be set!");
                if (para.input.mdp.md_type == "msst")
                {
                    if (para.input.mdp.msst_qmass <= 0)
                    {
                        ModuleBase::WARNING_QUIT("ReadInput", "msst_qmass must be greater than 0!");
                    }
                }
                if (para.input.esolver_type == "dp")
                {
                    if (access(para.input.mdp.pot_file.c_str(), 0) == -1)
                    {
                        ModuleBase::WARNING_QUIT("ReadInput", "Can not find DP model !");
                    }
                }
            }
            else if (calculation == "gen_bessel")
            {
                if (para.input.basis_type != "pw")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "to generate descriptors, please use pw basis");
                }
            }
            else if (calculation != "scf" && calculation != "relax" && calculation != "cell-relax"
                     && calculation != "test_memory" && calculation != "test_neighbour")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "check 'calculation' !");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("esolver_type");
        item.annotation = "the energy solver: ksdft, sdft, ofdft, tddft, lj, dp";
        read_sync_string(esolver_type);
        this->add_item(item);
    }
    {
        Input_Item item("ntype");
        item.annotation = "atom species number";
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.ntype <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ntype should be greater than 0.");
            }
        };
        read_sync_int(ntype);
        this->add_item(item);
    }
    {
        Input_Item item("nspin");
        item.annotation = "1: single spin; 2: up and down spin; 4: noncollinear spin";
        read_sync_int(nspin);
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nspin != 1 && para.input.nspin != 2 && para.input.nspin != 4)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nspin should be 1, 2 or 4.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("kspacing");
        item.annotation = "unit in 1/bohr, should be > 0, default is 0 which means read KPT file";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            if (count == 1)
            {
                para.input.kspacing[0] = para.input.kspacing[1] = para.input.kspacing[2] = doublevalue;
            }
            else if (count == 3)
            {
                para.input.kspacing[0] = std::stod(item.str_values[0]);
                para.input.kspacing[1] = std::stod(item.str_values[1]);
                para.input.kspacing[2] = std::stod(item.str_values[2]);
            }
            else
            {
                ModuleBase::WARNING_QUIT("ReadInput", "kspacing can only accept one or three double values.");
            }
        };
        sync_doublevec(kspacing, 3);
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            int kspacing_zero_num = 0;
            const std::vector<double>& kspacing = para.input.kspacing;
            for (int i = 0; i < 3; i++)
            {
                if (kspacing[i] < 0.0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "kspacing must > 0");
                }
                else if (kspacing[i] == 0.0)
                {
                    kspacing_zero_num++;
                }
            }
            if (kspacing_zero_num > 0 && kspacing_zero_num < 3)
            {
                std::cout << "kspacing: " << kspacing[0] << " " << kspacing[1] << " " << kspacing[2] << std::endl;
                ModuleBase::WARNING_QUIT("ReadInput", "kspacing must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("min_dist_coef");
        item.annotation = "factor related to the allowed minimum distance between two atoms";
        read_sync_double(min_dist_coef);
        this->add_item(item);
    }
    {
        Input_Item item("nbands");
        item.annotation = "number of bands";
        read_sync_int(nbands);
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nbands < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nbands should be greater than 0.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nbands_istate");
        item.annotation = "number of bands around Fermi level for get_pchg calulation";
        read_sync_int(nbands_istate);
        this->add_item(item);
    }
    {
        Input_Item item("bands_to_print");
        item.annotation = "specify the bands to be calculated in the get_pchg calculation";
        read_sync_string(bands_to_print);
        this->add_item(item);
    }
    {
        Input_Item item("symmetry");
        item.annotation = "the control of symmetry";
        read_sync_string(symmetry);
        autosetfuncs.push_back([](Parameter& para) {
            if (para.input.symmetry == "default")
            {
                if (para.input.gamma_only || para.input.calculation == "nscf" || para.input.calculation == "get_S"
                    || para.input.calculation == "get_pchg" || para.input.calculation == "get_wf")
                    para.input.symmetry = "0"; // if md or exx, symmetry will be force-set to 0 or -1 later
                else
                    para.input.symmetry = "1";
            }
        });
        this->add_item(item);
    }
    {
        Input_Item item("init_vel");
        item.annotation = "read velocity from STRU or not";
        read_sync_bool(init_vel);
        this->add_item(item);
    }
    {
        Input_Item item("symmetry_prec");
        item.annotation = "accuracy for symmetry";
        read_sync_double(symmetry_prec);
        this->add_item(item);
    }
    {
        Input_Item item("symmetry_autoclose");
        item.annotation = "whether to close symmetry automatically when error occurs in symmetry analysis";
        read_sync_bool(symmetry_autoclose);
        this->add_item(item);
    }
    {
        Input_Item item("nelec");
        item.annotation = "input number of electrons";
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nelec < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nelec should be greater than 0.");
            }
            if (para.input.nelec > 0 && para.input.nbands > 0 && para.input.nelec > 2 * para.input.nbands)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nelec > 2*nbnd , bands not enough!");
            }
        };
        read_sync_double(nelec);
        this->add_item(item);
    }
    {
        Input_Item item("nelec_delta");
        item.annotation = "change in the number of total electrons";
        read_sync_double(nelec_delta);
        this->add_item(item);
    }
    {
        Input_Item item("nupdown");
        item.annotation = "the difference number of electrons between spin-up and spin-down";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.input.nupdown = doublevalue;
            para.input.sup.two_fermi = true;
        };

        sync_double(nupdown);
        add_bool_bcast(sup.two_fermi);
        this->add_item(item);
    }
    {
        Input_Item item("out_mul");
        item.annotation = "mulliken charge or not";
        read_sync_bool(out_mul);
        this->add_item(item);
    }
    {
        Input_Item item("noncolin");
        item.annotation = "using non-collinear-spin";
        read_sync_bool(noncolin);
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            if (para.input.noncolin)
            {
                para.input.nspin = 4;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("lspinorb");
        item.annotation = "consider the spin-orbit interaction";
        read_sync_bool(lspinorb);
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            if (para.input.lspinorb)
            {
                para.input.nspin = 4;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("kpar");
        item.annotation = "devide all processors into kpar groups and k points will be distributed among";
        read_sync_int(kpar);
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.basis_type == "lcao")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "kpar > 1 has not been supported for lcao calculation.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("bndpar");
        item.annotation = "devide all processors into bndpar groups and bands will be distributed among each group";
        read_sync_int(bndpar);
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            if (para.input.esolver_type != "sdft")
                para.input.bndpar = 1;
            if (para.input.bndpar > GlobalV::NPROC)
                para.input.bndpar = GlobalV::NPROC;
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_freq_elec");
        item.annotation = "the frequency ( >= 0) of electronic iter to output charge density and wavefunction. 0: "
                          "output only when converged";
        read_sync_int(out_freq_elec);
        this->add_item(item);
    }
    {
        Input_Item item("dft_plus_dmft");
        item.annotation = "true:DFT+DMFT; false: standard DFT calcullation(default)";
        read_sync_bool(dft_plus_dmft);
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.basis_type != "lcao")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "WRONG ARGUMENTS OF basis_type, only lcao is support");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rpa");
        item.annotation = "true:generate output files used in rpa calculation; false:(default)";
        read_sync_bool(rpa);
        this->add_item(item);
    }
    {
        Input_Item item("printe");
        item.annotation = "Print out energy for each band for every printe steps";
        read_sync_int(printe);
        this->add_item(item);
    }
    {
        Input_Item item("mem_saver");
        item.annotation = "Only for nscf calculations. if set to 1, then a memory saving technique will be used for "
                          "many k point calculations.";
        read_sync_int(mem_saver);
        this->add_item(item);
    }
    {
        Input_Item item("diago_proc");
        item.annotation = "the number of procs used to do diagonalization";
        read_sync_int(diago_proc);
        autosetfuncs.push_back([](Parameter& para) {
            if (para.input.diago_proc > GlobalV::NPROC || para.input.diago_proc <= 0)
                para.input.diago_proc = GlobalV::NPROC;
        });

        this->add_item(item);
    }
    {
        Input_Item item("nbspline");
        item.annotation = "the order of B-spline basis";
        read_sync_int(nbspline);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_card");
        item.annotation = "input card for wannier functions";
        read_sync_string(wannier_card);
        this->add_item(item);
    }
    {
        Input_Item item("soc_lambda");
        item.annotation = "The fraction of averaged SOC pseudopotential is given by (1-soc_lambda)";
        read_sync_double(soc_lambda);
        this->add_item(item);
    }
    {
        Input_Item item("cal_force");
        item.annotation = "if calculate the force at the end of the electronic iteration";
        read_sync_bool(cal_force);
        this->add_item(item);
    }
    {
        Input_Item item("out_freq_ion");
        item.annotation = "the frequency ( >= 0 ) of ionic step to output charge density and wavefunction. 0: output "
                          "only when ion steps are finished";
        read_sync_int(out_freq_ion);
        this->add_item(item);
    }
    {
        Input_Item item("elpa_num_thread");
        item.annotation = "Number of threads need to use in elpa";
        read_sync_int(elpa_num_thread);
        this->add_item(item);
    }
    {
        Input_Item item("device");
        item.annotation = "the computing device for ABACUS";
        read_sync_string(device);
        this->add_item(item);
    }
    {
        Input_Item item("precision");
        item.annotation = "the computing precision for ABACUS";
        read_sync_string(precision);
        this->add_item(item);
    }
}

} // namespace ModuleIO