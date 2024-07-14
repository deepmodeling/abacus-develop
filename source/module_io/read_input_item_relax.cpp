#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_relax()
{
    {
        Input_Item item("ks_solver");
        item.annotation = "cg; dav; lapack; genelpa; scalapack_gvx; cusolver";
        read_sync_string(input.ks_solver);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.ks_solver == "default")
            {
                if (para.input.basis_type == "pw")
                {
                    para.input.ks_solver = "cg";
                    ModuleBase::GlobalFunc::AUTO_SET("ks_solver", "cg");
                }
                else if (para.input.basis_type == "lcao")
                {
                    if (para.input.device == "gpu")
                    {
                        para.input.ks_solver = "cusolver";
                        ModuleBase::GlobalFunc::AUTO_SET("ks_solver", "cusolver");
                    }
                    else
                    {
#ifdef __ELPA
                        para.input.ks_solver = "genelpa";
                        ModuleBase::GlobalFunc::AUTO_SET("ks_solver", "genelpa");
#else
#ifdef __MPI
                        para.input.ks_solver = "scalapack_gvx";
                        ModuleBase::GlobalFunc::AUTO_SET("ks_solver", "scalapack_gvx");
#else
                        para.input.ks_solver = "lapack";
                        ModuleBase::GlobalFunc::AUTO_SET("ks_solver", "lapack");
#endif
#endif
                    }
                }
            }
            if (para.input.towannier90)
            {
                if (para.input.basis_type == "lcao_in_pw")
                {
#ifdef __ELPA
                    para.input.ks_solver = "genelpa";
#else
#ifdef __MPI
                    para.input.ks_solver = "scalapack_gvx";
#else
                    para.input.ks_solver = "lapack";
#endif
#endif
                }
            };
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            const std::string& ks_solver = para.input.ks_solver;
            const std::vector<std::string> pw_solvers = {"cg", "dav", "bpcg", "dav_subspace"};
            const std::vector<std::string> lcao_solvers = {
                "genelpa",
                "lapack",
                "scalapack_gvx",
                "cusolver",
                "pexsi",
                "cg_in_lcao",
            };

            if (para.input.basis_type == "pw")
            {
                if (!find_str(pw_solvers, ks_solver))
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "ks_solver must be cg, dav, bpcg "
                                             "or dav_subspace for pw basis.");
                }
            }
            else if (para.input.basis_type == "lcao")
            {
                if (!find_str(lcao_solvers, ks_solver))
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "ks_solver must be genelpa, lapack, scalapack_gvx, "
                                             "cusolver, pexsi or "
                                             "cg_in_lcao for lcao basis.");
                }
                if (ks_solver == "cg_in_lcao")
                {
                    GlobalV::ofs_warning << "cg_in_lcao is under testing" << std::endl;
                }
                else if (ks_solver == "genelpa")
                {
#ifndef __MPI
                    ModuleBase::WARNING_QUIT("ReadInput", "genelpa can not be used for series version.");
#endif
#ifndef __ELPA
                    ModuleBase::WARNING_QUIT("Input",
                                             "Can not use genelpa if abacus is not compiled with "
                                             "ELPA. Please change "
                                             "ks_solver to scalapack_gvx.");
#endif
                }
                else if (ks_solver == "scalapack_gvx")
                {
#ifdef __MPI
                    GlobalV::ofs_warning << "scalapack_gvx is under testing" << std::endl;
#else
                    ModuleBase::WARNING_QUIT("ReadInput", "scalapack_gvx can not be used for series version.");
#endif
                }
                else if (ks_solver == "cusolver" || ks_solver == "cusolvermp")
                {
#ifndef __MPI
                    ModuleBase::WARNING_QUIT("ReadInput", "Cusolver can not be used for series version.");
#endif
                }
                else if (ks_solver == "pexsi")
                {
#ifdef __PEXSI
                    GlobalV::ofs_warning << " It's ok to use pexsi." << std::endl;
#else
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "Can not use PEXSI if abacus is not compiled with "
                                             "PEXSI. Please change "
                                             "ks_solver to scalapack_gvx.");
#endif
                }
            }
            else if (para.input.basis_type == "lcao_in_pw")
            {
                if (ks_solver != "lapack")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "LCAO in plane wave can only done with lapack.");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("scf_nmax");
        item.annotation = "number of electron iterations";
        read_sync_int(input.scf_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("relax_nmax");
        item.annotation = "number of ion iteration steps";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            const std::string& calculation = para.input.calculation;
            const std::vector<std::string> singlelist
                = {"scf", "nscf", "get_S", "get_pchg", "get_wf", "test_memory", "test_neighbour", "gen_bessel"};
            if (find_str(singlelist, calculation))
            {
                para.input.relax_nmax = 1;
            }
            else if (calculation == "relax" || calculation == "cell-relax")
            {
                if (!para.input.relax_nmax)
                {
                    para.input.relax_nmax = 50;
                }
            }
        };
        read_sync_int(input.relax_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("out_stru");
        item.annotation = "output the structure files after each ion step";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            const std::vector<std::string> offlist = {"nscf", "get_S", "get_pchg", "get_wf"};
            if (find_str(offlist, para.input.calculation))
            {
                para.input.out_stru = false;
            }
        };
        read_sync_bool(input.out_stru);
        this->add_item(item);
    }
    {
        Input_Item item("force_thr");
        item.annotation = "force threshold, unit: Ry/Bohr";
        // read_sync_double(input.force_thr);
        item.read_value = [](const Input_Item& item, Parameter& para) { para.input.force_thr = doublevalue; };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.force_thr == -1 && para.input.force_thr_ev == -1)
            {
                para.input.force_thr = 1.0e-3; // default value
                para.input.force_thr_ev = para.input.force_thr * 13.6058 / 0.529177;
            }
            else if (para.input.force_thr == -1 && para.input.force_thr_ev != -1)
            {
                para.input.force_thr = para.input.force_thr_ev / 13.6058 * 0.529177;
            }
            else
            {
                // if both force_thr and force_thr_ev are set, use force_thr
                ModuleBase::WARNING("ReadInput", "both force_thr and force_thr_ev are set, use force_thr");
                para.input.force_thr_ev = para.input.force_thr * 13.6058 / 0.529177;
            }
            para.input.mdp.force_thr = para.input.force_thr; // temperaory
        };
        sync_double(input.force_thr);
        add_double_bcast(input.mdp.force_thr);
        this->add_item(item);
    }
    {
        Input_Item item("force_thr_ev");
        item.annotation = "force threshold, unit: eV/Angstrom";
        item.read_value = [](const Input_Item& item, Parameter& para) { para.input.force_thr_ev = doublevalue; };
        sync_double(input.force_thr_ev);
        this->add_item(item);
    }
    {
        Input_Item item("force_thr_ev2");
        item.annotation = "force invalid threshold, unit: eV/Angstrom";
        read_sync_double(input.force_thr_ev2);
        this->add_item(item);
    }
    {
        Input_Item item("relax_cg_thr");
        item.annotation = "threshold for switching from cg to bfgs, unit: eV/Angstrom";
        read_sync_double(input.relax_cg_thr);
        this->add_item(item);
    }
    {
        Input_Item item("stress_thr");
        item.annotation = "stress threshold";
        read_sync_double(input.stress_thr);
        this->add_item(item);
    }
    {
        Input_Item item("press1");
        item.annotation = "target pressure, unit: KBar";
        read_sync_double(input.press1);
        this->add_item(item);
    }
    {
        Input_Item item("press2");
        item.annotation = "target pressure, unit: KBar";
        read_sync_double(input.press2);
        this->add_item(item);
    }
    {
        Input_Item item("press3");
        item.annotation = "target pressure, unit: KBar";
        read_sync_double(input.press3);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_w1");
        item.annotation = "wolfe condition 1 for bfgs";
        read_sync_double(input.relax_bfgs_w1);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_w2");
        item.annotation = "wolfe condition 2 for bfgs";
        read_sync_double(input.relax_bfgs_w2);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_rmax");
        item.annotation = "maximal trust radius, unit: Bohr";
        read_sync_double(input.relax_bfgs_rmax);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_rmin");
        item.annotation = "minimal trust radius, unit: Bohr";
        read_sync_double(input.relax_bfgs_rmin);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_init");
        item.annotation = "initial trust radius, unit: Bohr";
        read_sync_double(input.relax_bfgs_init);
        this->add_item(item);
    }
    {
        Input_Item item("cal_stress");
        item.annotation = "calculate the stress or not";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "md")
            {
                if (para.input.esolver_type == "lj" || para.input.esolver_type == "dp"
                    || para.input.mdp.md_type == "msst" || para.input.mdp.md_type == "npt")
                {
                    para.input.cal_stress = true;
                }
            }
            else if (para.input.calculation == "cell-relax")
            {
                para.input.cal_stress = true;
            }
            para.input.mdp.cal_stress = para.input.cal_stress; // temperaory
        };
        read_sync_bool(input.cal_stress);
        add_bool_bcast(input.mdp.cal_stress);
        this->add_item(item);
    }
    {
        Input_Item item("fixed_axes");
        item.annotation = "which axes are fixed";
        read_sync_string(input.fixed_axes);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if ((para.input.fixed_axes == "shape" || para.input.fixed_axes == "volume") && !para.input.relax_new)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "fixed shape and fixed volume only supported for relax_new = 1");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("fixed_ibrav");
        item.annotation = "whether to preseve lattice type during relaxation";
        read_sync_bool(input.fixed_ibrav);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.fixed_ibrav && !para.input.relax_new)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "fixed_ibrav only available for relax_new = 1");
            }
            if (para.input.latname == "none" && para.input.fixed_ibrav)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "to use fixed_ibrav, latname must be provided");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("fixed_atoms");
        item.annotation = "whether to preseve direct coordinates of atoms "
                          "during relaxation";
        read_sync_bool(input.fixed_atoms);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.fixed_atoms && para.input.calculation == "relax")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "fixed_atoms is not meant to be used for calculation = relax");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("relax_method");
        item.annotation = "cg; bfgs; sd; cg; cg_bfgs;";
        read_sync_string(input.relax_method);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            const std::vector<std::string> relax_methods = {"cg", "bfgs", "sd", "cg_bfgs"};
            if (!find_str(relax_methods, para.input.relax_method))
            {
                ModuleBase::WARNING_QUIT("ReadInput", "relax_method must be cg, bfgs, sd or cg_bfgs.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("relax_new");
        item.annotation = "whether to use the new relaxation method";
        read_sync_bool(input.relax_new);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.relax_new && para.input.relax_method != "cg")
            {
                para.input.relax_new = false;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("relax_scale_force");
        item.annotation = "controls the size of the first CG step if relax_new is true";
        read_sync_double(input.relax_scale_force);
        this->add_item(item);
    }
    {
        Input_Item item("out_level");
        item.annotation = "ie(for electrons); i(for ions);";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.out_level = strvalue;
            para.sys.out_md_control = true;
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (!para.sys.out_md_control && para.input.calculation == "md")
            {
                para.input.out_level = "m"; // zhengdy add 2019-04-07
            }
        };
        sync_string(input.out_level);
        this->add_item(item);
    }
    {
        Input_Item item("out_dm");
        item.annotation = ">0 output density matrix";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_pchg" || para.input.calculation == "get_wf")
            {
                para.input.out_dm = false;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.sys.gamma_only_local == false && para.input.out_dm)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_dm with k-point algorithm is not implemented yet.");
            }
        };
        read_sync_bool(input.out_dm);
        this->add_item(item);
    }
    {
        Input_Item item("out_dm1");
        item.annotation = ">0 output density matrix (multi-k points)";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_pchg" || para.input.calculation == "get_wf")
            {
                para.input.out_dm1 = false;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.sys.gamma_only_local == true && para.input.out_dm1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_dm1 is only for multi-k");
            }
        };
        read_sync_bool(input.out_dm1);
        this->add_item(item);
    }
    {
        Input_Item item("out_bandgap");
        item.annotation = "if true, print out bandgap";
        read_sync_bool(input.out_bandgap);
        this->add_item(item);
    }
    {
        Input_Item item("use_paw");
        item.annotation = "whether to use PAW in pw calculation";
        read_sync_bool(input.use_paw);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.use_paw)
            {
#ifndef USE_PAW
                ModuleBase::WARNING_QUIT("ReadInput", "to use PAW, compile with USE_PAW");
#endif
                if (para.input.basis_type != "pw")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "PAW is for pw basis only");
                }
                if (para.input.dft_functional == "default")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "dft_functional must be set when use_paw is true");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("deepks_out_labels");
        item.annotation = ">0 compute descriptor for deepks";
        read_sync_bool(input.deepks_out_labels);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_scf");
        item.annotation = ">0 add V_delta to Hamiltonian";
        read_sync_bool(input.deepks_scf);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_equiv");
        item.annotation = "whether to use equivariant version of DeePKS";
        read_sync_bool(input.deepks_equiv);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_bandgap");
        item.annotation = ">0 for bandgap label";
        read_sync_bool(input.deepks_bandgap);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_v_delta");
        item.annotation = ">0 for v_delta label. when output, 1 for v_delta_precalc, 2 for psialpha and grad_evdm ( can save memory )";
        read_sync_int(input.deepks_v_delta);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_out_unittest");
        item.annotation = "if set 1, prints intermediate quantities that shall "
                          "be used for making unit test";
        read_sync_bool(input.deepks_out_unittest);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_model");
        item.annotation = "file dir of traced pytorch model: 'model.ptg";
        read_sync_string(input.deepks_model);
        this->add_item(item);
    }
}
} // namespace ModuleIO