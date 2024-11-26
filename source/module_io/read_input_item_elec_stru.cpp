#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_elec_stru()
{
    // Electronic Structure
    {
        Input_Item item("ks_solver");
        item.annotation = "cg; dav; lapack; genelpa; elpa; scalapack_gvx; cusolver";
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
                "elpa",
                "lapack",
                "scalapack_gvx",
                "cusolver",
                "cusolvermp",
                "pexsi",
                "cg_in_lcao",
            };

            if (para.input.basis_type == "pw")
            {
                if (std::find(pw_solvers.begin(), pw_solvers.end(), ks_solver) == pw_solvers.end())
                {
                    const std::string warningstr = "For PW basis: " + nofound_str(pw_solvers, "ks_solver");
                    ModuleBase::WARNING_QUIT("ReadInput", warningstr);
                }
            }
            else if (para.input.basis_type == "lcao")
            {
                if (std::find(lcao_solvers.begin(), lcao_solvers.end(), ks_solver) == lcao_solvers.end())
                {
                    const std::string warningstr = "For LCAO basis: " + nofound_str(lcao_solvers, "ks_solver");
                    ModuleBase::WARNING_QUIT("ReadInput", warningstr);
                }
                if (ks_solver == "cg_in_lcao")
                {
                    GlobalV::ofs_warning << "cg_in_lcao is under testing" << std::endl;
                }
                else if (ks_solver == "genelpa")
                {
#ifndef __ELPA
                    ModuleBase::WARNING_QUIT("Input",
                                             "Can not use genelpa if abacus is not compiled with "
                                             "ELPA. Please change "
                                             "ks_solver to scalapack_gvx.");
#endif
                }
                else if (ks_solver == "elpa")
                {
#ifndef __ELPA
                    ModuleBase::WARNING_QUIT("Input",
                                             "Can not use elpa if abacus is not compiled with "
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
        Input_Item item("basis_type");
        item.annotation = "PW; LCAO in pw; LCAO";
        read_sync_string(input.basis_type);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.towannier90)
            {
                if (para.input.basis_type == "lcao_in_pw")
                {
                    para.input.basis_type = "lcao";
                }
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            const std::vector<std::string> basis_types = {"pw", "lcao_in_pw", "lcao"};
            if (std::find(basis_types.begin(), basis_types.end(), para.input.basis_type) == basis_types.end())
            {
                const std::string warningstr = nofound_str(basis_types, "basis_type");
                ModuleBase::WARNING_QUIT("ReadInput", warningstr);
            }
        };
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
        Input_Item item("nbands");
        item.annotation = "number of bands";
        read_sync_int(input.nbands);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nbands < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nbands should be greater than 0.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nelec");
        item.annotation = "input number of electrons";
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nelec < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nelec should be greater than 0.");
            }
            if (para.input.nelec > 0 && para.input.nbands > 0 && para.input.nelec > 2 * para.input.nbands)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nelec > 2*nbnd , bands not enough!");
            }
        };
        read_sync_double(input.nelec);
        this->add_item(item);
    }
    {
        Input_Item item("nelec_delta");
        item.annotation = "change in the number of total electrons";
        read_sync_double(input.nelec_delta);
        this->add_item(item);
    }
    {
        Input_Item item("nupdown");
        item.annotation = "the difference number of electrons between spin-up "
                          "and spin-down";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.nupdown = doublevalue;
            para.sys.two_fermi = true;
        };
        item.reset_value = [](const Input_Item&, Parameter& para) {
            if (para.input.nspin == 1)
            {
                para.sys.two_fermi = false;
            }
        };
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.nspin == 1 && para.input.nupdown != 0.0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nupdown mustn't have a non-zero value for spin-unpolarized calculations.");
            }
        };
        sync_double(input.nupdown);
        this->add_item(item);
    }
    {
        Input_Item item("dft_functional");
        item.annotation = "exchange correlation functional";
        read_sync_string(input.dft_functional);
        this->add_item(item);
    }
    {
        Input_Item item("xc_temperature");
        item.annotation = "temperature for finite temperature functionals";
        read_sync_double(input.xc_temperature);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_rcut");
        item.annotation = "default #exchange correlation functional";
        read_sync_double(input.pseudo_rcut);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_mesh");
        item.annotation = "0: use our own mesh to do radial renormalization; "
                          "1: use mesh as in QE";
        read_sync_bool(input.pseudo_mesh);
        this->add_item(item);
    }
    {
        Input_Item item("nspin");
        item.annotation = "1: single spin; 2: up and down spin; 4: noncollinear spin";
        read_sync_int(input.nspin);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.noncolin || para.input.lspinorb)
            {
                para.input.nspin = 4;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nspin != 1 && para.input.nspin != 2 && para.input.nspin != 4)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nspin should be 1, 2 or 4.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_nmax");
        item.annotation = "max iteration number for cg";
        read_sync_int(input.pw_diag_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_thr");
        item.annotation = "threshold for eigenvalues is cg electron iterations";
        read_sync_double(input.pw_diag_thr);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_S" && para.input.basis_type == "pw")
            {
                if (para.input.pw_diag_thr > 1.0e-3)
                {
                    para.input.pw_diag_thr = 1.0e-5;
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_ndim");
        item.annotation = "dimension of workspace for Davidson diagonalization";
        read_sync_int(input.pw_diag_ndim);
        this->add_item(item);
    }
    {
        Input_Item item("diago_cg_prec");
        item.annotation = "diago_cg_prec";
        read_sync_int(input.diago_cg_prec);
        this->add_item(item);
    }
    {
        Input_Item item("smearing_method");
        item.annotation = "type of smearing_method: gauss; fd; fixed; mp; mp2; mv";
        read_sync_string(input.smearing_method);
        this->add_item(item);
    }
    {
        Input_Item item("smearing_sigma");
        item.annotation = "energy range for smearing";
        read_sync_double(input.smearing_sigma);
        this->add_item(item);
    }
    {
        // Energy range for smearing,
        //`smearing_sigma` = 1/2 *kB* `smearing_sigma_temp`.
        Input_Item tmp_item("smearing_sigma_temp");
        tmp_item.read_value
            = [](const Input_Item& item, Parameter& para) { para.input.smearing_sigma = 3.166815e-6 * doublevalue; };
        // only to set smearing_sigma, so no need to write to output INPUT file
        // or bcast.
        this->add_item(tmp_item);
    }
    {
        Input_Item item("mixing_type");
        item.annotation = "plain; pulay; broyden";
        read_sync_string(input.mixing_mode);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_beta");
        item.annotation = "mixing parameter: 0 means no new charge";
        read_sync_double(input.mixing_beta);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.mixing_beta < 0.0)
            {
                if (para.input.nspin == 1)
                {
                    para.input.mixing_beta = 0.8;
                }
                else if (para.input.nspin == 2)
                {
                    para.input.mixing_beta = 0.4;
                    para.input.mixing_beta_mag = 1.6;
                    para.input.mixing_gg0_mag = 0.0;
                }
                else if (para.input.nspin == 4) // I will add this
                {
                    para.input.mixing_beta = 0.4;
                    para.input.mixing_beta_mag = 1.6;
                    para.input.mixing_gg0_mag = 0.0;
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("mixing_ndim");
        item.annotation = "mixing dimension in pulay or broyden";
        read_sync_int(input.mixing_ndim);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_restart");
        item.annotation = "threshold to restart mixing during SCF";
        read_sync_double(input.mixing_restart);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.sc_mag_switch == 1)
            {// for DeltaSpin calculation, the mixing_restart should be same as sc_scf_thr
                if(para.input.sc_scf_thr != 10.0)
                {
                    para.input.mixing_restart = para.input.sc_scf_thr;
                }
                else
                {// no mixing_restart until oscillation happen in PW base
                    para.input.mixing_restart = para.input.scf_thr / 10.0;
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0");
        item.annotation = "mixing parameter in kerker";
        read_sync_double(input.mixing_gg0);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_beta_mag");
        item.annotation = "mixing parameter for magnetic density";
        read_sync_double(input.mixing_beta_mag);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.mixing_beta_mag < 0.0)
            {
                if (para.input.nspin == 2 || para.input.nspin == 4)
                {
                    if (para.input.mixing_beta <= 0.4)
                    {
                        para.input.mixing_beta_mag = 4 * para.input.mixing_beta;
                    }
                    else
                    {
                        para.input.mixing_beta_mag = 1.6; // 1.6 can be discussed
                    }
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0_mag");
        item.annotation = "mixing parameter in kerker";
        read_sync_double(input.mixing_gg0_mag);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0_min");
        item.annotation = "the minimum kerker coefficient";
        read_sync_double(input.mixing_gg0_min);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_angle");
        item.annotation = "angle mixing parameter for non-colinear calculations";
        read_sync_double(input.mixing_angle);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_tau");
        item.annotation = "whether to mix tau in mGGA calculation";
        read_sync_bool(input.mixing_tau);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_dftu");
        item.annotation = "whether to mix locale in DFT+U calculation";
        read_sync_bool(input.mixing_dftu);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_dmr");
        item.annotation = "whether to mix real-space density matrix";
        read_sync_bool(input.mixing_dmr);
        this->add_item(item);
    }
    {
        Input_Item item("gamma_only");
        item.annotation = "Only for localized orbitals set and gamma point. If "
                          "set to 1, a fast algorithm is used";
        read_sync_bool(input.gamma_only);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.basis_type == "pw" && para.input.gamma_only) 
            {
                para.input.gamma_only = false;   
                GlobalV::ofs_warning << " WARNING : gamma_only has not been implemented for pw yet" << std::endl;
                GlobalV::ofs_warning << "gamma_only is not supported in the pw model" << std::endl;
                GlobalV::ofs_warning << " the INPUT parameter gamma_only has been reset to 0" << std::endl;
                GlobalV::ofs_warning << " and a new KPT is generated with gamma point as the only k point"<< std::endl;
                GlobalV::ofs_warning << " Auto generating k-points file: " << para.input.kpoint_file << std::endl;
                std::ofstream ofs(para.input.kpoint_file.c_str());
                ofs << "K_POINTS" << std::endl;
                ofs << "0" << std::endl;
                ofs << "Gamma" << std::endl;
                ofs << "1 1 1 0 0 0" << std::endl;
                ofs.close();
            }
            if (para.input.basis_type == "lcao" && para.input.gamma_only)
            {
                if (para.input.nspin == 4)
                {
                    ModuleBase::WARNING_QUIT("NOTICE", "nspin=4(soc or noncollinear-spin) does not support gamma only calculation");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("scf_nmax");
        item.annotation = "number of electron iterations";
        read_sync_int(input.scf_nmax);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "nscf")
            {
                para.input.scf_nmax = 1;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("scf_thr");
        item.annotation = "charge density error";
        read_sync_double(input.scf_thr);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.scf_thr == -1.0)
            {
                if (para.input.basis_type == "lcao" || para.input.basis_type == "lcao_in_pw")
                {
                    para.input.scf_thr = 1.0e-7;
                }
                else if (para.input.basis_type == "pw" && para.input.calculation != "nscf")
                {
                    para.input.scf_thr = 1.0e-9;
                }
                else if (para.input.basis_type == "pw" && para.input.calculation == "nscf")
                {
                    para.input.scf_thr = 1.0e-6;
                    // In NSCF calculation, the diagonalization threshold is set
                    // to 0.1*scf/nelec. In other words, the scf_thr is used to
                    // control diagonalization convergence threthod in NSCF. In
                    // this case, the default 1.0e-9 is too strict. renxi
                    // 20230908
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("scf_ene_thr");
        item.annotation = "total energy error threshold";
        read_sync_double(input.scf_ene_thr);
        this->add_item(item);
    }
    {
        Input_Item item("scf_os_stop");
        item.annotation = "whether to stop scf when oscillation is detected";
        read_sync_bool(input.scf_os_stop);
        this->add_item(item);
    }
    {
        Input_Item item("scf_os_thr");
        item.annotation = "charge density threshold for oscillation";
        read_sync_double(input.scf_os_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.scf_os_thr >= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "scf_os_thr should be negative");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("scf_os_ndim");
        item.annotation = "number of old iterations used for oscillation detection";
        read_sync_int(input.scf_os_ndim);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.scf_os_ndim <= 0) // default value
            {
                para.input.scf_os_ndim = para.input.mixing_ndim;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("scf_thr_type");
        item.annotation = "type of the criterion of scf_thr, 1: reci drho for "
                          "pw, 2: real drho for lcao";
        read_sync_int(input.scf_thr_type);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.scf_thr_type == -1)
            {
                if (para.input.basis_type == "lcao" || para.input.basis_type == "lcao_in_pw")
                {
                    para.input.scf_thr_type = 2;
                }
                else if (para.input.basis_type == "pw")
                {
                    para.input.scf_thr_type = 1;
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("lspinorb");
        item.annotation = "consider the spin-orbit interaction";
        read_sync_bool(input.lspinorb);
        this->add_item(item);
    }
    {
        Input_Item item("noncolin");
        item.annotation = "using non-collinear-spin";
        read_sync_bool(input.noncolin);
        this->add_item(item);
    }
    {
        Input_Item item("soc_lambda");
        item.annotation = "The fraction of averaged SOC pseudopotential is "
                          "given by (1-soc_lambda)";
        read_sync_double(input.soc_lambda);
        this->add_item(item);
    }

    // LCAO
    {
        Input_Item item("nb2d");
        item.annotation = "matrix 2d division";
        read_sync_int(input.nb2d);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nb2d < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nb2d should be greater than 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("lmaxmax");
        item.annotation = "maximum of l channels used";
        read_sync_int(input.lmaxmax);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_ecut");
        item.annotation = "energy cutoff for LCAO";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.lcao_ecut == 0 && para.input.basis_type == "lcao")
            {
                para.input.lcao_ecut = para.input.ecutwfc;
                ModuleBase::GlobalFunc::AUTO_SET("lcao_ecut", para.input.ecutwfc);
            }
        };
        read_sync_double(input.lcao_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_dk");
        item.annotation = "delta k for 1D integration in LCAO";
        read_sync_double(input.lcao_dk);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_dr");
        item.annotation = "delta r for 1D integration in LCAO";
        read_sync_double(input.lcao_dr);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_rmax");
        item.annotation = "max R for 1D two-center integration table";
        read_sync_double(input.lcao_rmax);
        this->add_item(item);
    }
    {
        Input_Item item("search_radius");
        item.annotation = "input search radius (Bohr)";
        read_sync_double(input.search_radius);
        this->add_item(item);
    }
    {
        Input_Item item("search_pbc");
        item.annotation = "input periodic boundary condition";
        read_sync_bool(input.search_pbc);
        this->add_item(item);
    }
    {
        Input_Item item("bx");
        item.annotation = "division of an element grid in FFT grid along x";
        read_sync_int(input.bx);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.bx > 10)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "bx should be no more than 10");
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.basis_type == "pw" || para.input.basis_type == "lcao_in_pw"
                || para.input.calculation == "get_wf")
            {
                para.input.bx = 1;
                para.input.by = 1;
                para.input.bz = 1;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("by");
        item.annotation = "division of an element grid in FFT grid along y";
        read_sync_int(input.by);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.by > 10)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "by should be no more than 10");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("bz");
        item.annotation = "division of an element grid in FFT grid along z";
        read_sync_int(input.bz);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.bz > 10)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "bz should be no more than 10");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("elpa_num_thread");
        item.annotation = "Number of threads need to use in elpa";
        read_sync_int(input.elpa_num_thread);
        this->add_item(item);
    }
    {
        Input_Item item("num_stream");
        item.annotation = "the nstream in compute the LCAO with CUDA";
        read_sync_int(input.nstream);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_ecut");
        item.annotation = "energy cutoff for spherical bessel functions(Ry)";
        read_sync_string(input.bessel_nao_ecut);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.bessel_nao_ecut == "default")
            {
                para.input.bessel_nao_ecut = std::to_string(para.input.ecutwfc);
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (std::stod(para.input.bessel_nao_ecut) < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "bessel_nao_ecut must >= 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_tolerence");
        item.annotation = "tolerence for spherical bessel root";
        read_sync_double(input.bessel_nao_tolerence);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_rcut");
        item.annotation = "radial cutoff for spherical bessel functions(a.u.)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.bessel_nao_rcuts.push_back(std::stod(item.str_values[i]));
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            for(auto rcut: para.input.bessel_nao_rcuts)
            {
                if (rcut < 0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "bessel_nao_rcut must >= 0");
                }
            }
        };
        sync_doublevec(input.bessel_nao_rcuts, para.input.bessel_nao_rcuts.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_smooth");
        item.annotation = "spherical bessel smooth or not";
        read_sync_bool(input.bessel_nao_smooth);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_sigma");
        item.annotation = "spherical bessel smearing_sigma";
        read_sync_double(input.bessel_nao_sigma);
        this->add_item(item);
    }
}
} // namespace ModuleIO
