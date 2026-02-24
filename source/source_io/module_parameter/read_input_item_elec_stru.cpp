#include "source_base/global_function.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_elec_stru()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("basis_type");
        item.annotation = "PW; LCAO in pw; LCAO";
        item.category = "Electronic structure";
        item.type = "String";
        item.description = R"(Choose the basis set.
* pw: Using plane-wave basis set only.
* lcao: Using localized atomic orbital sets.
* lcao_in_pw: Expand the localized atomic set in plane-wave basis, non-self-consistent field calculation not tested.)";
        item.default_value = "pw";
        item.unit = "";
        item.availability = "";
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
    // Electronic Structure
    {
        Input_Item item("ks_solver");
        item.annotation = "cg; dav; lapack; genelpa; elpa; scalapack_gvx; cusolver";
        item.category = "Electronic structure";
        item.type = "String";
        item.description = R"(Choose the diagonalization methods for the Hamiltonian matrix expanded in a certain basis set.

For plane-wave basis,

* cg: The conjugate-gradient (CG) method.
* bpcg: The BPCG method, which is a block-parallel Conjugate Gradient (CG) method, typically exhibits higher acceleration in a GPU environment.
* dav: The Davidson algorithm.
* dav_subspace: The Davidson algorithm without orthogonalization operation, this method is the most recommended for efficiency. pw_diag_ndim can be set to 2 for this method.

For numerical atomic orbitals basis,

* lapack: Use LAPACK to diagonalize the Hamiltonian, only used for serial version
* genelpa: Use GEN-ELPA to diagonalize the Hamiltonian.
* scalapack_gvx: Use Scalapack to diagonalize the Hamiltonian.
* cusolver: Use CUSOLVER to diagonalize the Hamiltonian, at least one GPU is needed.
* cusolvermp: Use CUSOLVER to diagonalize the Hamiltonian, supporting multi-GPU devices. Note that you should set the number of MPI processes equal to the number of GPUs.
* elpa: The ELPA solver supports both CPU and GPU. By setting the device to GPU, you can launch the ELPA solver with GPU acceleration (provided that you have installed a GPU-supported version of ELPA, which requires you to manually compile and install ELPA, and the ABACUS should be compiled with -DUSE_ELPA=ON and -DUSE_CUDA=ON). The ELPA solver also supports multi-GPU acceleration.

If you set ks_solver=genelpa for basis_type=pw, the program will stop with an error message:

``text genelpa can not be used with plane wave basis. ``

Then the user has to correct the input file and restart the calculation.)";
        item.default_value = "";
        item.unit = "";
        item.availability = "";
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
                    std::string warningstr;
#ifndef __MPI
                    ModuleBase::WARNING_QUIT("ReadInput", "Cusolver can not be used for series version.");
#endif
#ifndef __CUDA
                    warningstr = "ks_solver is set to " + ks_solver + " but ABACUS is built with CPU only!\n"
                    + " Please rebuild ABACUS with GPU support or change the ks_solver.";
                    ModuleBase::WARNING_QUIT("ReadInput", warningstr);
#endif
                    if( ks_solver == "cusolvermp")
                    {
#ifndef __CUSOLVERMP
                    warningstr = "ks_solver is set to cusolvermp, but ABACUS is not built with cusolvermp support\n"
                    " Please rebuild ABACUS with cusolvermp support or change the ks_solver.";
                    ModuleBase::WARNING_QUIT("ReadInput", warningstr);
#endif
                    }
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
        Input_Item item("nbands");
        item.annotation = "number of bands";
        item.category = "Electronic structure";
        item.type = "Integer";
        item.description = "The number of Kohn-Sham orbitals to calculate. It is recommended to setup this value, especially when smearing techniques are utilized, more bands should be included.";
        item.default_value = "";
        item.unit = "";
        item.availability = "";
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
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = R"(* 0.0: The total number of electrons will be calculated by the sum of valence electrons (i.e. assuming neutral system).
* >0.0: this denotes the total number of electrons in the system. Must be less than 2*nbands.)";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "";
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
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "The total number of electrons will be calculated by nelec+nelec_delta.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.nelec_delta);
        this->add_item(item);
    }
    {
        Input_Item item("nupdown");
        item.annotation = "the difference number of electrons between spin-up "
                          "and spin-down";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = R"(* 0.0: no constrain apply to system.
* >0.0: The different number of electrons between spin-up and spin-down channels. The range of value must be in [-nelec ~ nelec]. It is one type of constrainted DFT method, two Fermi energies will be calculated.)";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "";
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
        add_bool_bcast(sys.two_fermi);
        this->add_item(item);
    }
    {
        Input_Item item("dft_functional");
        item.annotation = "exchange correlation functional";
        item.category = "Electronic structure";
        item.type = "String";
        item.description = R"(In our package, the XC functional can either be set explicitly using the dft_functional keyword in INPUT file. If dft_functional is not specified, ABACUS will use the xc functional indicated in the pseudopotential file. On the other hand, if dft_functional is specified, it will overwrite the functional from pseudopotentials and performs calculation with whichever functional the user prefers. We further offer two ways of supplying exchange-correlation functional. The first is using 'short-hand' names. A complete list of 'short-hand' expressions can be found in the source code. Supported density functionals are:
* LDA functionals
  * LDA (equivalent with PZ and SLAPZNOGXNOGC), PWLDA
* GGA functionals
  * PBE (equivalent with SLAPWPBXPBC), PBESOL, REVPBE, WC, BLYP, BP(referred to BP86), PW91, HCTH, OLYP, BLYP_LR
* meta-GGA functionals
  * SCAN (require LIBXC)
* Hybrid functionals
  * PBE0, HF
  * If LIBXC is available, additional short-hand names of hybrid functionals are supported: HSE(referred to HSE06), B3LYP, LC_PBE, LC_WPBE, LRC_WPBE, LRC_WPBEH, CAM_PBEH, WP22, CWP22, MULLER (equivalent with POWER)
* Hybrid meta-GGA functionals
  * SCAN0 (require LIBXC)

The other way is only available when compiling with LIBXC, and it allows for supplying exchange-correlation functionals as combinations of LIBXC keywords for functional components, joined by a plus sign, for example, dft_functional='LDA_X_1D_EXPONENTIAL+LDA_C_1D_CSC'.)";
        item.default_value = "Used the same as DFT functional as specified in the pseudopotential files.";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.dft_functional);
        this->add_item(item);
    }
    {
        Input_Item item("xc_temperature");
        item.annotation = "temperature for finite temperature functionals";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "Specifies temperature when using temperature-dependent XC functionals (KSDT and so on).";
        item.default_value = "0.0";
        item.unit = "Ry";
        item.availability = "";
        read_sync_double(input.xc_temperature);
        this->add_item(item);
    }
    {
        Input_Item item("xc_exch_ext");
        item.annotation = "placeholder for xcpnet exchange functional";
        item.category = "Electronic structure";
        item.type = "Integer Real ...";
        item.description = "Customized parameterization on the exchange part of XC functional. The first value should be the LibXC ID of the original functional, and latter values are external parameters. Default values are those of Perdew-Burke-Ernzerhof (PBE) functional. For more information on LibXC ID of functionals, please refer to LibXC. For parameters of functionals of interest, please refer to the source code of LibXC, such as PBE functional interface in LibXC: gga_x_pbe.c."
                          "\n\n[NOTE] Solely setting this keyword will take no effect on XC functionals. One should also set "
                          "dft_functional to the corresponding functional to apply the customized parameterization. "
                          "Presently this feature can only support parameterization on one exchange functional.";
        item.default_value = "101 0.8040 0.2195149727645171";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.xc_exch_ext.resize(item.get_size());
            std::transform(item.str_values.begin(), item.str_values.end(),
                           para.input.xc_exch_ext.begin(),
                           [](const std::string& str) { return std::stod(str); });
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            // at least one value should be set
            if (para.input.xc_exch_ext.empty())
            {
                ModuleBase::WARNING_QUIT("ReadInput", "xc_exch_ext should not be empty.");
            }
            // the first value is actually an integer, not a double
            const double libxc_id_dbl = para.input.xc_exch_ext[0];
            if (std::abs(libxc_id_dbl - std::round(libxc_id_dbl)) > 1.0e-6)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                    "The first parameter (libxc id) can never be a float number");
            }
            // the first value is a positive integer
            if (libxc_id_dbl < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                    "The first parameter (libxc id) should be a positive integer");
            }
        };
        sync_doublevec(input.xc_exch_ext,
                       para.input.xc_exch_ext.size(),
                       0.0);
        this->add_item(item);
    }
    {
        Input_Item item("xc_corr_ext");
        item.annotation = "placeholder for xcpnet exchange functional";
        item.category = "Electronic structure";
        item.type = "Integer Real ...";
        item.description = "Customized parameterization on the correlation part of XC functional. The first value should be the LibXC ID of the original functional, and latter values are external parameters. Default values are those of Perdew-Burke-Ernzerhof (PBE) functional. For more information on LibXC ID of functionals, please refer to LibXC. For parameters of functionals of interest, please refer to the source code of LibXC, such as PBE functional interface in LibXC: gga_c_pbe.c."
                          "\n\n[NOTE] Solely setting this keyword will take no effect on XC functionals. One should also set "
                          "dft_functional to the corresponding functional to apply the customized parameterization. "
                          "Presently this feature can only support parameterization on one correlation functional.";
        item.default_value = "130 0.06672455060314922 0.031090690869654895034 1.0";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.xc_corr_ext.resize(item.get_size());
            std::transform(item.str_values.begin(), item.str_values.end(),
                           para.input.xc_corr_ext.begin(),
                           [](const std::string& str) { return std::stod(str); });
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            // at least one value should be set
            if (para.input.xc_corr_ext.empty())
            {
                ModuleBase::WARNING_QUIT("ReadInput", "xc_corr_ext should not be empty.");
            }
            // the first value is actually an integer, not a double
            const double libxc_id_dbl = para.input.xc_corr_ext[0];
            if (std::abs(libxc_id_dbl - std::round(libxc_id_dbl)) > 1.0e-6)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                    "The first parameter (libxc id) can never be a float number");
            }
            // the first value is a positive integer
            if (libxc_id_dbl < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                    "The first parameter (libxc id) should be a positive integer");
            }
        };
        sync_doublevec(input.xc_corr_ext,
                       para.input.xc_corr_ext.size(),
                       0.0);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_rcut");
        item.annotation = "default #exchange correlation functional";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "Cut-off of radial integration for pseudopotentials.";
        item.default_value = "15";
        item.unit = "Bohr";
        item.availability = "";
        read_sync_double(input.pseudo_rcut);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_mesh");
        item.annotation = "0: use our own mesh to do radial renormalization; "
                          "1: use mesh as in QE";
        item.category = "Electronic structure";
        item.type = "Boolean";
        item.description = R"(* 0: Use a mesh for radial integration of pseudopotentials.
* 1: Use the mesh that is consistent with quantum espresso)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pseudo_mesh);
        this->add_item(item);
    }
    {
        Input_Item item("nspin");
        item.annotation = "1: single spin; 2: up and down spin; 4: noncollinear spin";
        item.category = "Electronic structure";
        item.type = "Integer";
        item.description = R"(The number of spin components of wave functions.
* 1: Spin degeneracy
* 2: Collinear spin polarized.
* 4: For the case of noncollinear polarized, nspin will be automatically set to 4 without being specified by the user.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("smearing_method");
        item.annotation = "type of smearing_method: gauss; fd; fixed; mp; mp2; mv";
        item.category = "Electronic structure";
        item.type = "String";
        item.description = R"(It indicates which occupation and smearing method is used in the calculation.
* fixed: fixed occupations (available for non-coductors only)
* gauss or gaussian: Gaussian smearing method.
* mp: methfessel-paxton smearing method; recommended for metals.
* mp2: 2-nd methfessel-paxton smearing method; recommended for metals.
* mv or cold: marzari-vanderbilt smearing method.
* fd: Fermi-Dirac smearing method: and smearing_sigma below is the temperature (in Ry).)";
        item.default_value = "gauss";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.smearing_method);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            const std::vector<std::string> methods = {"gauss", "gaussian",
                                                      "fd", "fermi-dirac",
                                                      "fixed",
                                                      "mp", "mp2", "mp3"
                                                      "marzari-vanderbilt", "cold", "mv"};
            if (std::find(methods.begin(), methods.end(), para.input.smearing_method) == methods.end())
            {
                const std::string warningstr = nofound_str(methods, "smearing_method");
                ModuleBase::WARNING_QUIT("ReadInput", warningstr);
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("smearing_sigma");
        item.annotation = "energy range for smearing";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "Energy range for smearing.";
        item.default_value = "0.015";
        item.unit = "Ry";
        item.availability = "";
        read_sync_double(input.smearing_sigma);
        this->add_item(item);
    }
    {
        // Energy range for smearing,
        //`smearing_sigma` = 1/2 *kB* `smearing_sigma_temp`.
        // NOTE: Use 'item' as the variable name for automatic documentation generation.
        Input_Item item("smearing_sigma_temp");
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "Energy range for smearing, smearing_sigma = 1/2 kB smearing_sigma_temp.";
        item.default_value = "2 * smearing_sigma / kB.";
        item.unit = "K";
        item.availability = "";
        item.read_value
            = [](const Input_Item& item, Parameter& para) { para.input.smearing_sigma = 3.166815e-6 * doublevalue; };
        // only to set smearing_sigma, so no need to write to output INPUT file
        // or bcast.
        this->add_item(item);
    }
    {
        Input_Item item("mixing_type");
        item.annotation = "plain; pulay; broyden";
        item.category = "Electronic structure";
        item.type = "String";
        item.description = R"(Charge mixing methods.
* plain: Just simple mixing.
* pulay: Standard Pulay method. P. Pulay Chemical Physics Letters, (1980)
* broyden: Simplified modified Broyden method. D.D. Johnson Physical Review B (1988)

In general, the convergence of the Broyden method is slightly faster than that of the Pulay method.)";
        item.default_value = "broyden";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.mixing_mode);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_beta");
        item.annotation = "mixing parameter: 0 means no new charge";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = R"(In general, the formula of charge mixing can be written as rho_new = rho_old + mixing_beta * drho, where rho_new represents the new charge density after charge mixing, rho_old represents the charge density in previous step, drho is obtained through various mixing methods, and mixing_beta is set by this parameter. A lower value of 'mixing_beta' results in less influence of drho on rho_new, making the self-consistent field (SCF) calculation more stable. However, it may require more steps to achieve convergence. We recommend the following options:
* 0.8: nspin=1
* 0.4: nspin=2 and nspin=4
* 0: keep charge density unchanged, usually used for restarting with init_chg=file or testing.
* 0.1 or less: if convergence of SCF calculation is difficult to reach, please try 0 < mixing_beta < 0.1.

Note: For low-dimensional large systems, the setup of mixing_beta=0.1, mixing_ndim=20, and mixing_gg0=1.0 usually works well.)";
        item.default_value = "0.8 for nspin=1, 0.4 for nspin=2 and nspin=4.";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("mixing_beta_mag");
        item.annotation = "mixing parameter for magnetic density";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "Mixing parameter of magnetic density.";
        item.default_value = "4*mixing_beta, but the maximum value is 1.6.";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("mixing_ndim");
        item.annotation = "mixing dimension in pulay or broyden";
        item.category = "Electronic structure";
        item.type = "Integer";
        item.description = R"(It indicates the mixing dimensions in Pulay or Broyden. Pulay and Broyden method use the density from previous mixing_ndim steps and do a charge mixing based on this density.

For systems that are difficult to converge, one could try increasing the value of 'mixing_ndim' to enhance the stability of the self-consistent field (SCF) calculation.)";
        item.default_value = "8";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.mixing_ndim);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_restart");
        item.annotation = "threshold to restart mixing during SCF";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "If the density difference between input and output drho is smaller than mixing_restart, SCF will restart at next step which means SCF will restart by using output charge density from perivos iteration as input charge density directly, and start a new mixing. Notice that mixing_restart will only take effect once in one SCF.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("mixing_dmr");
        item.annotation = "whether to mix real-space density matrix";
        item.category = "Electronic structure";
        item.type = "Boolean";
        item.description = "At n-th iteration which is calculated by drho<mixing_restart, SCF will start a mixing for real-space density matrix by using the same coefficiences as the mixing of charge density.";
        item.default_value = "false";
        item.unit = "";
        item.availability = "Only for mixing_restart>=0.0";
        read_sync_bool(input.mixing_dmr);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0");
        item.annotation = "mixing parameter in kerker";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = R"(Whether to perfom Kerker scaling for charge density.
* >0: The high frequency wave vectors will be suppressed by multiplying a scaling factor. Setting mixing_gg0 = 1.0 is normally a good starting point. Kerker preconditioner will be automatically turned off if mixing_beta <= 0.1.
* 0: No Kerker scaling is performed.

For systems that are difficult to converge, particularly metallic systems, enabling Kerker scaling may aid in achieving convergence.)";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.mixing_gg0);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0_mag");
        item.annotation = "mixing parameter in kerker";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "Whether to perfom Kerker preconditioner of magnetic density. Note: we do not recommand to open Kerker preconditioner of magnetic density unless the system is too hard to converge.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.mixing_gg0_mag);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0_min");
        item.annotation = "the minimum kerker coefficient";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "The minimum kerker coefficient.";
        item.default_value = "0.1";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.mixing_gg0_min);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_angle");
        item.annotation = "angle mixing parameter for non-colinear calculations";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = R"(Normal broyden mixing can give the converged result for a given magnetic configuration. If one is not interested in the energies of a given magnetic configuration but wants to determine the ground state by relaxing the magnetic moments' directions, one cannot rely on the standard Broyden mixing algorithm. To enhance the ability to find correct magnetic configuration for non-colinear calculations, ABACUS implements a promising mixing method proposed by J. Phys. Soc. Jpn. 82 (2013) 114706. Here, mixing_angle is the angle mixing parameter. In fact, only mixing_angle=1.0 is implemented currently.
* <=0: Normal broyden mixing
* >0: Angle mixing for the modulus with mixing_angle=1.0)";
        item.default_value = "-10.0";
        item.unit = "";
        item.availability = "Only relevant for non-colinear calculations nspin=4.";
        read_sync_double(input.mixing_angle);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_tau");
        item.annotation = "whether to mix tau in mGGA calculation";
        item.category = "Electronic structure";
        item.type = "Boolean";
        item.description = R"(Whether to mix the kinetic energy density.
* True: The kinetic energy density will also be mixed. It seems for general cases, SCF converges fine even without this mixing. However, if there is difficulty in converging SCF for meta-GGA, it might be helpful to turn this on.
* False: The kinetic energy density will not be mixed.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "Only relevant for meta-GGA calculations.";
        read_sync_bool(input.mixing_tau);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_dftu");
        item.annotation = "whether to mix locale in DFT+U calculation";
        item.category = "Electronic structure";
        item.type = "Boolean";
        item.description = R"(Whether to mix the occupation matrices.
* True: The occupation matrices will also be mixed by plain mixing. From experience this is not very helpful if the +U calculation does not converge.
* False: The occupation matrices will not be mixed.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "Only relevant for DFT+U calculations.";
        read_sync_bool(input.mixing_dftu);
        this->add_item(item);
    }
    {
        Input_Item item("gamma_only");
        item.annotation = "Only for localized orbitals set and gamma point. If "
                          "set to 1, a fast algorithm is used";
        item.category = "Electronic structure";
        item.type = "Boolean";
        item.description = R"(Whether to use gamma_only algorithm.
* 0: more than one k-point is used and the ABACUS is slower compared to the gamma only algorithm.
* 1: ABACUS uses gamma only, the algorithm is faster and you don't need to specify the k-points file.

Note: If gamma_only is set to 1, the KPT file will be overwritten. So make sure to turn off gamma_only for multi-k calculations.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "Only used in localized orbitals set";
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
                    ModuleBase::WARNING_QUIT("NOTICE", "nspin=4 (soc or noncollinear-spin) does not support gamma\n only calculation");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("scf_nmax");
        item.annotation = "number of electron iterations";
        item.category = "Electronic structure";
        item.type = "Integer";
        item.description = "This variable indicates the maximal iteration number for electronic iterations.";
        item.default_value = "100";
        item.unit = "";
        item.availability = "";
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
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "It's the density threshold for electronic iteration. It represents the charge density error between two sequential densities from electronic iterations. Usually for local orbitals, usually 1e-6 may be accurate enough.";
        item.default_value = "1.0e-9 (plane-wave basis), or 1.0e-7 (localized atomic orbital basis).";
        item.unit = "Ry if scf_thr_type=1, dimensionless if scf_thr_type=2";
        item.availability = "";
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
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "It's the energy threshold for electronic iteration. It represents the total energy error between two sequential densities from electronic iterations.";
        item.default_value = "-1.0. If the user does not set this parameter, it will not take effect.";
        item.unit = "eV";
        item.availability = "";
        read_sync_double(input.scf_ene_thr);
        this->add_item(item);
    }
    {
        Input_Item item("scf_thr_type");
        item.annotation = "type of the criterion of scf_thr, 1: reci drho for "
                          "pw, 2: real drho for lcao";
        item.category = "Electronic structure";
        item.type = "Integer";
        item.description = R"(Choose the calculation method of convergence criterion.
* 1: the criterion is defined in reciprocal space, which is used in SCF of PW basis with unit Ry.
* 2: the criterion is defined in real space, where is the number of electron, which is used in SCF of LCAO with unit dimensionless.)";
        item.default_value = "1 (plane-wave basis), or 2 (localized atomic orbital basis).";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("scf_os_stop");
        item.annotation = "whether to stop scf when oscillation is detected";
        item.category = "Electronic structure";
        item.type = "Boolean";
        item.description = R"(For systems that are difficult to converge, the SCF process may exhibit oscillations in charge density, preventing further progress toward the specified convergence criteria and resulting in continuous oscillation until the maximum number of steps is reached; this greatly wastes computational resources. To address this issue, this function allows ABACUS to terminate the SCF process early upon detecting oscillations, thus reducing subsequent meaningless calculations. The detection of oscillations is based on the slope of the logarithm of historical drho values. To this end, Least Squares Method is used to calculate the slope of the logarithmically taken drho for the previous scf_os_ndim iterations. If the calculated slope is larger than scf_os_thr, stop the SCF.

* 0: The SCF will continue to run regardless of whether there is oscillation or not.
* 1: If the calculated slope is larger than scf_os_thr, stop the SCF.)";
        item.default_value = "false";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.scf_os_stop);
        this->add_item(item);
    }
    {
        Input_Item item("scf_os_thr");
        item.annotation = "charge density threshold for oscillation";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = "The slope threshold to determine if the SCF is stuck in a charge density oscillation. If the calculated slope is larger than scf_os_thr, stop the SCF.";
        item.default_value = "-0.01";
        item.unit = "";
        item.availability = "";
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
        item.category = "Electronic structure";
        item.type = "Integer";
        item.description = "To determine the number of old iterations' drho used in slope calculations.";
        item.default_value = "mixing_ndim";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("sc_os_ndim");
        item.annotation = "number of old iterations used for oscillation detection, for Spin-Constrained DFT";
        item.category = "Electronic structure";
        item.type = "Integer";
        item.description = "To determine the number of old iterations to judge oscillation, it occured, more accurate lambda with DeltaSpin method would be calculated, only for PW base.";
        item.default_value = "5";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.sc_os_ndim);
        this->add_item(item);
    }
    {
        Input_Item item("lspinorb");
        item.annotation = "consider the spin-orbit interaction";
        item.category = "Electronic structure";
        item.type = "Boolean";
        item.description = R"(Whether to consider spin-orbit coupling (SOC) effect in the calculation.
* True: Consider spin-orbit coupling effect. When enabled:
  * nspin is automatically set to 4 (noncollinear spin representation)
  * Symmetry is automatically disabled (SOC breaks inversion symmetry)
  * Requires full-relativistic pseudopotentials with has_so=true in the UPF header
* False: Do not consider spin-orbit coupling effect.
* Common Error: "no soc upf used for lspinorb calculation" - ensure you are using full-relativistic pseudopotentials)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.lspinorb);
        this->add_item(item);
    }
    {
        Input_Item item("noncolin");
        item.annotation = "using non-collinear-spin";
        item.category = "Electronic structure";
        item.type = "Boolean";
        item.description = R"(Whether to allow non-collinear magnetic moments, where magnetization can point in arbitrary directions (x, y, z components) rather than being constrained to the z-axis.
* True: Allow non-collinear polarization. When enabled:
  * nspin is automatically set to 4
  * Wave function dimension is doubled (npol=2), and the number of occupied states is doubled
  * Charge density has 4 components (Pauli spin matrices)
  * Cannot be used with gamma_only=true
  * Can be combined with lspinorb=true for SOC effects with non-collinear magnetism
* False: Do not allow non-collinear polarization (magnetization constrained to z-axis).
* Relationship with lspinorb:
  * noncolin=0, lspinorb=1: SOC with z-axis magnetism only (for non-magnetic materials with SOC)
  * noncolin=1, lspinorb=0: Non-collinear magnetism without SOC
  * noncolin=1, lspinorb=1: Both non-collinear magnetism and SOC)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.noncolin);
        this->add_item(item);
    }
    {
        Input_Item item("soc_lambda");
        item.annotation = "The fraction of SOC based on scalar relativity (SR) of the pseudopotential";
        item.category = "Electronic structure";
        item.type = "Real";
        item.description = R"(Modulates the strength of spin-orbit coupling effect. Sometimes, for some real materials, both scalar-relativistic and full-relativistic pseudopotentials cannot describe the exact spin-orbit coupling. Artificial modulation may help in such cases.

soc_lambda, which has value range [0.0, 1.0], is used to modulate SOC effect:
* soc_lambda 0.0: Scalar-relativistic case (no SOC)
* soc_lambda 1.0: Full-relativistic case (full SOC)
* Intermediate values: Partial-relativistic SOC (interpolation between scalar and full)

Use case: When experimental or high-level theoretical results suggest that the SOC effect is weaker or stronger than what full-relativistic pseudopotentials predict, you can adjust this parameter to match the target behavior.)";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "Only works when lspinorb=true";
        read_sync_double(input.soc_lambda);
        this->add_item(item);
    }
    {
        Input_Item item("dfthalf_type");
        item.annotation = "DFT-1/2 type, 0:off; 1:shell DFT-1/2";
        item.category = "Electronic structure";
        item.type = "Integer";
        item.description = "DFT-1/2 type:\n* 0: DFT-1/2 is off.\n* 1: Shell DFT-1/2 method is used.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.dfthalf_type);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_thr");
        item.annotation = "threshold for eigenvalues is cg electron iterations";
        item.category = "Plane wave related variables";
        item.type = "Real";
        item.description = "Only used when you use ks_solver = cg/dav/dav_subspace/bpcg. It indicates the threshold for the first electronic iteration, from the second iteration the pw_diag_thr will be updated automatically. For nscf calculations with planewave basis set, pw_diag_thr should be <= 1e-3.";
        item.default_value = "0.01";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pw_diag_thr);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_s" && para.input.basis_type == "pw")
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
        Input_Item item("diago_smooth_ethr");
        item.annotation = "smooth ethr for iter methods";
        item.category = "Plane wave related variables";
        item.type = "Boolean";
        item.description = "If TRUE, the smooth threshold strategy, which applies a larger threshold (10e-5) for the empty states, will be implemented in the diagonalization methods. (This strategy should not affect total energy, forces, and other ground-state properties, but computational efficiency will be improved.) If FALSE, the smooth threshold strategy will not be applied.";
        item.default_value = "false";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.diago_smooth_ethr);
        this->add_item(item);
    }
    {
        Input_Item item("use_k_continuity");
        item.annotation = "whether to use k-point continuity for initializing wave functions";
        item.category = "Plane wave related variables";
        item.type = "Boolean";
        item.description = "If TRUE, the wavefunctions at k-point will be initialized from the converged wavefunctions at the nearest k-point, which can speed up the SCF convergence. Only works for PW basis.";
        item.default_value = "false";
        item.unit = "";
        item.availability = "Used only for plane wave basis set.";
        read_sync_bool(input.use_k_continuity);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.use_k_continuity && para.input.basis_type != "pw") {
                ModuleBase::WARNING_QUIT("ReadInput", "use_k_continuity only works for PW basis");
            }
            if (para.input.use_k_continuity && para.input.calculation == "nscf") {
                ModuleBase::WARNING_QUIT("ReadInput", "use_k_continuity cannot work for NSCF calculation");
            }
            if (para.input.use_k_continuity && para.input.nspin == 2) {
                ModuleBase::WARNING_QUIT("ReadInput", "use_k_continuity cannot work for spin-polarized calculation");
            }
            if (para.input.use_k_continuity && para.input.esolver_type == "sdft") {
                ModuleBase::WARNING_QUIT("ReadInput", "use_k_continuity cannot work for SDFT calculation");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_nmax");
        item.annotation = "max iteration number for cg";
        item.category = "Plane wave related variables";
        item.type = "Integer";
        item.description = "Only useful when you use ks_solver = cg/dav/dav_subspace/bpcg. It indicates the maximal iteration number for cg/david/dav_subspace/bpcg method.";
        item.default_value = "40";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pw_diag_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_ndim");
        item.annotation = "dimension of workspace for Davidson diagonalization";
        item.category = "Plane wave related variables";
        item.type = "Integer";
        item.description = "Only useful when you use ks_solver = dav or ks_solver = dav_subspace. It indicates dimension of workspace(number of wavefunction packets, at least 2 needed) for the Davidson method. A larger value may yield a smaller number of iterations in the algorithm but uses more memory and more CPU time in subspace diagonalization.";
        item.default_value = "4";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pw_diag_ndim);
        this->add_item(item);
    }
    {
        Input_Item item("diago_cg_prec");
        item.annotation = "diago_cg_prec";
        item.category = "Plane wave related variables";
        item.type = "Integer";
        item.description = "Preconditioner type for conjugate gradient diagonalization method.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.diago_cg_prec);
        this->add_item(item);
    }

    // LCAO
    {
        Input_Item item("nb2d");
        item.annotation = "matrix 2d division";
        item.category = "System variables";
        item.type = "Integer";
        item.description = "In LCAO calculations, the Hamiltonian and overlap matrices are distributed across 2D processor grid. This parameter controls the 2D block size for distribution.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
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
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Integer";
        item.description = "If not equals to 2, then the maximum l channels on LCAO is set to lmaxmax. If 2, then the number of l channels will be read from the LCAO data sets. Normally no input should be supplied for this variable so that it is kept as its default.";
        item.default_value = "2.";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.lmaxmax);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_ecut");
        item.annotation = "energy cutoff for LCAO";
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Real";
        item.description = "Energy cutoff (in Ry) for two-center integrals in LCAO. The two-center integration table are obtained via a k space integral whose upper limit is about sqrt(lcao_ecut).";
        item.default_value = "ecutwfc";
        item.unit = "";
        item.availability = "";
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
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Real";
        item.description = "the interval of k points for two-center integrals. The two-center integration table are obtained via a k space integral on a uniform grid with spacing lcao_dk.";
        item.default_value = "0.01";
        item.unit = "Bohr";
        item.availability = "";
        read_sync_double(input.lcao_dk);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_dr");
        item.annotation = "delta r for 1D integration in LCAO";
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Real";
        item.description = "r spacing of the integration table of two-center integrals.";
        item.default_value = "0.01";
        item.unit = "Bohr";
        item.availability = "";
        read_sync_double(input.lcao_dr);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_rmax");
        item.annotation = "max R for 1D two-center integration table";
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Real";
        item.description = "Maximum distance for the two-center integration table.";
        item.default_value = "30";
        item.unit = "Bohr";
        item.availability = "";
        read_sync_double(input.lcao_rmax);
        this->add_item(item);
    }
    {
        Input_Item item("search_radius");
        item.annotation = "input search radius (Bohr)";
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Real";
        item.description = "Searching radius in finding the neighbouring atoms. By default the radius will be automatically determined by the cutoffs of orbitals and nonlocal beta projectors.";
        item.default_value = "-1";
        item.unit = "Bohr";
        item.availability = "";
        read_sync_double(input.search_radius);
        this->add_item(item);
    }
    {
        Input_Item item("bx");
        item.annotation = "division of an element grid in FFT grid along x";
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Integer";
        item.description = "In the matrix operation of grid integral, bx/by/bz grids (in x, y, z directions) are treated as a whole as a matrix element. A different value will affect the calculation speed. The default is 0, which means abacus will automatically calculate these values.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
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
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Integer";
        item.description = "In the matrix operation of grid integral, bx/by/bz grids (in x, y, z directions) are treated as a whole as a matrix element. A different value will affect the calculation speed. The default is 0, which means abacus will automatically calculate these values.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
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
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Integer";
        item.description = "In the matrix operation of grid integral, bx/by/bz grids (in x, y, z directions) are treated as a whole as a matrix element. A different value will affect the calculation speed. The default is 0, which means abacus will automatically calculate these values.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
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
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Integer";
        item.description = "Number of threads used in one elpa calculation.\n\nIf the number is below 0 or 0 or beyond the max number of threads, all elpa calculation will be using all mpi threads";
        item.default_value = "-1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.elpa_num_thread);
        this->add_item(item);
    }
    {
        Input_Item item("num_stream");
        item.annotation = "the nstream in compute the LCAO with CUDA";
        item.category = "Numerical atomic orbitals related variables";
        item.type = "Integer";
        item.description = "The number of CUDA streams used in LCAO calculations with GPU acceleration.";
        item.default_value = "4";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.nstream);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_ecut");
        item.annotation = "energy cutoff for spherical bessel functions(Ry)";
        item.category = "NAOs";
        item.type = "String";
        item.description = "\"Energy cutoff\" (in Ry) of spherical Bessel functions. The number of spherical Bessel functions that constitute the radial parts of NAOs is determined by sqrt(bessel_nao_ecut)*bessel_nao_rcut/.";
        item.default_value = "ecutwfc";
        item.unit = "";
        item.availability = "";
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
        item.category = "NAOs";
        item.type = "Real";
        item.description = "Tolerance when searching for the zeros of spherical Bessel functions.";
        item.default_value = "1.0e-12";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.bessel_nao_tolerence);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_rcut");
        item.annotation = "radial cutoff for spherical bessel functions(a.u.)";
        item.category = "NAOs";
        item.type = "Real";
        item.description = "Cutoff radius (in Bohr) and the common node of spherical Bessel functions used to construct the NAOs.";
        item.default_value = "6.0";
        item.unit = "";
        item.availability = "";
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
        item.category = "NAOs";
        item.type = "Boolean";
        item.description = "If True, NAOs will be smoothed near the cutoff radius. See bessel_nao_rcut and bessel_nao_sigma for parameters.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.bessel_nao_smooth);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_nao_sigma");
        item.annotation = "spherical bessel smearing_sigma";
        item.category = "NAOs";
        item.type = "Real";
        item.description = "Smoothing range (in Bohr). See also bessel_nao_smooth.";
        item.default_value = "0.1";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.bessel_nao_sigma);
        this->add_item(item);
    }
}
} // namespace ModuleIO
