#include "source_base/global_function.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

#include <algorithm>
#include <cstring>
#include <iostream>

namespace ModuleIO
{
void ReadInput::item_others()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    // non-collinear spin-constrained
    {
        Input_Item item("sc_mag_switch");
        item.annotation = "switch to control spin-constrained DFT";
        item.category = "Spin-Constrained DFT";
        item.type = "Boolean";
        item.description = "Switch to control spin-constrained DFT calculation";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.sc_mag_switch);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_mag_switch)
            {
//                ModuleBase::WARNING_QUIT("ReadInput",
//                                         "This feature is not stable yet and might lead to "
//                                         "erroneous results.\n"
//                                         " Please wait for the official release version.");
                // if (para.input.nspin != 4 && para.input.nspin != 2)
                // {
                //     ModuleBase::WARNING_QUIT("ReadInput", "nspin must be 2 or
                //     4 when sc_mag_switch > 0");
                // }
                // if (para.input.calculation != "scf")
                // {
                //     ModuleBase::WARNING_QUIT("ReadInput", "calculation must
                //     be scf when sc_mag_switch > 0");
                // }
                // if (para.input.nupdown > 0.0)
                // {
                //     ModuleBase::WARNING_QUIT("ReadInput", "nupdown should not
                //     be set when sc_mag_switch > 0");
                // }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("decay_grad_switch");
        item.annotation = "switch to control gradient break condition";
        item.category = "Spin-Constrained DFT";
        item.type = "Boolean";
        item.description = "Switch to control gradient break condition in spin-constrained DFT";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.decay_grad_switch);
        this->add_item(item);
    }
    {
        Input_Item item("sc_thr");
        item.annotation = "Convergence criterion of spin-constrained iteration (RMS) in uB";
        item.category = "Spin-Constrained DFT";
        item.type = "Real";
        item.description = "Convergence criterion of spin-constrained iteration (RMS) in uB";
        item.default_value = "1.0e-6";
        item.unit = "uB";
        item.availability = "sc_mag_switch is true";
        read_sync_double(input.sc_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_thr < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sc_thr must >= 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nsc");
        item.annotation = "Maximal number of spin-constrained iteration";
        item.category = "Spin-Constrained DFT";
        item.type = "Integer";
        item.description = "Maximal number of spin-constrained iteration";
        item.default_value = "100";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_int(input.nsc);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nsc <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nsc must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nsc_min");
        item.annotation = "Minimum number of spin-constrained iteration";
        item.category = "Spin-Constrained DFT";
        item.type = "Integer";
        item.description = "Minimum number of spin-constrained iteration";
        item.default_value = "2";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_int(input.nsc_min);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nsc_min <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nsc_min must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("sc_scf_nmin");
        item.annotation = "Minimum number of outer scf loop before "
                          "initializing lambda loop";
        item.category = "Spin-Constrained DFT";
        item.type = "Integer";
        item.description = "Minimum number of outer scf loop before initializing lambda loop";
        item.default_value = "2";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_int(input.sc_scf_nmin);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_scf_nmin < 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sc_scf_nmin must >= 2");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("alpha_trial");
        item.annotation = "Initial trial step size for lambda in eV/uB^2";
        item.category = "Spin-Constrained DFT";
        item.type = "Real";
        item.description = "Initial trial step size for lambda in eV/uB^2";
        item.default_value = "0.01";
        item.unit = "eV/uB^2";
        item.availability = "sc_mag_switch is true";
        read_sync_double(input.alpha_trial);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.alpha_trial <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "alpha_trial must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("sccut");
        item.annotation = "Maximal step size for lambda in eV/uB";
        item.category = "Spin-Constrained DFT";
        item.type = "Real";
        item.description = "Maximal step size for lambda in eV/uB";
        item.default_value = "3.0";
        item.unit = "eV/uB";
        item.availability = "sc_mag_switch is true";
        read_sync_double(input.sccut);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sccut <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sccut must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("sc_drop_thr");
        item.annotation = "Convergence criterion ratio of lambda iteration in Spin-constrained DFT";
        item.category = "Spin-Constrained DFT";
        item.type = "Real";
        item.description = "Convergence criterion ratio of lambda iteration in Spin-constrained DFT";
        item.default_value = "1.0e-2";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_double(input.sc_drop_thr);
        this->add_item(item);
    }
    {
        Input_Item item("sc_scf_thr");
        item.annotation = "Density error threshold for inner loop of spin-constrained SCF";
        item.category = "Spin-Constrained DFT";
        item.type = "Real";
        item.description = "Density error threshold for inner loop of spin-constrained SCF";
        item.default_value = "1.0e-4";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_double(input.sc_scf_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_scf_thr <= 0.0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sc_scf_thr must > 0.0");
            }
        };
        this->add_item(item);
    }

    // Quasiatomic Orbital analysis
    {
        Input_Item item("qo_switch");
        item.annotation = "switch to control quasiatomic orbital analysis";
        item.category = "Quasiatomic Orbital (QO) analysis";
        item.type = "Boolean";
        item.description = "Whether to let ABACUS output QO analysis required files";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.qo_switch);
        this->add_item(item);
    }
    {
        Input_Item item("qo_basis");
        item.annotation = "type of QO basis function: hydrogen: hydrogen-like "
                          "basis, pswfc: read basis from pseudopotential";
        item.category = "Quasiatomic Orbital (QO) analysis";
        item.type = "String";
        item.description = R"(Type of QO basis function:
* hydrogen: hydrogen-like basis
* pswfc: read basis from pseudopotential
* szv: single-zeta valence basis)";
        item.default_value = "szv";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.qo_basis);
        this->add_item(item);
    }
    {
        Input_Item item("qo_strategy");
        item.annotation = "strategy to generate generate radial orbitals";
        item.category = "Quasiatomic Orbital (QO) analysis";
        item.type = "Vector of String (1 or n values where n is the number of atomic types)";
        item.description = "Strategy to generate radial orbitals for QO analysis. For hydrogen: energy-valence, for pswfc and szv: all";
        item.default_value = "for hydrogen: energy-valence, for pswfc and szv: all";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.qo_strategy.push_back(item.str_values[i]);
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.qo_strategy.size() != para.input.ntype)
            {
                if (para.input.qo_strategy.size() == 1)
                {
                    para.input.qo_strategy.resize(para.input.ntype, para.input.qo_strategy[0]);
                }
                else
                {
                    std::string default_strategy;
                    if (para.input.qo_basis == "hydrogen")
                    {
                        default_strategy = "energy-valence";
                    }
                    else if ((para.input.qo_basis == "pswfc") || (para.input.qo_basis == "szv"))
                    {
                        default_strategy = "all";
                    }
                    else
                    {
                        ModuleBase::WARNING_QUIT("ReadInput",
                                                 "When setting default values for qo_strategy, "
                                                 "unexpected/unknown "
                                                 "qo_basis is found. Please check it.");
                    }
                    para.input.qo_strategy.resize(para.input.ntype, default_strategy);
                }
            }
        };
        sync_stringvec(input.qo_strategy, para.input.ntype, "all");
        this->add_item(item);
    }
    {
        Input_Item item("qo_screening_coeff");
        item.annotation = "rescale the shape of radial orbitals";
        item.category = "Quasiatomic Orbital (QO) analysis";
        item.type = "Vector of Real (n values where n is the number of atomic types; 1 value allowed for qo_basis=pswfc)";
        item.description = "The screening coefficient for each atom type to rescale the shape of radial orbitals";
        item.default_value = "0.1";
        item.unit = "Bohr^-1";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.qo_screening_coeff.push_back(std::stod(item.str_values[i]));
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (!item.is_read())
            {
                return;
            }
            if (para.input.qo_screening_coeff.size() != para.input.ntype)
            {
                if (para.input.qo_basis == "pswfc")
                {
                    double default_screening_coeff
                        = (para.input.qo_screening_coeff.size() == 1) ? para.input.qo_screening_coeff[0] : 0.1;
                    para.input.qo_screening_coeff.resize(para.input.ntype, default_screening_coeff);
                }
                else
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "qo_screening_coeff should have the same number of "
                                             "elements as ntype");
                }
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            for (auto screen_coeff: para.input.qo_screening_coeff)
            {
                if (screen_coeff < 0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "screening coefficient must >= 0 "
                                             "to tune the pswfc decay");
                }
                if (std::fabs(screen_coeff) < 1e-6)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "every low screening coefficient might yield very high "
                                             "computational cost");
                }
            }
        };
        sync_doublevec(input.qo_screening_coeff, para.input.ntype, 0.1);
        this->add_item(item);
    }
    {
        Input_Item item("qo_thr");
        item.annotation = "accuracy for evaluating cutoff radius of QO basis function";
        item.category = "Quasiatomic Orbital (QO) analysis";
        item.type = "Real";
        item.description = "The convergence threshold determining the cutoff of generated orbital. Lower threshold will yield orbital with larger cutoff radius.";
        item.default_value = "1.0e-6";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.qo_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.qo_thr > 1e-6)
            {
                ModuleBase::WARNING("ReadInput",
                                    "too high the convergence threshold might "
                                    "yield unacceptable result");
            }
        };
        this->add_item(item);
    }

    // PEXSI
    {
        Input_Item item("pexsi_npole");
        item.annotation = "Number of poles in expansion";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "The number of poles used in the pole expansion method, should be a even number.";
        item.default_value = "40";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_npole);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_inertia");
        item.annotation = "Whether inertia counting is used at the very "
                          "beginning of PEXSI process";
        item.category = "PEXSI";
        item.type = "Boolean";
        item.description = "Whether inertia counting is used at the very beginning.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pexsi_inertia);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nmax");
        item.annotation = "Maximum number of PEXSI iterations after each "
                          "inertia counting procedure";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "Maximum number of PEXSI iterations after each inertia counting procedure.";
        item.default_value = "80";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_comm");
        item.annotation = "Whether to construct PSelInv communication pattern";
        item.category = "PEXSI";
        item.type = "Boolean";
        item.description = "Whether to construct PSelInv communication pattern.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pexsi_comm);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_storage");
        item.annotation = "Storage space used by the Selected Inversion "
                          "algorithm for symmetric matrices";
        item.category = "PEXSI";
        item.type = "Boolean";
        item.description = "Whether to use symmetric storage space used by the Selected Inversion algorithm for symmetric matrices.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pexsi_storage);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_ordering");
        item.annotation = "Ordering strategy for factorization and selected inversion";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "Ordering strategy for factorization and selected inversion. 0: Parallel ordering using ParMETIS, 1: Sequential ordering using METIS, 2: Multiple minimum degree ordering";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_ordering);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_row_ordering");
        item.annotation = "Row permutation strategy for factorization and "
                          "selected inversion, 0: NoRowPerm, 1: LargeDiag";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "Row permutation strategy for factorization and selected inversion, 0: No row permutation, 1: Make the diagonal entry of the matrix larger than the off-diagonal entries.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_row_ordering);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nproc");
        item.annotation = "Number of processors for parmetis";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "Number of processors for PARMETIS. Only used if pexsi_ordering == 0.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_nproc);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_symm");
        item.annotation = "Matrix symmetry";
        item.category = "PEXSI";
        item.type = "Boolean";
        item.description = "Whether the matrix is symmetric.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pexsi_symm);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_trans");
        item.annotation = "Whether to transpose";
        item.category = "PEXSI";
        item.type = "Boolean";
        item.description = "Whether to factorize the transpose of the matrix.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pexsi_trans);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_method");
        item.annotation = "pole expansion method, 1: Cauchy Contour Integral, "
                          "2: Moussa optimized method";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "The pole expansion method to be used. 1 for Cauchy Contour Integral method, 2 for Moussa optimized method.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_method);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nproc_pole");
        item.annotation = "Number of processes used by each pole";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "The point parallelizaion of PEXSI. Recommend two points parallelization.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_nproc_pole);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_temp");
        item.annotation = "Temperature, in the same unit as H";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Temperature in Fermi-Dirac distribution, in Ry, should have the same effect as the smearing sigma when smearing method is set to Fermi-Dirac.";
        item.default_value = "0.015";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_temp);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_gap");
        item.annotation = "Spectral gap";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Spectral gap, this can be set to be 0 in most cases.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_gap);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_delta_e");
        item.annotation = "An upper bound for the spectral radius of S^{-1} H";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Upper bound for the spectral radius of S^{-1}H.";
        item.default_value = "20";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_delta_e);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_lower");
        item.annotation = "Initial guess of lower bound for mu";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Initial guess of lower bound for mu.";
        item.default_value = "-10";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu_lower);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_upper");
        item.annotation = "Initial guess of upper bound for mu";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Initial guess of upper bound for mu.";
        item.default_value = "10";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu_upper);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu");
        item.annotation = "Initial guess for mu (for the solver)";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Initial guess for mu (for the solver).";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_thr");
        item.annotation = "Stopping criterion in terms of the chemical "
                          "potential for the inertia counting procedure";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Stopping criterion in terms of the chemical potential for the inertia counting procedure.";
        item.default_value = "0.05";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu_thr);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_expand");
        item.annotation = "If the chemical potential is not in the initial "
                          "interval, the interval is expanded by "
                          "muInertiaExpansion";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "If the chemical potential is not in the initial interval, the interval is expanded by this value.";
        item.default_value = "0.3";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu_expand);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_guard");
        item.annotation = "Safe guard criterion in terms of the chemical potential to "
                          "reinvoke the inertia counting procedure";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Safe guard criterion in terms of the chemical potential to reinvoke the inertia counting procedure.";
        item.default_value = "0.2";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu_guard);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_elec_thr");
        item.annotation = "Stopping criterion of the PEXSI iteration in terms "
                          "of the number of electrons compared to "
                          "numElectronExact";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Stopping criterion of the PEXSI iteration in terms of the number of electrons compared to numElectronExact.";
        item.default_value = "0.001";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_elec_thr);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_zero_thr");
        item.annotation = "if the absolute value of matrix element is less "
                          "than ZERO_Limit, it will be considered as 0";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "if the absolute value of CCS matrix element is less than this value, it will be considered as zero.";
        item.default_value = "1e-10";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_zero_thr);
        this->add_item(item);
    }

    // Only for Test
    {
        Input_Item item("out_alllog");
        item.annotation = "output information for each processor, when parallel";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print information into individual logs from all ranks in an MPI run.\n* True: Information from each rank will be written into individual files named OUT.{calculation}_{suffix}/running_${calculation}.log.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_alllog);
        this->add_item(item);
    }
    {
        Input_Item item("nurse");
        item.annotation = "for coders";
        item.category = "Variables useful for debugging";
        item.type = "Integer";
        item.description = "Debugging flag for developers";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.nurse);
        this->add_item(item);
    }
    {
        Input_Item item("t_in_h");
        item.annotation = "calculate the kinetic energy or not";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to include kinetic term in obtaining the Hamiltonian matrix.
* 0: No.
* 1: Yes.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.t_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vl_in_h");
        item.annotation = "calculate the local potential or not";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to include local pseudopotential term in obtaining the Hamiltonian matrix.
* 0: No.
* 1: Yes.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.vl_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vnl_in_h");
        item.annotation = "calculate the nonlocal potential or not";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to include non-local pseudopotential term in obtaining the Hamiltonian matrix.
* 0: No.
* 1: Yes.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.vnl_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vh_in_h");
        item.annotation = "calculate the hartree potential or not";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to include Hartree potential term in obtaining the Hamiltonian matrix.
* 0: No.
* 1: Yes.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.vh_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vion_in_h");
        item.annotation = "calculate the local ionic potential or not";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to include local ionic potential term in obtaining the Hamiltonian matrix.
* 0: No.
* 1: Yes.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.vion_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("test_force");
        item.annotation = "test the force";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to output the detailed components in forces.
* 0: No.
* 1: Yes.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.test_force);
        this->add_item(item);
    }
    {
        Input_Item item("test_stress");
        item.annotation = "test the stress";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to output the detailed components in stress.
* 0: No.
* 1: Yes.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.test_stress);
        this->add_item(item);
    }
    {
        Input_Item item("test_skip_ewald");
        item.annotation = "whether to skip ewald";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to skip the calculation of the ewald energy.
* 0: No.
* 1: Yes.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.test_skip_ewald);
        this->add_item(item);
    }
    {
        Input_Item item("ri_hartree_benchmark");
        item.annotation = "whether to use the RI approximation for the Hartree term in LR-TDDFT for benchmark (with FHI-aims/ABACUS read-in style)";
        item.category = "Linear Response TDDFT";
        item.type = "String";
        item.description = "Whether to use the RI approximation for the Hartree term in LR-TDDFT for benchmark (with FHI-aims/ABACUS read-in style)";
        item.default_value = "none";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.ri_hartree_benchmark);
        this->add_item(item);
    }
    {
        Input_Item item("aims_nbasis");
        item.annotation = "the number of basis functions for each atom type used in FHI-aims (for benchmark)";
        item.category = "Linear Response TDDFT";
        item.type = "A number(ntype) of Integers";
        item.description = "Atomic basis set size for each atom type (with the same order as in STRU) in FHI-aims.";
        item.default_value = "{} (empty list, where ABACUS use its own basis set size)";
        item.unit = "";
        item.availability = "ri_hartree_benchmark = aims";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.aims_nbasis.push_back(std::stod(item.str_values[i]));
            }
            };
        sync_intvec(input.aims_nbasis, para.input.aims_nbasis.size(), 0);
        this->add_item(item);
    }

    // RDMFT, added by jghan, 2024-10-16
    {
        Input_Item item("rdmft");
        item.annotation = "whether to perform rdmft calculation, default is false";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Boolean";
        item.description = "Whether to perform rdmft calculation (reduced density matrix funcional theory)";
        item.default_value = "false";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.rdmft);
        this->add_item(item);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && para.input.nspin == 4)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "rdmft is not available for nspin = 4");
            }
        };
    }
    {
        Input_Item item("rdmft_power_alpha");
        item.annotation = "the alpha parameter of power-functional, g(occ_number) = occ_number^alpha"
                          " used in exx-type functionals such as muller and power";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "The alpha parameter of power-functional(or other exx-type/hybrid functionals) which used in RDMFT, g(occ_number) = occ_number^alpha";
        item.default_value = "0.656";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.rdmft_power_alpha);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if( para.input.dft_functional == "hf" || para.input.dft_functional == "pbe0" )
            {
                para.input.rdmft_power_alpha = 1.0;
            }
            else if( para.input.dft_functional == "muller" )
            {
                para.input.rdmft_power_alpha = 0.5;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if( (para.input.rdmft_power_alpha < 0) || (para.input.rdmft_power_alpha > 1) )
            {
                ModuleBase::WARNING_QUIT("ReadInput", "rdmft_power_alpha should be greater than 0.0 and less than 1.0");
            }
        };
        this->add_item(item);
    }

    // EXX PW by rhx0820, 2025-03-10
    {
        Input_Item item("exxace");
        item.annotation = "whether to perform ace calculation in exxpw";
        item.category = "Exact Exchange (PW)";
        item.type = "Boolean";
        item.description = R"(Whether to use the ACE method (https://doi.org/10.1021/acs.jctc.6b00092) to accelerate the calculation the Fock exchange matrix. Should be set to true most of the time.
* True: Use the ACE method to calculate the Fock exchange operator.
* False: Use the traditional method to calculate the Fock exchange operator.)";
        item.default_value = "True";
        item.unit = "";
        item.availability = "exx_separate_loop==True.";
        read_sync_bool(input.exxace);
        this->add_item(item);
    }
    {
        Input_Item item("exx_gamma_extrapolation");
        item.annotation = "whether to perform gamma extrapolation in exxpw";
        item.category = "Exact Exchange (PW)";
        item.type = "Boolean";
        item.description = "Whether to use the gamma point extrapolation method to calculate the Fock exchange operator. See https://doi.org/10.1103/PhysRevB.79.205114 for details. Should be set to true most of the time.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.exx_gamma_extrapolation);
        this->add_item(item);
    }
    {
        Input_Item item("ecutexx");
        item.annotation = "energy cutoff for exx calculation, Ry";
        item.category = "Exact Exchange (PW)";
        item.type = "Real";
        item.description = "The energy cutoff for EXX (Fock) exchange operator in plane wave basis calculations. Reducing ecutexx below ecutrho may significantly accelerate EXX computations. This speed improvement comes with a reduced numerical accuracy in the exchange energy calculation.";
        item.default_value = "same as ecutrho";
        item.unit = "Ry";
        item.availability = "";
        read_sync_double(input.ecutexx);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.ecutexx < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ecutexx must >= 0");
            }
        };
        this->add_item(item);
    }

    {
        Input_Item item("exx_thr_type");
        item.annotation = "threshold type for exx outer loop, energy or density";
        item.category = "Exact Exchange (PW)";
        item.type = "String";
        item.description = R"(The type of threshold used to judge whether the outer loop has converged in the separate loop EXX calculation.
* energy: use the change of exact exchange energy to judge convergence.
* density: if the change of charge density difference between two successive outer loop iterations is seen as converged according to scf_thr, then the outer loop is seen as converged.)";
        item.default_value = "density";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.exx_thr_type);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            std::string thr_type = para.input.exx_thr_type;
            std::transform(thr_type.begin(), thr_type.end(), thr_type.begin(), ::tolower);
            if (thr_type != "energy" && thr_type != "density")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exx_thr_type should be energy or density");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_ene_thr");
        item.annotation = "threshold for exx outer loop when exx_thr_type = energy";
        item.category = "Exact Exchange (PW)";
        item.type = "Real";
        item.description = "The threshold for the change of exact exchange energy to judge convergence of the outer loop in the separate loop EXX calculation.";
        item.default_value = "1e-5";
        item.unit = "Ry";
        item.availability = "exx_thr_type==energy";
        read_sync_double(input.exx_ene_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.exx_ene_thr <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exx_ene_thr must > 0");
            }
        };
        this->add_item(item);
    }
}
} // namespace ModuleIO
