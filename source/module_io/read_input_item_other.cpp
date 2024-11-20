#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

#include <algorithm>
#include <cstring>
#include <iostream>

namespace ModuleIO
{
void ReadInput::item_others()
{
    // non-collinear spin-constrained
    {
        Input_Item item("sc_mag_switch");
        item.annotation = "switch to control spin-constrained DFT";
        read_sync_bool(input.sc_mag_switch);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_mag_switch)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                                         "This feature is not stable yet and might lead to "
                                         "erroneous results.\n"
                                         " Please wait for the official release version.");
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
        read_sync_bool(input.decay_grad_switch);
        this->add_item(item);
    }
    {
        Input_Item item("sc_thr");
        item.annotation = "Convergence criterion of spin-constrained iteration (RMS) in uB";
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
        read_sync_double(input.sc_drop_thr);
        this->add_item(item);
    }
    {
        Input_Item item("sc_scf_thr");
        item.annotation = "Density error threshold for inner loop of spin-constrained SCF";
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
        read_sync_bool(input.qo_switch);
        this->add_item(item);
    }
    {
        Input_Item item("qo_basis");
        item.annotation = "type of QO basis function: hydrogen: hydrogen-like "
                          "basis, pswfc: read basis from pseudopotential";
        read_sync_string(input.qo_basis);
        this->add_item(item);
    }
    {
        Input_Item item("qo_thr");
        item.annotation = "accuracy for evaluating cutoff radius of QO basis function";
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
    {
        Input_Item item("qo_strategy");
        item.annotation = "strategy to generate generate radial orbitals";
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

    // PEXSI
    {
        Input_Item item("pexsi_npole");
        item.annotation = "Number of poles in expansion";
        read_sync_int(input.pexsi_npole);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_inertia");
        item.annotation = "Whether inertia counting is used at the very "
                          "beginning of PEXSI process";
        read_sync_bool(input.pexsi_inertia);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nmax");
        item.annotation = "Maximum number of PEXSI iterations after each "
                          "inertia counting procedure";
        read_sync_int(input.pexsi_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_comm");
        item.annotation = "Whether to construct PSelInv communication pattern";
        read_sync_bool(input.pexsi_comm);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_storage");
        item.annotation = "Storage space used by the Selected Inversion "
                          "algorithm for symmetric matrices";
        read_sync_bool(input.pexsi_storage);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_ordering");
        item.annotation = "Ordering strategy for factorization and selected inversion";
        read_sync_int(input.pexsi_ordering);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_row_ordering");
        item.annotation = "Row permutation strategy for factorization and "
                          "selected inversion, 0: NoRowPerm, 1: LargeDiag";
        read_sync_int(input.pexsi_row_ordering);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nproc");
        item.annotation = "Number of processors for parmetis";
        read_sync_int(input.pexsi_nproc);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_symm");
        item.annotation = "Matrix symmetry";
        read_sync_bool(input.pexsi_symm);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_trans");
        item.annotation = "Whether to transpose";
        read_sync_bool(input.pexsi_trans);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_method");
        item.annotation = "pole expansion method, 1: Cauchy Contour Integral, "
                          "2: Moussa optimized method";
        read_sync_int(input.pexsi_method);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nproc_pole");
        item.annotation = "Number of processes used by each pole";
        read_sync_int(input.pexsi_nproc_pole);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_temp");
        item.annotation = "Temperature, in the same unit as H";
        read_sync_double(input.pexsi_temp);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_gap");
        item.annotation = "Spectral gap";
        read_sync_double(input.pexsi_gap);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_delta_e");
        item.annotation = "An upper bound for the spectral radius of S^{-1} H";
        read_sync_double(input.pexsi_delta_e);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_lower");
        item.annotation = "Initial guess of lower bound for mu";
        read_sync_double(input.pexsi_mu_lower);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_upper");
        item.annotation = "Initial guess of upper bound for mu";
        read_sync_double(input.pexsi_mu_upper);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu");
        item.annotation = "Initial guess for mu (for the solver)";
        read_sync_double(input.pexsi_mu);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_thr");
        item.annotation = "Stopping criterion in terms of the chemical "
                          "potential for the inertia counting procedure";
        read_sync_double(input.pexsi_mu_thr);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_expand");
        item.annotation = "If the chemical potential is not in the initial "
                          "interval, the interval is expanded by "
                          "muInertiaExpansion";
        read_sync_double(input.pexsi_mu_expand);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_guard");
        item.annotation = "Safe guard criterion in terms of the chemical potential to "
                          "reinvoke the inertia counting procedure";
        read_sync_double(input.pexsi_mu_guard);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_elec_thr");
        item.annotation = "Stopping criterion of the PEXSI iteration in terms "
                          "of the number of electrons compared to "
                          "numElectronExact";
        read_sync_double(input.pexsi_elec_thr);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_zero_thr");
        item.annotation = "if the absolute value of matrix element is less "
                          "than ZERO_Limit, it will be considered as 0";
        read_sync_double(input.pexsi_zero_thr);
        this->add_item(item);
    }

    // Only for Test
    {
        Input_Item item("out_alllog");
        item.annotation = "output information for each processor, when parallel";
        read_sync_bool(input.out_alllog);
        this->add_item(item);
    }
    {
        Input_Item item("nurse");
        item.annotation = "for coders";
        read_sync_int(input.nurse);
        this->add_item(item);
    }
    {
        Input_Item item("t_in_h");
        item.annotation = "calculate the kinetic energy or not";
        read_sync_bool(input.t_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vl_in_h");
        item.annotation = "calculate the local potential or not";
        read_sync_bool(input.vl_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vnl_in_h");
        item.annotation = "calculate the nonlocal potential or not";
        read_sync_bool(input.vnl_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vh_in_h");
        item.annotation = "calculate the hartree potential or not";
        read_sync_bool(input.vh_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vion_in_h");
        item.annotation = "calculate the local ionic potential or not";
        read_sync_bool(input.vion_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("test_force");
        item.annotation = "test the force";
        read_sync_bool(input.test_force);
        this->add_item(item);
    }
    {
        Input_Item item("test_stress");
        item.annotation = "test the stress";
        read_sync_bool(input.test_stress);
        this->add_item(item);
    }
    {
        Input_Item item("test_skip_ewald");
        item.annotation = "whether to skip ewald";
        read_sync_bool(input.test_skip_ewald);
        this->add_item(item);
    }
    {
        Input_Item item("ri_hartree_benchmark");
        item.annotation = "whether to use the RI approximation for the Hartree term in LR-TDDFT for benchmark (with FHI-aims/ABACUS read-in style)";
        read_sync_string(input.ri_hartree_benchmark);
        this->add_item(item);
    }
    {
        Input_Item item("aims_nbasis");
        item.annotation = "the number of basis functions for each atom type used in FHI-aims (for benchmark)";
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
}
} // namespace ModuleIO