#include "module_base/constants.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_exx()
{
    // EXX
    {
        Input_Item item("exx_hybrid_alpha");
        item.annotation = "fraction of Fock exchange in hybrid functionals";
        read_sync_string(input.exx_hybrid_alpha);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.exx_hybrid_alpha == "default")
            {
                std::string& dft_functional = para.input.dft_functional;
                std::string dft_functional_lower = dft_functional;
                std::transform(dft_functional.begin(), dft_functional.end(), dft_functional_lower.begin(), tolower);
                if (dft_functional_lower == "hf")
                {
                    para.input.exx_hybrid_alpha = "1";
                }
                else if (dft_functional_lower == "pbe0" || dft_functional_lower == "hse"
                         || dft_functional_lower == "scan0")
                {
                    para.input.exx_hybrid_alpha = "0.25";
                }
                else
                { // no exx in scf, but will change to non-zero in
                    // postprocess like rpa
                    para.input.exx_hybrid_alpha = "0";
                }
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            const double exx_hybrid_alpha_value = std::stod(para.input.exx_hybrid_alpha);
            if (exx_hybrid_alpha_value < 0 || exx_hybrid_alpha_value > 1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "must 0 <= exx_hybrid_alpha <= 1");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_hse_omega");
        item.annotation = "range-separation parameter in HSE functional";
        read_sync_double(input.exx_hse_omega);
        this->add_item(item);
    }
    {
        Input_Item item("exx_separate_loop");
        item.annotation = "if 1, a two-step method is employed, else it will "
                          "start with a GGA-Loop, and then Hybrid-Loop";
        read_sync_bool(input.exx_separate_loop);
        this->add_item(item);
    }
    {
        Input_Item item("exx_hybrid_step");
        item.annotation = "the maximal electronic iteration number in the "
                          "evaluation of Fock exchange";
        read_sync_int(input.exx_hybrid_step);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.exx_hybrid_step <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exx_hybrid_step must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_mixing_beta");
        item.annotation = "mixing_beta for outer-loop when exx_separate_loop=1";
        read_sync_double(input.exx_mixing_beta);
        this->add_item(item);
    }
    {
        Input_Item item("exx_lambda");
        item.annotation = "used to compensate for divergence points at G=0 in "
                          "the evaluation of Fock exchange using "
                          "lcao_in_pw method";
        read_sync_double(input.exx_lambda);
        this->add_item(item);
    }
    {
        Input_Item item("exx_real_number");
        item.annotation = "exx calculated in real or complex";
        read_sync_string(input.exx_real_number);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.exx_real_number == "default")
            {  // to run through here, the default value of para.input.exx_real_number should be "default"
                if (para.input.gamma_only)
                {
                    para.input.exx_real_number = "1";
                }
                else
                {
                    para.input.exx_real_number = "0";
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_pca_threshold");
        item.annotation = "threshold to screen on-site ABFs in exx";
        read_sync_double(input.exx_pca_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_c_threshold");
        item.annotation = "threshold to screen C matrix in exx";
        read_sync_double(input.exx_c_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_v_threshold");
        item.annotation = "threshold to screen C matrix in exx";
        read_sync_double(input.exx_v_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_dm_threshold");
        item.annotation = "threshold to screen density matrix in exx";
        read_sync_double(input.exx_dm_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_schwarz_threshold");
        item.annotation = "threshold to screen exx using Cauchy-Schwartz inequality";
        read_sync_double(input.exx_schwarz_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_cauchy_threshold");
        item.annotation = "threshold to screen exx using Cauchy-Schwartz inequality";
        read_sync_double(input.exx_cauchy_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_c_grad_threshold");
        item.annotation = "threshold to screen nabla C matrix in exx";
        read_sync_double(input.exx_c_grad_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_v_grad_threshold");
        item.annotation = "threshold to screen nabla V matrix in exx";
        read_sync_double(input.exx_v_grad_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_c_grad_r_threshold");
        item.annotation = "threshold to screen nabla C matrix in exx";
        read_sync_double(input.exx_c_grad_r_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_v_grad_r_threshold");
        item.annotation = "threshold to screen nabla V matrix in exx";
        read_sync_double(input.exx_v_grad_r_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_cauchy_force_threshold");
        item.annotation = "threshold to screen exx force using Cauchy-Schwartz inequality";
        read_sync_double(input.exx_cauchy_force_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_cauchy_stress_threshold");
        item.annotation = "threshold to screen exx stress using Cauchy-Schwartz inequality";
        read_sync_double(input.exx_cauchy_stress_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_ccp_rmesh_times");
        item.annotation = "how many times larger the radial mesh required for "
                          "calculating Columb potential is to that "
                          "of atomic orbitals";
        read_sync_string(input.exx_ccp_rmesh_times);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.exx_ccp_rmesh_times == "default")
            {   // to run through here, the default value of para.input.exx_ccp_rmesh_times should be "default"
                std::string& dft_functional = para.input.dft_functional;
                std::string dft_functional_lower = dft_functional;
                std::transform(dft_functional.begin(), dft_functional.end(), dft_functional_lower.begin(), tolower);
                if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0" || dft_functional_lower == "scan0")
                {
                    para.input.exx_ccp_rmesh_times = "5";
                }
                else if (dft_functional_lower == "hse")
                {
                    para.input.exx_ccp_rmesh_times = "1.5";
                }
                else
                { // no exx in scf
                    para.input.exx_ccp_rmesh_times = "1";
                }
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (std::stod(para.input.exx_ccp_rmesh_times) < 1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exx_ccp_rmesh_times must >= 1");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_distribute_type");
        item.annotation = "exx_distribute_type";
        read_sync_string(input.exx_distribute_type);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.exx_distribute_type != "htime" && para.input.exx_distribute_type != "kmeans2"
                && para.input.exx_distribute_type != "kmeans1" && para.input.exx_distribute_type != "order")
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                                         "exx_distribute_type must be htime or "
                                         "kmeans2 or kmeans1 or order");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_opt_orb_lmax");
        item.annotation = "the maximum l of the spherical Bessel functions for opt ABFs";
        read_sync_int(input.exx_opt_orb_lmax);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.exx_opt_orb_lmax < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exx_opt_orb_lmax must >= 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_opt_orb_ecut");
        item.annotation = "the cut-off of plane wave expansion for opt ABFs";
        read_sync_double(input.exx_opt_orb_ecut);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.exx_opt_orb_ecut < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exx_opt_orb_ecut must >= 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_opt_orb_tolerence");
        item.annotation = "the threshold when solving for the zeros of "
                          "spherical Bessel functions for opt ABFs";
        read_sync_double(input.exx_opt_orb_tolerence);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.exx_opt_orb_tolerence < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exx_opt_orb_tolerence must >= 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_symmetry_realspace");
        item.annotation = "whether to reduce real-space sector in Hexx calculation";
        read_sync_bool(input.exx_symmetry_realspace);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.symmetry != "1") { para.input.exx_symmetry_realspace = false; }
            };
        this->add_item(item);
    }
    {
        Input_Item item("rpa_ccp_rmesh_times");
        item.annotation = "how many times larger the radial mesh required for "
                          "calculating Columb potential is to that "
                          "of atomic orbitals";
        read_sync_double(input.rpa_ccp_rmesh_times);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rpa_ccp_rmesh_times < 1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "rpa_ccp_rmesh_times must >= 1");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_ri_cv");
        item.annotation = "Whether to output the coefficient tensor C and ABFs-representation Coulomb matrix V";
        read_sync_bool(input.out_ri_cv);
        this->add_item(item);
    }
}
void ReadInput::item_dftu()
{
    // dft+u
    {
        Input_Item item("dft_plus_u");
        item.annotation = "DFT+U correction method";
        read_sync_int(input.dft_plus_u);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            bool all_minus1 = true;
            for (auto& val: para.input.orbital_corr)
            {
                if (val != -1)
                {
                    all_minus1 = false;
                    break;
                }
            }
            if (all_minus1)
            {
                if (para.input.dft_plus_u != 0)
                {
                    para.input.dft_plus_u = 0;
                    ModuleBase::WARNING("ReadInput", "No atoms are correlated, DFT+U is closed!!!");
                }
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            const Input_para& input = para.input;
            if (input.dft_plus_u != 0)
            {
                if (input.basis_type != "lcao")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "WRONG ARGUMENTS OF basis_type, only lcao is support");
                }
                if (input.ks_solver != "genelpa" && input.ks_solver != "scalapack_gvx" && input.ks_solver != "default")
                {
                    std::cout << " You'are using " << input.ks_solver << std::endl;
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "WRONG ARGUMENTS OF ks_solver in DFT+U routine, only "
                                             "genelpa and scalapack_gvx are supported ");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("dft_plus_dmft");
        item.annotation = "true:DFT+DMFT; false: standard DFT calcullation(default)";
        read_sync_bool(input.dft_plus_dmft);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.basis_type != "lcao" && para.input.dft_plus_dmft)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "DFT+DMFT is only supported for lcao calculation.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("yukawa_lambda");
        item.annotation = "default:0.0";
        read_sync_double(input.yukawa_lambda);
        this->add_item(item);
    }
    {
        Input_Item item("yukawa_potential");
        item.annotation = "default: false";
        read_sync_bool(input.yukawa_potential);
        this->add_item(item);
    }
    {
        Input_Item item("uramping");
        item.annotation = "increasing U values during SCF";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.uramping_eV = doublevalue;
            para.sys.uramping = para.input.uramping_eV / ModuleBase::Ry_to_eV;
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            bool all_minus1 = true;
            for (auto& val: para.input.orbital_corr)
            {
                if (val != -1)
                {
                    all_minus1 = false;
                    break;
                }
            }
            if (all_minus1)
            {
                if (para.sys.uramping != 0.0)
                {
                    para.sys.uramping = 0.0;
                    ModuleBase::WARNING("ReadInput", "No atoms are correlated, U-ramping is closed!!!");
                }
            }
        };
        sync_double(input.uramping_eV);
        this->add_item(item);
    }
    {
        Input_Item item("omc");
        item.annotation = "the mode of occupation matrix control";
        read_sync_int(input.omc);
        this->add_item(item);
    }
    {
        Input_Item item("onsite_radius");
        item.annotation = "radius of the sphere for onsite projection (Bohr)";
        read_sync_double(input.onsite_radius);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if ((para.input.dft_plus_u == 1 || para.input.sc_mag_switch) && para.input.onsite_radius == 0.0)
            {
                // autoset onsite_radius to 3.0 as default, this default value comes from the systematic atomic magnetism test
                para.input.onsite_radius = 3.0;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("hubbard_u");
        item.annotation = "Hubbard Coulomb interaction parameter U(ev)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.hubbard_u_eV.push_back(std::stod(item.str_values[i]));
                para.sys.hubbard_u.push_back(para.input.hubbard_u_eV[i] / ModuleBase::Ry_to_eV);
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (!item.is_read())
            {
                return;
            }
            if (para.sys.hubbard_u.size() != para.input.ntype)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                                         "hubbard_u should have the same "
                                         "number of elements as ntype");
            }
            for (auto& value: para.sys.hubbard_u)
            {
                if (value < -1.0e-3)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "WRONG ARGUMENTS OF hubbard_u");
                }
            }
        };
        sync_doublevec(input.hubbard_u_eV, para.input.ntype, 0.0);
        add_doublevec_bcast(sys.hubbard_u, para.input.ntype, 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("orbital_corr");
        item.annotation = "which correlated orbitals need corrected ; d:2 "
                          ",f:3, do not need correction:-1";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.orbital_corr.push_back(std::stoi(item.str_values[i]));
            }
        };

        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (!item.is_read())
            {
                return;
            }
            if (para.input.orbital_corr.size() != para.input.ntype)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                                         "orbital_corr should have the same "
                                         "number of elements as ntype");
            }
            for (auto& val: para.input.orbital_corr)
            {
                if (val < -1 || val > 3)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "WRONG ARGUMENTS OF orbital_corr");
                }
            }
        };
        sync_intvec(input.orbital_corr, para.input.ntype, -1);
        this->add_item(item);
    }
}
} // namespace ModuleIO