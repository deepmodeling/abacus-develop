#include <cstring>

#include <algorithm>
#include <iostream>

#include "module_base/constants.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_others()
{
    // 10. Electric field and dipole correction
    {
        Input_Item item("efield_flag");
        item.annotation = "add electric field";
        read_sync_bool(input.efield_flag);
        this->add_item(item);
    }
    {
        Input_Item item("dip_cor_flag");
        item.annotation = "dipole correction";
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.dip_cor_flag && !para.input.efield_flag)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dipole correction is not active if efield_flag=false !");
            }
        };
        read_sync_bool(input.dip_cor_flag);
        this->add_item(item);
    }
    {
        Input_Item item("efield_dir");
        item.annotation = "the direction of the electric field or dipole correction";
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.gate_flag && para.input.efield_flag && !para.input.dip_cor_flag)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                                         "gate field cannot be used with "
                                         "efield if dip_cor_flag=false !");
            }
        };
        read_sync_int(input.efield_dir);
        this->add_item(item);
    }
    {
        Input_Item item("efield_pos_max");
        item.annotation = "position of the maximum of the saw-like potential "
                          "along crystal axis efield_dir";
        read_sync_double(input.efield_pos_max);
        this->add_item(item);
    }
    {
        Input_Item item("efield_pos_dec");
        item.annotation = "zone in the unit cell where the saw-like potential decreases";
        read_sync_double(input.efield_pos_dec);
        this->add_item(item);
    }
    {
        Input_Item item("efield_amp");
        item.annotation = "amplitude of the electric field";
        read_sync_double(input.efield_amp);
        this->add_item(item);
    }

    // 11. Gate field
    {
        Input_Item item("gate_flag");
        item.annotation = "compensating charge or not";
        read_sync_bool(input.gate_flag);
        this->add_item(item);
    }
    {
        Input_Item item("zgate");
        item.annotation = "position of charged plate";
        read_sync_double(input.zgate);
        this->add_item(item);
    }
    {
        Input_Item item("relax");
        item.annotation = "allow relaxation along the specific direction";
        read_sync_bool(input.relax);
        this->add_item(item);
    }
    {
        Input_Item item("block");
        item.annotation = "add a block potential or not";
        read_sync_bool(input.block);
        this->add_item(item);
    }
    {
        Input_Item item("block_down");
        item.annotation = "low bound of the block";
        read_sync_double(input.block_down);
        this->add_item(item);
    }
    {
        Input_Item item("block_up");
        item.annotation = "high bound of the block";
        read_sync_double(input.block_up);
        this->add_item(item);
    }
    {
        Input_Item item("block_height");
        item.annotation = "height of the block";
        read_sync_double(input.block_height);
        this->add_item(item);
    }

    // 12. Test
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
        Input_Item item("colour");
        item.annotation = "for coders, make their live colourful";
        read_sync_bool(input.colour);
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

    // 13. vdW Correction
    {
        Input_Item item("vdw_method");
        item.annotation = "the method of calculating vdw (none ; d2 ; d3_0 ; d3_bj";
        read_sync_string(input.vdw_method);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_s6");
        item.annotation = "scale parameter of d2/d3_0/d3_bj";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.vdw_s6 == "default")
            {
                if (para.input.vdw_method == "d2")
                {
                    para.input.vdw_s6 = "0.75";
                }
                else if (para.input.vdw_method == "d3_0" || para.input.vdw_method == "d3_bj")
                {
                    para.input.vdw_s6 = "1.0";
                }
            }
        };
        read_sync_string(input.vdw_s6);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_s8");
        item.annotation = "scale parameter of d3_0/d3_bj";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.vdw_s8 == "default")
            {
                if (para.input.vdw_method == "d3_0")
                {
                    para.input.vdw_s8 = "0.722";
                }
                else if (para.input.vdw_method == "d3_bj")
                {
                    para.input.vdw_s8 = "0.7875";
                }
            }
        };
        read_sync_string(input.vdw_s8);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_a1");
        item.annotation = "damping parameter of d3_0/d3_bj";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.vdw_a1 == "default")
            {
                if (para.input.vdw_method == "d3_0")
                {
                    para.input.vdw_a1 = "1.217";
                }
                else if (para.input.vdw_method == "d3_bj")
                {
                    para.input.vdw_a1 = "0.4289";
                }
            }
        };
        read_sync_string(input.vdw_a1);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_a2");
        item.annotation = "damping parameter of d3_bj";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.vdw_a2 == "default")
            {
                if (para.input.vdw_method == "d3_0")
                {
                    para.input.vdw_a2 = "1.0";
                }
                else if (para.input.vdw_method == "d3_bj")
                {
                    para.input.vdw_a2 = "4.4407";
                }
            }
        };
        read_sync_string(input.vdw_a2);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_d");
        item.annotation = "damping parameter of d2";
        read_sync_double(input.vdw_d);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_abc");
        item.annotation = "third-order term?";
        read_sync_bool(input.vdw_abc);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_c6_file");
        item.annotation = "filename of C6";
        read_sync_string(input.vdw_C6_file);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_c6_unit");
        item.annotation = "unit of C6, Jnm6/mol or eVA6";
        read_sync_string(input.vdw_C6_unit);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if ((para.input.vdw_C6_unit != "Jnm6/mol") && (para.input.vdw_C6_unit != "eVA6"))
            {
                ModuleBase::WARNING_QUIT("ReadInput", "vdw_C6_unit must be Jnm6/mol or eVA6");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("vdw_r0_file");
        item.annotation = "filename of R0";
        read_sync_string(input.vdw_R0_file);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_r0_unit");
        item.annotation = "unit of R0, A or Bohr";
        read_sync_string(input.vdw_R0_unit);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if ((para.input.vdw_R0_unit != "A") && (para.input.vdw_R0_unit != "Bohr"))
            {
                ModuleBase::WARNING_QUIT("ReadInput", "vdw_R0_unit must be A or Bohr");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("vdw_cutoff_type");
        item.annotation = "expression model of periodic structure, radius or period";
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.vdw_cutoff_type != "radius" && para.input.vdw_cutoff_type != "period")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "vdw_cutoff_type must be radius or period");
            }
        };
        read_sync_string(input.vdw_cutoff_type);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_cutoff_radius");
        item.annotation = "radius cutoff for periodic structure";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.vdw_cutoff_radius == "default")
            {
                if (para.input.vdw_method == "d2")
                {
                    para.input.vdw_cutoff_radius = "56.6918";
                }
                else if (para.input.vdw_method == "d3_0" || para.input.vdw_method == "d3_bj")
                {
                    para.input.vdw_cutoff_radius = "95";
                }
                else
                {
                    para.input.vdw_cutoff_radius = "0";
                }
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (std::stod(para.input.vdw_cutoff_radius) <= 0 && para.input.vdw_method != "none")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "vdw_cutoff_radius <= 0 is not allowd");
            }
        };
        read_sync_string(input.vdw_cutoff_radius);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_radius_unit");
        item.annotation = "unit of radius cutoff for periodic structure";
        read_sync_string(input.vdw_radius_unit);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if ((para.input.vdw_radius_unit != "A") && (para.input.vdw_radius_unit != "Bohr"))
            {
                ModuleBase::WARNING_QUIT("ReadInput", "vdw_radius_unit must be A or Bohr");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("vdw_cn_thr");
        item.annotation = "radius cutoff for cn";
        read_sync_double(input.vdw_cn_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.vdw_cn_thr <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "vdw_cn_thr <= 0 is not allowd");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("vdw_cn_thr_unit");
        item.annotation = "unit of cn_thr, Bohr or Angstrom";
        read_sync_string(input.vdw_cn_thr_unit);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if ((para.input.vdw_cn_thr_unit != "A") && (para.input.vdw_cn_thr_unit != "Bohr"))
            {
                ModuleBase::WARNING_QUIT("ReadInput", "vdw_cn_thr_unit must be A or Bohr");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("vdw_cutoff_period");
        item.annotation = "periods of periodic structure";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            if (count == 3)
            {
                para.input.vdw_cutoff_period[0] = std::stoi(item.str_values[0]);
                para.input.vdw_cutoff_period[1] = std::stoi(item.str_values[1]);
                para.input.vdw_cutoff_period[2] = std::stoi(item.str_values[2]);
            }
            else
            {
                ModuleBase::WARNING_QUIT("ReadInput", "vdw_cutoff_period should have 3 values");
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.vdw_cutoff_period[0] <= 0 || para.input.vdw_cutoff_period[1] <= 0
                || para.input.vdw_cutoff_period[2] <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "vdw_cutoff_period should be positive");
            }
        };
        item.get_final_value = [](Input_Item& item, const Parameter& para) {
            item.final_value << para.input.vdw_cutoff_period[0] << " " << para.input.vdw_cutoff_period[1] << " "
                             << para.input.vdw_cutoff_period[2];
        };
#ifdef __MPI
        bcastfuncs.push_back(
            [](Parameter& para) { Parallel_Common::bcast_int((int*)&para.input.vdw_cutoff_period, 3); });
#endif
        this->add_item(item);
    }

    // 14. exx
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
                if (dft_functional_lower == "hf") {
                    para.input.exx_hybrid_alpha = "1";
                } else if (dft_functional_lower == "pbe0" || dft_functional_lower == "hse"
                         || dft_functional_lower == "scan0") {
                    para.input.exx_hybrid_alpha = "0.25";
                } else { // no exx in scf, but will change to non-zero in
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
            {
                if (para.input.gamma_only) {
                    para.input.exx_real_number = "1";
                } else {
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
            {
                std::string& dft_functional = para.input.dft_functional;
                std::string dft_functional_lower = dft_functional;
                std::transform(dft_functional.begin(), dft_functional.end(), dft_functional_lower.begin(), tolower);
                if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0" || dft_functional_lower == "scan0") {
                    para.input.exx_ccp_rmesh_times = "5";
                } else if (dft_functional_lower == "hse") {
                    para.input.exx_ccp_rmesh_times = "1.5";
                } else { // no exx in scf
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

    // 16. tddft
    {
        Input_Item item("td_force_dt");
        item.annotation = "time of force change";
        read_sync_double(input.td_force_dt);
        this->add_item(item);
    }
    {
        Input_Item item("td_vext");
        item.annotation = "add extern potential or not";
        read_sync_bool(input.td_vext);
        this->add_item(item);
    }
    {
        Input_Item item("td_vext_dire");
        item.annotation = "extern potential direction";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_vext_dire = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_vext_dire);
        this->add_item(item);
    }
    {
        Input_Item item("out_dipole");
        item.annotation = "output dipole or not";
        read_sync_bool(input.out_dipole);
        this->add_item(item);
    }
    {
        Input_Item item("out_efield");
        item.annotation = "output dipole or not";
        read_sync_bool(input.out_efield);
        this->add_item(item);
    }
    {
        Input_Item item("out_current");
        item.annotation = "output current or not";
        read_sync_bool(input.out_current);
        this->add_item(item);
    }
    {
        Input_Item item("out_current_k");
        item.annotation = "output current for each k";
        read_sync_bool(input.out_current_k);
        this->add_item(item);
    }
    {
        Input_Item item("out_vecpot");
        item.annotation = "output TDDFT vector potential or not";
        read_sync_bool(input.out_vecpot);
        this->add_item(item);
    }
    {
        Input_Item item("init_vecpot_file");
        item.annotation = "init vector potential through file or not";
        read_sync_bool(input.init_vecpot_file);
        this->add_item(item);
    }
    {
        Input_Item item("td_print_eij");
        item.annotation = "print eij or not";
        read_sync_double(input.td_print_eij);
        this->add_item(item);
    }
    {
        Input_Item item("td_edm");
        item.annotation = "the method to calculate the energy density matrix";
        read_sync_int(input.td_edm);
        this->add_item(item);
    }
    {
        Input_Item item("td_propagator");
        item.annotation = "method of propagator";
        read_sync_int(input.propagator);
        this->add_item(item);
    }
    {
        Input_Item item("td_stype");
        item.annotation = "type of electric field in space domain";
        read_sync_int(input.td_stype);
        this->add_item(item);
    }
    {
        Input_Item item("td_ttype");
        item.annotation = "type of electric field in time domain";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_ttype = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_ttype);
        this->add_item(item);
    }
    {
        Input_Item item("td_tstart");
        item.annotation = " number of steps where electric field starts";
        read_sync_int(input.td_tstart);
        this->add_item(item);
    }
    {
        Input_Item item("td_tend");
        item.annotation = "number of steps where electric field ends";
        read_sync_int(input.td_tend);
        this->add_item(item);
    }
    {
        Input_Item item("td_lcut1");
        item.annotation = "cut1 of interval in length gauge";
        read_sync_double(input.td_lcut1);
        this->add_item(item);
    }
    {
        Input_Item item("td_lcut2");
        item.annotation = "cut2 of interval in length gauge";
        read_sync_double(input.td_lcut2);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_freq");
        item.annotation = "frequency (freq) of Gauss type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_freq = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_gauss_freq);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_phase");
        item.annotation = "phase of Gauss type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_phase = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_gauss_phase);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_sigma");
        item.annotation = "sigma of Gauss type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_sigma = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_gauss_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_t0");
        item.annotation = "step number of time center (t0) of Gauss type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_t0 = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_gauss_t0);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_amp");
        item.annotation = "amplitude of Gauss type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_amp = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_gauss_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_freq");
        item.annotation = "frequency of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_freq = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trape_freq);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_phase");
        item.annotation = "phase of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_phase = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trape_phase);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t1");
        item.annotation = "t1 of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_t1 = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trape_t1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t2");
        item.annotation = "t2 of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_t2 = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trape_t2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t3");
        item.annotation = "t3 of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_t3 = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trape_t3);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_amp");
        item.annotation = "amplitude of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_amp = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trape_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_freq1");
        item.annotation = "frequency 1 of Trigonometric type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_freq1 = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trigo_freq1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_freq2");
        item.annotation = "frequency 2 of Trigonometric type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_freq2 = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trigo_freq2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_phase1");
        item.annotation = "phase 1 of Trigonometric type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_phase1 = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trigo_phase1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_phase2");
        item.annotation = "phase 2 of Trigonometric type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_phase2 = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trigo_phase2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_amp");
        item.annotation = "amplitude of Trigonometric type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_amp = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_trigo_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_heavi_t0");
        item.annotation = "t0 of Heaviside type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_heavi_t0 = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_heavi_t0);
        this->add_item(item);
    }
    {
        Input_Item item("td_heavi_amp");
        item.annotation = "amplitude of Heaviside type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_heavi_amp = longstring(item.str_values, item.get_size());
        };
        sync_string(input.td_heavi_amp);
        this->add_item(item);
    }
    {
        Input_Item item("ocp");
        item.annotation = "change occupation or not";
        read_sync_bool(input.ocp);
        this->add_item(item);
    }
    {
        Input_Item item("ocp_set");
        item.annotation = "set occupation";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.ocp_set = longstring(item.str_values, item.get_size());
        };
        sync_string(input.ocp_set);
        this->add_item(item);
    }

    // 17. berry_wannier
    {
        Input_Item item("berry_phase");
        item.annotation = "calculate berry phase or not";
        read_sync_bool(input.berry_phase);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.berry_phase)
            {
                if (para.input.basis_type != "pw" && para.input.basis_type != "lcao")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "calculate berry phase, please "
                                             "set basis_type = pw or lcao");
                }
                if (para.input.calculation != "nscf")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "calculate berry phase, please set calculation = nscf");
                }
                if (!(para.input.gdir == 1 || para.input.gdir == 2 || para.input.gdir == 3))
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "calculate berry phase, please set gdir = 1 or 2 or 3");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("gdir");
        item.annotation = "calculate the polarization in the direction of the "
                          "lattice vector";
        read_sync_int(input.gdir);
        this->add_item(item);
    }
    {
        Input_Item item("towannier90");
        item.annotation = "use wannier90 code interface or not";
        read_sync_bool(input.towannier90);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.towannier90)
            {
                if (para.input.calculation != "nscf")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "to use towannier90, please set calculation = nscf");
                }
                if (para.input.nspin == 2)
                {
                    if (para.input.wannier_spin != "up" && para.input.wannier_spin != "down")
                    {
                        ModuleBase::WARNING_QUIT("ReadInput",
                                                 "to use towannier90, please set wannier_spin = up "
                                                 "or down");
                    }
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nnkpfile");
        item.annotation = "the wannier90 code nnkp file name";
        read_sync_string(input.nnkpfile);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_spin");
        item.annotation = "calculate spin in wannier90 code interface";
        read_sync_string(input.wannier_spin);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_method");
        item.annotation = "different implementation methods under Lcao basis set";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            /*
                       Developer's notes: on the repair of lcao_in_pw

                       lcao_in_pw is a special basis_type, for scf calculation,
                      it follows workflow of pw, but for nscf the toWannier90
                      calculation, the interface is in ESolver_KS_LCAO_elec,
                       therefore lcao_in_pw for towannier90 calculation follows
                      lcao.

                       In the future lcao_in_pw will have its own ESolver.

                       2023/12/22 use new psi_initializer to expand numerical
                      atomic orbitals, ykhuang
                   */
            if (para.input.towannier90 && para.input.basis_type == "lcao_in_pw")
            {
                para.input.wannier_method = 1;
            }
        };
        read_sync_int(input.wannier_method);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_mmn");
        item.annotation = "output .mmn file or not";
        read_sync_bool(input.out_wannier_mmn);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_amn");
        item.annotation = "output .amn file or not";
        read_sync_bool(input.out_wannier_amn);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_unk");
        item.annotation = "output UNK. file or not";
        read_sync_bool(input.out_wannier_unk);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_eig");
        item.annotation = "output .eig file or not";
        read_sync_bool(input.out_wannier_eig);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_wvfn_formatted");
        item.annotation = "output UNK. file in text format or in binary format";
        read_sync_bool(input.out_wannier_wvfn_formatted);
        this->add_item(item);
    }

    // 18. imlicit_solvation
    {
        Input_Item item("imp_sol");
        item.annotation = "calculate implicit solvation correction or not";
        read_sync_bool(input.imp_sol);
        this->add_item(item);
    }
    {
        Input_Item item("eb_k");
        item.annotation = "the relative permittivity of the bulk solvent";
        read_sync_double(input.eb_k);
        this->add_item(item);
    }
    {
        Input_Item item("tau");
        item.annotation = "the effective surface tension parameter";
        read_sync_double(input.tau);
        this->add_item(item);
    }
    {
        Input_Item item("sigma_k");
        item.annotation = "the width of the diffuse cavity";
        read_sync_double(input.sigma_k);
        this->add_item(item);
    }
    {
        Input_Item item("nc_k");
        item.annotation = "the cut-off charge density";
        read_sync_double(input.nc_k);
        this->add_item(item);
    }

    // 19. OFDFT
    {
        Input_Item item("of_kinetic");
        item.annotation = "kinetic energy functional, such as tf, vw, wt";
        read_sync_string(input.of_kinetic);
        this->add_item(item);
    }
    {
        Input_Item item("of_method");
        item.annotation = "optimization method used in OFDFT, including cg1, "
                          "cg2, tn (default)";
        read_sync_string(input.of_method);
        this->add_item(item);
    }
    {
        Input_Item item("of_conv");
        item.annotation = "the convergence criterion, potential, energy (default), or both";
        read_sync_string(input.of_conv);
        this->add_item(item);
    }
    {
        Input_Item item("of_tole");
        item.annotation = "tolerance of the energy change (in Ry) for "
                          "determining the convergence, default=2e-6 Ry";
        read_sync_double(input.of_tole);
        this->add_item(item);
    }
    {
        Input_Item item("of_tolp");
        item.annotation = "tolerance of potential for determining the "
                          "convergence, default=1e-5 in a.u.";
        read_sync_double(input.of_tolp);
        this->add_item(item);
    }
    {
        Input_Item item("of_tf_weight");
        item.annotation = "weight of TF KEDF";
        read_sync_double(input.of_tf_weight);
        this->add_item(item);
    }
    {
        Input_Item item("of_vw_weight");
        item.annotation = "weight of vW KEDF";
        read_sync_double(input.of_vw_weight);
        this->add_item(item);
    }
    {
        Input_Item item("of_wt_alpha");
        item.annotation = "parameter alpha of WT KEDF";
        read_sync_double(input.of_wt_alpha);
        this->add_item(item);
    }
    {
        Input_Item item("of_wt_beta");
        item.annotation = "parameter beta of WT KEDF";
        read_sync_double(input.of_wt_beta);
        this->add_item(item);
    }
    {
        Input_Item item("of_wt_rho0");
        item.annotation = "the average density of system, used in WT KEDF, in Bohr^-3";
        read_sync_double(input.of_wt_rho0);
        this->add_item(item);
    }
    {
        Input_Item item("of_hold_rho0");
        item.annotation = "If set to 1, the rho0 will be fixed even if the "
                          "volume of system has changed, it will be "
                          "set to 1 automaticly if of_wt_rho0 is not zero";
        read_sync_bool(input.of_hold_rho0);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.of_wt_rho0 != 0)
            {
                para.input.of_hold_rho0 = true; // sunliang add 2022-06-17
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("of_lkt_a");
        item.annotation = "parameter a of LKT KEDF";
        read_sync_double(input.of_lkt_a);
        this->add_item(item);
    }
    {
        Input_Item item("of_full_pw");
        item.annotation = "If set to 1, ecut will be ignored when collect "
                          "planewaves, so that all planewaves will be used";
        read_sync_bool(input.of_full_pw);
        this->add_item(item);
    }
    {
        Input_Item item("of_full_pw_dim");
        item.annotation = "If of_full_pw = true, dimention of FFT is "
                          "testricted to be (0) either odd or even; (1) odd "
                          "only; (2) even only";
        read_sync_int(input.of_full_pw_dim);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (!para.input.of_full_pw)
            {
                para.input.of_full_pw_dim = 0; // sunliang add 2022-08-31
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("of_read_kernel");
        item.annotation = "If set to 1, the kernel of WT KEDF will be filled "
                          "from file of_kernel_file, not from "
                          "formula. Only usable for WT KEDF";
        read_sync_bool(input.of_read_kernel);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.of_kinetic != "wt")
            {
                para.input.of_read_kernel = false; // sunliang add 2022-09-12
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("of_kernel_file");
        item.annotation = "The name of WT kernel file.";
        read_sync_string(input.of_kernel_file);
        this->add_item(item);
    }

    // 20. dft+u
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
            if (para.input.dft_plus_u == 1 && para.input.onsite_radius == 0.0)
            {
                // autoset onsite_radius to 5.0 as default
                para.input.onsite_radius = 5.0;
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
            if (!item.is_read()) {
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
            if (!item.is_read()) {
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

    // 21. spherical bessel
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
            if (count > 0) {
                para.sys.bessel_nao_rcut = para.input.bessel_nao_rcuts[0]; // also compatible with
}
                                                                                 // old input file
            para.sys.nrcut = count;
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.sys.bessel_nao_rcut < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "bessel_nao_rcut must >= 0");
            }
        };
        sync_doublevec(input.bessel_nao_rcuts, para.sys.nrcut, 0.0);
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
    {
        Input_Item item("bessel_descriptor_lmax");
        item.annotation = "lmax used in generating spherical bessel functions";
        read_sync_int(input.bessel_descriptor_lmax);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_ecut");
        item.annotation = "energy cutoff for spherical bessel functions(Ry)";
        read_sync_string(input.bessel_descriptor_ecut);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.bessel_descriptor_ecut == "default")
            {
                para.input.bessel_descriptor_ecut = std::to_string(para.input.ecutwfc);
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (std::stod(para.input.bessel_descriptor_ecut) < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "bessel_descriptor_ecut must >= 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_tolerence");
        item.annotation = "tolerence for spherical bessel root";
        read_sync_double(input.bessel_descriptor_tolerence);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_rcut");
        item.annotation = "radial cutoff for spherical bessel functions(a.u.)";
        read_sync_double(input.bessel_descriptor_rcut);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.bessel_descriptor_rcut < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "bessel_descriptor_rcut must >= 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_smooth");
        item.annotation = "spherical bessel smooth or not";
        read_sync_bool(input.bessel_descriptor_smooth);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_sigma");
        item.annotation = "sphereical bessel smearing_sigma";
        read_sync_double(input.bessel_descriptor_sigma);
        this->add_item(item);
    }

    // 22. non-collinear spin-constrained
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
        Input_Item item("sc_file");
        item.annotation = "file name for parameters used in non-collinear "
                          "spin-constrained DFT (json format)";
        read_sync_string(input.sc_file);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_mag_switch)
            {
                const std::string ss = "test -f " + para.input.sc_file;
                if (system(ss.c_str()))
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "sc_file does not exist");
                }
            }
        };
        this->add_item(item);
    }

    // 23. Quasiatomic Orbital analysis
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
                    if (para.input.qo_basis == "hydrogen") {
                        default_strategy = "energy-valence";
                    } else if ((para.input.qo_basis == "pswfc") || (para.input.qo_basis == "szv")) {
                        default_strategy = "all";
                    } else
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
            if (!item.is_read()) {
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

    // 24. PEXSI
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

    // 25. Linear Response
    {
        Input_Item item("lr_nstates");
        item.annotation = "the number of 2-particle states to be solved";
        read_sync_int(input.lr_nstates);
        this->add_item(item);
    }
    {
        Input_Item item("nocc");
        item.annotation = "the number of occupied orbitals to form the 2-particle basis ( <= nelec/2)";
        read_sync_int(input.nocc);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            const int nocc_default = std::max(static_cast<int>(para.input.nelec + 1) / 2, para.input.nbands);
            if (para.input.nocc <= 0 || para.input.nocc > nocc_default) { para.input.nocc = nocc_default; }
            };
        this->add_item(item);
    }
    {
        Input_Item item("nvirt");
        item.annotation = "the number of virtual orbitals to form the 2-particle basis (nocc + nvirt <= nbands)";
        read_sync_int(input.nvirt);
        this->add_item(item);
    }
    {
        Input_Item item("xc_kernel");
        item.annotation = "exchange correlation (XC) kernel for LR-TDDFT";
        read_sync_string(input.xc_kernel);
        this->add_item(item);
    }
    {
        Input_Item item("lr_solver");
        item.annotation = "the eigensolver for LR-TDDFT";
        read_sync_string(input.lr_solver);
        this->add_item(item);
    }
    {
        Input_Item item("lr_thr");
        item.annotation = "convergence threshold of the LR-TDDFT eigensolver";
        read_sync_double(input.lr_thr);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_lr");
        item.annotation = "whether to output the eigenvectors (excitation amplitudes) in the particle-hole basis";
        read_sync_bool(input.out_wfc_lr);
        this->add_item(item);
    }
    {
        Input_Item item("abs_wavelen_range");
        item.annotation = "the range of wavelength(nm) to output the absorption spectrum ";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.abs_wavelen_range.push_back(std::stod(item.str_values[i]));
            }
            };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            auto& awr = para.input.abs_wavelen_range;
            if (awr.size() < 2) { ModuleBase::WARNING_QUIT("ReadInput", "abs_wavelen_range must have two values"); }
            };
        sync_doublevec(input.abs_wavelen_range, 2, 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("abs_broadening");
        item.annotation = "the broadening (eta) for LR-TDDFT absorption spectrum";
        read_sync_double(input.abs_broadening);
        this->add_item(item);
    }
}
} // namespace ModuleIO