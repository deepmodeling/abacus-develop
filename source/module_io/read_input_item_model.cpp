#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_model()
{
    // Electric field and dipole correction
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

    // Gate field
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

    // imlicit_solvation
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

    // vdW Correction
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
                // else if (para.input.vdw_method == "d3_0" || para.input.vdw_method == "d3_bj")
                // {
                //     para.input.vdw_s6 = "1.0";
                // }
            }
        };
        read_sync_string(input.vdw_s6);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_s8");
        item.annotation = "scale parameter of d3_0/d3_bj";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            // if (para.input.vdw_s8 == "default")
            // {
            //     if (para.input.vdw_method == "d3_0")
            //     {
            //         para.input.vdw_s8 = "0.722";
            //     }
            //     else if (para.input.vdw_method == "d3_bj")
            //     {
            //         para.input.vdw_s8 = "0.7875";
            //     }
            // }
        };
        read_sync_string(input.vdw_s8);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_a1");
        item.annotation = "damping parameter of d3_0/d3_bj";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            // if (para.input.vdw_a1 == "default")
            // {
            //     if (para.input.vdw_method == "d3_0")
            //     {
            //         para.input.vdw_a1 = "1.217";
            //     }
            //     else if (para.input.vdw_method == "d3_bj")
            //     {
            //         para.input.vdw_a1 = "0.4289";
            //     }
            // }
        };
        read_sync_string(input.vdw_a1);
        this->add_item(item);
    }
    {
        Input_Item item("vdw_a2");
        item.annotation = "damping parameter of d3_bj";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            // if (para.input.vdw_a2 == "default")
            // {
            //     if (para.input.vdw_method == "d3_0")
            //     {
            //         para.input.vdw_a2 = "1.0";
            //     }
            //     else if (para.input.vdw_method == "d3_bj")
            //     {
            //         para.input.vdw_a2 = "4.4407";
            //     }
            // }
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
}
} // namespace ModuleIO