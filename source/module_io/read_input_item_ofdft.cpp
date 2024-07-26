
#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_ofdft()
{
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
}
} // namespace ModuleIO