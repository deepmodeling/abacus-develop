
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
        item.check_value = [](const Input_Item& item, const Parameter& para) {
#ifndef __MLKEDF
            if (para.input.of_kinetic == "ml" || para.input.of_kinetic == "mpn" || para.input.of_kinetic == "cpn5")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ML KEDF is not supported.");
            }
#endif
            if (para.input.of_kinetic != "tf" && para.input.of_kinetic != "vw" && para.input.of_kinetic != "wt"
                && para.input.of_kinetic != "lkt" && para.input.of_kinetic != "tf+" 
                && para.input.of_kinetic != "ml" && para.input.of_kinetic != "mpn" && para.input.of_kinetic != "cpn5")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "of_kinetic must be tf, vw, tf+, wt, lkt, ml, mpn, or cpn5");
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            // Set the default parameters for MPN or CPN5 KEDF
            if (para.input.of_kinetic == "mpn")
            {
                para.input.of_kinetic = "ml";

                para.input.of_ml_feg = 3;
                para.input.of_ml_nkernel = 1;
                para.input.of_ml_kernel = {1};
                para.input.of_ml_kernel_scaling = {1.0};
                para.input.of_ml_yukawa_alpha = {1.0};
                para.input.of_ml_gamma = false;
                para.input.of_ml_p = false;
                para.input.of_ml_q = false;
                para.input.of_ml_tanhp = true;
                para.input.of_ml_tanhq = false;
                para.input.of_ml_chi_p = 0.2;
                para.input.of_ml_chi_q = 0.1;
                para.input.of_ml_gammanl = {0};
                para.input.of_ml_pnl = {0};
                para.input.of_ml_qnl = {0};
                para.input.of_ml_xi = {0};
                para.input.of_ml_tanhxi = {1};
                para.input.of_ml_tanhxi_nl = {1};
                para.input.of_ml_tanh_pnl = {0};
                para.input.of_ml_tanh_qnl = {0};
                para.input.of_ml_tanhp_nl = {1};
                para.input.of_ml_tanhq_nl = {0};
                para.input.of_ml_chi_xi = {1.0};
                para.input.of_ml_chi_pnl = {0.2};
                para.input.of_ml_chi_qnl = {0.1};
            }

            if (para.input.of_kinetic == "cpn5")
            {
                para.input.of_kinetic = "ml";

                para.input.of_ml_feg = 3;
                para.input.of_ml_nkernel = 5;
                para.input.of_ml_kernel = {1, 1, 1, 1, 1};
                para.input.of_ml_kernel_scaling = {2.0, 1.5, 1.0, 0.75, 0.5};
                para.input.of_ml_yukawa_alpha = {1.0, 1.0, 1.0, 1.0, 1.0};
                para.input.of_ml_gamma = false;
                para.input.of_ml_p = false;
                para.input.of_ml_q = false;
                para.input.of_ml_tanhp = true;
                para.input.of_ml_tanhq = false;
                para.input.of_ml_chi_p = 0.2;
                para.input.of_ml_chi_q = 0.1;
                para.input.of_ml_gammanl = {0, 0, 0, 0, 0};
                para.input.of_ml_pnl = {0, 0, 0, 0, 0};
                para.input.of_ml_qnl = {0, 0, 0, 0, 0};
                para.input.of_ml_xi = {0, 0, 0, 0, 0};
                para.input.of_ml_tanhxi = {1, 1, 1, 1, 1};
                para.input.of_ml_tanhxi_nl = {1, 1, 1, 1, 1};
                para.input.of_ml_tanh_pnl = {0, 0, 0, 0, 0};
                para.input.of_ml_tanh_qnl = {0, 0, 0, 0, 0};
                para.input.of_ml_tanhp_nl = {1, 1, 1, 1, 1};
                para.input.of_ml_tanhq_nl = {0, 0, 0, 0, 0};
                para.input.of_ml_chi_xi = {0.6, 0.8, 1.0, 1.5, 3.0};
                para.input.of_ml_chi_pnl = {0.2, 0.2, 0.2, 0.2, 0.2};
                para.input.of_ml_chi_qnl = {0.1, 0.1, 0.1, 0.1, 0.1};
            }
        };
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
    {
        Input_Item item("of_ml_gene_data");
        item.annotation = "Generate training data or not";
        read_sync_bool(input.of_ml_gene_data);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_device");
        item.annotation = "Run NN on GPU or CPU";
        read_sync_string(input.of_ml_device);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_feg");
        item.annotation = "The Free Electron Gas limit: 0: no, 3: yes";
        read_sync_int(input.of_ml_feg);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_nkernel");
        item.annotation = "Number of kernels";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.of_ml_nkernel > 0)
            {
                reset_vector(para.input.of_ml_gammanl, para.input.of_ml_nkernel, 0);
                reset_vector(para.input.of_ml_pnl, para.input.of_ml_nkernel, 0);
                reset_vector(para.input.of_ml_qnl, para.input.of_ml_nkernel, 0);
                reset_vector(para.input.of_ml_xi, para.input.of_ml_nkernel, 0);
                reset_vector(para.input.of_ml_tanhxi, para.input.of_ml_nkernel, 0);
                reset_vector(para.input.of_ml_tanhxi_nl, para.input.of_ml_nkernel, 0);
                reset_vector(para.input.of_ml_tanh_pnl, para.input.of_ml_nkernel, 0);
                reset_vector(para.input.of_ml_tanh_qnl, para.input.of_ml_nkernel, 0);
                reset_vector(para.input.of_ml_tanhp_nl, para.input.of_ml_nkernel, 0);
                reset_vector(para.input.of_ml_tanhq_nl, para.input.of_ml_nkernel, 0);
                reset_vector(para.input.of_ml_chi_xi, para.input.of_ml_nkernel, 1.0);
                reset_vector(para.input.of_ml_chi_pnl, para.input.of_ml_nkernel, 1.0);
                reset_vector(para.input.of_ml_chi_qnl, para.input.of_ml_nkernel, 1.0);
                reset_vector(para.input.of_ml_kernel, para.input.of_ml_nkernel, 1);
                reset_vector(para.input.of_ml_kernel_scaling, para.input.of_ml_nkernel, 1.0);
                reset_vector(para.input.of_ml_yukawa_alpha, para.input.of_ml_nkernel, 1.0);
                std::string none = "none";
                reset_vector(para.input.of_ml_kernel_file, para.input.of_ml_nkernel, none);
            }
        };
        read_sync_int(input.of_ml_nkernel);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_kernel");
        item.annotation = "Type of kernel, 1 for wt, 2 for yukawa, and 3 for TKK";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_kernel);
        };
        sync_intvec(input.of_ml_kernel, para.input.of_ml_kernel.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_kernel_scaling");
        item.annotation = "Scaling parameter of kernel, w(r-r') = scaling^3 * w(scaling (r-r'))";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_kernel_scaling);
        };
        sync_doublevec(input.of_ml_kernel_scaling, para.input.of_ml_kernel_scaling.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_yukawa_alpha");
        item.annotation = "Parameter alpha of yukawa kernel";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_yukawa_alpha);
        };
        sync_doublevec(input.of_ml_yukawa_alpha, para.input.of_ml_yukawa_alpha.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_kernel_file");
        item.annotation = "The file of TKK";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            { 
                para.input.of_ml_kernel_file.push_back(item.str_values[i]);
            }
        };
        sync_stringvec(input.of_ml_kernel_file, para.input.of_ml_kernel_file.size(), "");
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_gamma");
        item.annotation = "Descriptor: gamma = (rho / rho0)^(1/3)";
        read_sync_bool(input.of_ml_gamma);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_p");
        item.annotation = "Descriptor: p = |nabla rho|^2 / [2 (3 pi^2)^(1/3) rho^(4/3)]^2";
        read_sync_bool(input.of_ml_p);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_q");
        item.annotation = "Descriptor: q = nabla^2 rho / [4 (3 pi^2)^(2/3) rho^(5/3)]";
        read_sync_bool(input.of_ml_q);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhp");
        item.annotation = "Descriptor: tanhp = tanh(chi_p * p)";
        read_sync_bool(input.of_ml_tanhp);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhq");
        item.annotation = "Descriptor: tanhq = tanh(chi_q * q)";
        read_sync_bool(input.of_ml_tanhq);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_chi_p");
        item.annotation = "Hyperparameter: tanhp = tanh(chi_p * p)";
        read_sync_double(input.of_ml_chi_p);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_chi_q");
        item.annotation = "Hyperparameter: tanhq = tanh(chi_q * q)";
        read_sync_double(input.of_ml_chi_q);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_gammanl");
        item.annotation = "Descriptor: gammanl = int{gamma(r') * w(r-r') dr'}";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_gammanl);
        };
        sync_intvec(input.of_ml_gammanl, para.input.of_ml_gammanl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_pnl");
        item.annotation = "Descriptor: pnl = int{p(r') * w(r-r') dr'}";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_pnl);
        };
        sync_intvec(input.of_ml_pnl, para.input.of_ml_pnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_qnl");
        item.annotation = "Descriptor: qnl = int{q(r') * w(r-r') dr'}";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_qnl);
        };
        sync_intvec(input.of_ml_qnl, para.input.of_ml_qnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_xi");
        item.annotation = "Descriptor: xi = int{rho(r')^(1/3) * w(r-r') dr'} / rho^(1/3)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_xi);
        };
        sync_intvec(input.of_ml_xi, para.input.of_ml_xi.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhxi");
        item.annotation = "Descriptor: tanhxi = tanh(chi_xi * xi)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanhxi);
        };
        sync_intvec(input.of_ml_tanhxi, para.input.of_ml_tanhxi.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhxi_nl");
        item.annotation = "Descriptor: tanhxi_nl = int{tanhxi(r') * w(r-r') dr'}";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanhxi_nl);
        };
        sync_intvec(input.of_ml_tanhxi_nl, para.input.of_ml_tanhxi_nl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanh_pnl");
        item.annotation = "Descriptor: tanh_pnl = tanh(chi_pnl * pnl)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanh_pnl);
        };
        sync_intvec(input.of_ml_tanh_pnl, para.input.of_ml_tanh_pnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanh_qnl");
        item.annotation = "Descriptor: tanh_qnl = tanh(chi_qnl * qnl)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanh_qnl);
        };
        sync_intvec(input.of_ml_tanh_qnl, para.input.of_ml_tanh_qnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhp_nl");
        item.annotation = "Descriptor: tanhp_nl = int{tanhp(r') * w(r-r') dr'}";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanhp_nl);
        };
        sync_intvec(input.of_ml_tanhp_nl, para.input.of_ml_tanhp_nl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhq_nl");
        item.annotation = "Descriptor: tanhq_nl = int{tanhq(r') * w(r-r') dr'}";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanhq_nl);
        };
        sync_intvec(input.of_ml_tanhq_nl, para.input.of_ml_tanhq_nl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_chi_xi");
        item.annotation = "Hyperparameter: tanhpxi = tanh(chi_xi * xi)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_chi_xi);
        };
        sync_doublevec(input.of_ml_chi_xi, para.input.of_ml_chi_xi.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_chi_pnl");
        item.annotation = "Hyperparameter: tanh_pnl = tanh(chi_pnl * pnl)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_chi_pnl);
        };
        sync_doublevec(input.of_ml_chi_pnl, para.input.of_ml_chi_pnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_chi_qnl");
        item.annotation = "Hyperparameter: tanh_qnl = tanh(chi_qnl * qnl)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_chi_qnl);
        };
        sync_doublevec(input.of_ml_chi_qnl, para.input.of_ml_chi_qnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_local_test");
        item.annotation = "Test: read in the density, and output the F and Pauli potential";
        read_sync_bool(input.of_ml_local_test);
        this->add_item(item);
    }
}
} // namespace ModuleIO