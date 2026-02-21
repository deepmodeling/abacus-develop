
#include "source_base/global_function.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_ofdft()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("of_kinetic");
        item.annotation = "kinetic energy functional, such as tf, vw, wt";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "String";
        item.description = R"(Kinetic energy functional type:
* tf: Thomas-Fermi (TF) functional
* vw: von Weizsacker (vW) functional
* tf+: TF + vW functional
* wt: Wang-Teter (WT) functional
* xwm: XWM functional
* lkt: Luo-Karasiev-Trickey (LKT) functional
* ml: Machine learning KEDF
* mpn: MPN KEDF (automatically sets ml parameters)
* cpn5: CPN5 KEDF (automatically sets ml parameters))";
        item.default_value = "wt";
        item.unit = "";
        item.availability = "OFDFT";
        item.check_value = [](const Input_Item& item, const Parameter& para) {
#ifndef __MLALGO
            if (para.input.of_kinetic == "ml" || para.input.of_kinetic == "mpn" || para.input.of_kinetic == "cpn5")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "Error: ML KEDF requires ENABLE_MLALGO option.\n "
                                                      "Please enable ENABLE_MLALGO during compilation to use this feature.");
            }
#endif
            if (para.input.of_kinetic != "tf" && para.input.of_kinetic != "vw" && para.input.of_kinetic != "wt"
                && para.input.of_kinetic != "xwm" && para.input.of_kinetic != "lkt" && para.input.of_kinetic != "tf+" 
                && para.input.of_kinetic != "ml" && para.input.of_kinetic != "mpn" && para.input.of_kinetic != "cpn5")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "of_kinetic must be tf, vw, tf+, wt, xwm, lkt, ml, mpn, or cpn5");
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
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "String";
        item.description = R"(The optimization method used in OFDFT.
* cg1: Polak-Ribiere. Standard CG algorithm.
* cg2: Hager-Zhang (generally faster than cg1).
* tn: Truncated Newton algorithm.)";
        item.default_value = "tn";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_string(input.of_method);
        this->add_item(item);
    }
    {
        Input_Item item("of_conv");
        item.annotation = "the convergence criterion, potential, energy (default), or both";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "String";
        item.description = R"(Criterion used to check the convergence of OFDFT.
* energy: Total energy changes less than of_tole.
* potential: The norm of potential is less than of_tolp.
* both: Both energy and potential must satisfy the convergence criterion.)";
        item.default_value = "energy";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_string(input.of_conv);
        this->add_item(item);
    }
    {
        Input_Item item("of_tole");
        item.annotation = "tolerance of the energy change (in Ry) for "
                          "determining the convergence, default=2e-6 Ry";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Real";
        item.description = "Tolerance of the energy change for determining the convergence.";
        item.default_value = "2e-6";
        item.unit = "Ry";
        item.availability = "OFDFT";
        read_sync_double(input.of_tole);
        this->add_item(item);
    }
    {
        Input_Item item("of_tolp");
        item.annotation = "tolerance of potential for determining the "
                          "convergence, default=1e-5 in a.u.";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Real";
        item.description = "Tolerance of potential for determining the convergence.";
        item.default_value = "1e-5";
        item.unit = "Ry";
        item.availability = "OFDFT";
        read_sync_double(input.of_tolp);
        this->add_item(item);
    }
    {
        Input_Item item("of_tf_weight");
        item.annotation = "weight of TF KEDF";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Real";
        item.description = "Weight of TF KEDF (kinetic energy density functional).";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "OFDFT with of_kinetic=tf, tf+, wt, xwm";
        read_sync_double(input.of_tf_weight);
        this->add_item(item);
    }
    {
        Input_Item item("of_vw_weight");
        item.annotation = "weight of vW KEDF";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Real";
        item.description = "Weight of vW KEDF (kinetic energy density functional).";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "OFDFT with of_kinetic=vw, tf+, wt, lkt, xwm";
        read_sync_double(input.of_vw_weight);
        this->add_item(item);
    }
    {
        Input_Item item("of_wt_alpha");
        item.annotation = "parameter alpha of WT KEDF";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Real";
        item.description = "Parameter alpha of WT KEDF (kinetic energy density functional).";
        item.default_value = "";
        item.unit = "";
        item.availability = "OFDFT with of_kinetic=wt";
        read_sync_double(input.of_wt_alpha);
        this->add_item(item);
    }
    {
        Input_Item item("of_wt_beta");
        item.annotation = "parameter beta of WT KEDF";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Real";
        item.description = "Parameter beta of WT KEDF (kinetic energy density functional).";
        item.default_value = "";
        item.unit = "";
        item.availability = "OFDFT with of_kinetic=wt";
        read_sync_double(input.of_wt_beta);
        this->add_item(item);
    }
    {
        Input_Item item("of_wt_rho0");
        item.annotation = "the average density of system, used in WT KEDF, in Bohr^-3";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Real";
        item.description = "The average density of system.";
        item.default_value = "0.0";
        item.unit = "Bohr^-3";
        item.availability = "OFDFT with of_kinetic=wt";
        read_sync_double(input.of_wt_rho0);
        this->add_item(item);
    }
    {
        Input_Item item("of_hold_rho0");
        item.annotation = "If set to 1, the rho0 will be fixed even if the "
                          "volume of system has changed, it will be "
                          "set to 1 automaticly if of_wt_rho0 is not zero";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Boolean";
        item.description = R"(Whether to fix the average density rho0.
* True: rho0 will be fixed even if the volume of system has changed, it will be set to True automatically if of_wt_rho0 is not zero.
* False: rho0 will change if volume of system has changed.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "OFDFT with of_kinetic=wt";
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
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Real";
        item.description = "Parameter a of LKT KEDF (kinetic energy density functional).";
        item.default_value = "1.3";
        item.unit = "";
        item.availability = "OFDFT with of_kinetic=lkt";
        read_sync_double(input.of_lkt_a);
        this->add_item(item);
    }
    {
        Input_Item item("of_xwm_rho_ref");
        item.annotation = "The reference density of XWM KEDF";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Real";
        item.description = "Reference charge density for XWM kinetic energy functional. If set to 0, the program will use average charge density.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "OFDFT with of_kinetic=xwm";
        read_sync_double(input.of_xwm_rho_ref);
        this->add_item(item);
    }
    {
        Input_Item item("of_xwm_kappa");
        item.annotation = "The parameter kappa of XWM KEDF";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Real";
        item.description = "Parameter for XWM kinetic energy functional. See PHYSICAL REVIEW B 100, 205132 (2019) for optimal values.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "OFDFT with of_kinetic=xwm";
        read_sync_double(input.of_xwm_kappa);
        this->add_item(item);
    }
    {
        Input_Item item("of_read_kernel");
        item.annotation = "If set to 1, the kernel of WT KEDF will be filled "
                          "from file of_kernel_file, not from "
                          "formula. Only usable for WT KEDF";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Boolean";
        item.description = R"(Whether to read in the kernel file.
* True: The kernel of WT KEDF (kinetic energy density functional) will be filled from the file specified by of_kernel_file.
* False: The kernel of WT KEDF (kinetic energy density functional) will be filled from formula.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "OFDFT with of_kinetic=wt";
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
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "String";
        item.description = "The name of WT kernel file.";
        item.default_value = "WTkernel.txt";
        item.unit = "";
        item.availability = "OFDFT with of_read_kernel=True";
        read_sync_string(input.of_kernel_file);
        this->add_item(item);
    }
    {
        Input_Item item("of_full_pw");
        item.annotation = "If set to 1, ecut will be ignored when collect "
                          "planewaves, so that all planewaves will be used";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Boolean";
        item.description = R"(Whether to use full planewaves.
* True: Ecut will be ignored while collecting planewaves, so that all planewaves will be used in FFT.
* False: Only use the planewaves inside ecut, the same as KSDFT.)";
        item.default_value = "True";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_bool(input.of_full_pw);
        this->add_item(item);
    }
    {
        Input_Item item("of_full_pw_dim");
        item.annotation = "If of_full_pw = true, dimention of FFT is "
                          "testricted to be (0) either odd or even; (1) odd "
                          "only; (2) even only";
        item.category = "OFDFT: orbital free density functional theory";
        item.type = "Integer";
        item.description = R"(Specify the parity of FFT dimensions.
* 0: either odd or even.
* 1: odd only.
* 2: even only.

Note: Even dimensions may cause slight errors in FFT. It should be ignorable in ofdft calculation, but it may make Cardinal B-spline interpolation unstable, so please set of_full_pw_dim = 1 if nbspline != -1.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT with of_full_pw = True";
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
        Input_Item item("of_ml_gene_data");
        item.annotation = "Generate training data or not";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Boolean";
        item.description = "Controls the generation of machine learning training data. When enabled, training data in .npy format will be saved in the directory OUT.${suffix}/.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "Used only for KSDFT with plane wave basis";
        read_sync_bool(input.of_ml_gene_data);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_device");
        item.annotation = "Run NN on GPU or CPU";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "String";
        item.description = R"(Run Neural Network on GPU or CPU.
* cpu: CPU
* gpu: GPU)";
        item.default_value = "cpu";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_string(input.of_ml_device);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_feg");
        item.annotation = "The Free Electron Gas limit: 0: no, 3: yes";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Integer";
        item.description = R"(The method to incorporate the Free Electron Gas (FEG) limit.
* 0: Do not incorporate the FEG limit.
* 1: Incorporate the FEG limit by translation.
* 3: Incorporate the FEG limit by nonlinear transformation using softplus function.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_int(input.of_ml_feg);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_nkernel");
        item.annotation = "Number of kernels";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Integer";
        item.description = "Number of kernel functions.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "OFDFT";
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
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = R"(Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the type of the i-th kernel function.
* 1: Wang-Teter kernel function.
* 2: Modified Yukawa function, and alpha is specified by of_ml_yukawa_alpha.
* 3: Truncated kinetic kernel (TKK), the file containing TKK is specified by of_ml_kernel_file.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_kernel);
        };
        sync_intvec(input.of_ml_kernel, para.input.of_ml_kernel.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_kernel_scaling");
        item.annotation = "Scaling parameter of kernel, w(r-r') = scaling^3 * w(scaling (r-r'))";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Real";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the RECIPROCAL of scaling parameter of the i-th kernel function.";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_kernel_scaling);
        };
        sync_doublevec(input.of_ml_kernel_scaling, para.input.of_ml_kernel_scaling.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_yukawa_alpha");
        item.annotation = "Parameter alpha of yukawa kernel";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Real";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the parameter alpha of i-th kernel function. ONLY used for Yukawa kernel function.";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_yukawa_alpha);
        };
        sync_doublevec(input.of_ml_yukawa_alpha, para.input.of_ml_yukawa_alpha.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_kernel_file");
        item.annotation = "The file of TKK";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of String";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the file containing the i-th kernel function. ONLY used for TKK.";
        item.default_value = "none";
        item.unit = "";
        item.availability = "OFDFT";
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
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Boolean";
        item.description = "Local descriptor: gamma = (rho / rho0)^(1/3).";
        item.default_value = "False";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_bool(input.of_ml_gamma);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_p");
        item.annotation = "Descriptor: p = |nabla rho|^2 / [2 (3 pi^2)^(1/3) rho^(4/3)]^2";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Boolean";
        item.description = "Semi-local descriptor: p = |nabla rho|^2 / [2 (3 pi^2)^(1/3) rho^(4/3)]^2.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_bool(input.of_ml_p);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_q");
        item.annotation = "Descriptor: q = nabla^2 rho / [4 (3 pi^2)^(2/3) rho^(5/3)]";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Boolean";
        item.description = "Semi-local descriptor: q = nabla^2 rho / [4 (3 pi^2)^(2/3) rho^(5/3)].";
        item.default_value = "False";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_bool(input.of_ml_q);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhp");
        item.annotation = "Descriptor: tanhp = tanh(chi_p * p)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Boolean";
        item.description = "Semi-local descriptor: tanhp = tanh(chi_p * p).";
        item.default_value = "False";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_bool(input.of_ml_tanhp);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhq");
        item.annotation = "Descriptor: tanhq = tanh(chi_q * q)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Boolean";
        item.description = "Semi-local descriptor: tanhq = tanh(chi_q * q).";
        item.default_value = "False";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_bool(input.of_ml_tanhq);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_chi_p");
        item.annotation = "Hyperparameter: tanhp = tanh(chi_p * p)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Real";
        item.description = "Hyperparameter chi_p: tanhp = tanh(chi_p * p).";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_double(input.of_ml_chi_p);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_chi_q");
        item.annotation = "Hyperparameter: tanhq = tanh(chi_q * q)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Real";
        item.description = "Hyperparameter chi_q: tanhq = tanh(chi_q * q).";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_double(input.of_ml_chi_q);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_gammanl");
        item.annotation = "Descriptor: gammanl = int{gamma(r') * w(r-r') dr'}";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor gammanl defined by the i-th kernel function.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_gammanl);
        };
        sync_intvec(input.of_ml_gammanl, para.input.of_ml_gammanl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_pnl");
        item.annotation = "Descriptor: pnl = int{p(r') * w(r-r') dr'}";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor pnl defined by the i-th kernel function.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_pnl);
        };
        sync_intvec(input.of_ml_pnl, para.input.of_ml_pnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_qnl");
        item.annotation = "Descriptor: qnl = int{q(r') * w(r-r') dr'}";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor qnl defined by the i-th kernel function.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_qnl);
        };
        sync_intvec(input.of_ml_qnl, para.input.of_ml_qnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_xi");
        item.annotation = "Descriptor: xi = int{rho(r')^(1/3) * w(r-r') dr'} / rho^(1/3)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor xi defined by the i-th kernel function.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_xi);
        };
        sync_intvec(input.of_ml_xi, para.input.of_ml_xi.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhxi");
        item.annotation = "Descriptor: tanhxi = tanh(chi_xi * xi)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanhxi defined by the i-th kernel function.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanhxi);
        };
        sync_intvec(input.of_ml_tanhxi, para.input.of_ml_tanhxi.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhxi_nl");
        item.annotation = "Descriptor: tanhxi_nl = int{tanhxi(r') * w(r-r') dr'}";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanhxi_nl defined by the i-th kernel function.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanhxi_nl);
        };
        sync_intvec(input.of_ml_tanhxi_nl, para.input.of_ml_tanhxi_nl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanh_pnl");
        item.annotation = "Descriptor: tanh_pnl = tanh(chi_pnl * pnl)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanh_pnl defined by the i-th kernel function.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanh_pnl);
        };
        sync_intvec(input.of_ml_tanh_pnl, para.input.of_ml_tanh_pnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanh_qnl");
        item.annotation = "Descriptor: tanh_qnl = tanh(chi_qnl * qnl)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanh_qnl defined by the i-th kernel function.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanh_qnl);
        };
        sync_intvec(input.of_ml_tanh_qnl, para.input.of_ml_tanh_qnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhp_nl");
        item.annotation = "Descriptor: tanhp_nl = int{tanhp(r') * w(r-r') dr'}";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanhp_nl defined by the i-th kernel function.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanhp_nl);
        };
        sync_intvec(input.of_ml_tanhp_nl, para.input.of_ml_tanhp_nl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_tanhq_nl");
        item.annotation = "Descriptor: tanhq_nl = int{tanhq(r') * w(r-r') dr'}";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Integer";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanhq_nl defined by the i-th kernel function.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_tanhq_nl);
        };
        sync_intvec(input.of_ml_tanhq_nl, para.input.of_ml_tanhq_nl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_chi_xi");
        item.annotation = "Hyperparameter: tanhpxi = tanh(chi_xi * xi)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Real";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the hyperparameter chi_xi of non-local descriptor tanhxi defined by the i-th kernel function.";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_chi_xi);
        };
        sync_doublevec(input.of_ml_chi_xi, para.input.of_ml_chi_xi.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_chi_pnl");
        item.annotation = "Hyperparameter: tanh_pnl = tanh(chi_pnl * pnl)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Real";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the hyperparameter chi_pnl of non-local descriptor tanh_pnl defined by the i-th kernel function.";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_chi_pnl);
        };
        sync_doublevec(input.of_ml_chi_pnl, para.input.of_ml_chi_pnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_chi_qnl");
        item.annotation = "Hyperparameter: tanh_qnl = tanh(chi_qnl * qnl)";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Vector of Real";
        item.description = "Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the hyperparameter chi_qnl of non-local descriptor tanh_qnl defined by the i-th kernel function.";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "OFDFT";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.of_ml_chi_qnl);
        };
        sync_doublevec(input.of_ml_chi_qnl, para.input.of_ml_chi_qnl.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("of_ml_local_test");
        item.annotation = "Test: read in the density, and output the F and Pauli potential";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Boolean";
        item.description = "FOR TEST. Read in the density, and output the F and Pauli potential.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "OFDFT";
        read_sync_bool(input.of_ml_local_test);
        this->add_item(item);
    }
    {
        Input_Item item("ml_exx");
        item.annotation = "Use ML EXX or not";
        item.category = "ML-KEDF: machine learning based kinetic energy density functional for OFDFT";
        item.type = "Boolean";
        item.description = "Whether to use machine learning based exact exchange (ML-EXX).";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.ml_exx);
        this->add_item(item);
    }
}
} // namespace ModuleIO