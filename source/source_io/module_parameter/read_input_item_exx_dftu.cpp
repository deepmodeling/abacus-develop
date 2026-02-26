#include "source_base/constants.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_exx()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    // EXX
    {
        Input_Item item("exx_fock_alpha");
        item.annotation = "fraction of Fock exchange 1/r in hybrid functionals";
        item.category = "Exact Exchange (Common)";
        item.type = "Real";
        item.description = R"(Fraction of full-ranged Fock exchange $1/r$ in range-separated hybrid functionals.)";
        item.default_value = "see hybrid_func_params";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para)
        {
            para.input.exx_fock_alpha = item.str_values;
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) 
        {
            if (para.input.exx_fock_alpha.size()==1 && para.input.exx_fock_alpha[0]=="default")
            {
                std::string& dft_functional = para.input.dft_functional;
                std::string dft_functional_lower = dft_functional;
                std::transform(dft_functional.begin(), dft_functional.end(), dft_functional_lower.begin(), tolower);
                if (dft_functional_lower == "hf" ||
                    dft_functional_lower == "lc_pbe" || dft_functional_lower == "lc_wpbe" ||
                    dft_functional_lower == "lrc_wpbe" || dft_functional_lower == "lrc_wpbeh" ||
                    dft_functional_lower == "muller" || dft_functional_lower == "power"      // added by jghan 2024-07-06
                    || dft_functional_lower == "wp22" )
                {
                    para.input.exx_fock_alpha = {"1"};
                }
                else if (dft_functional_lower == "pbe0" || dft_functional_lower == "scan0")
                {
                    para.input.exx_fock_alpha = {"0.25"};
                }
                else if (dft_functional_lower == "b3lyp")
                {
                    para.input.exx_fock_alpha = {"0.2"};
                }
                else if (dft_functional_lower == "cam_pbeh")
                {
                    para.input.exx_fock_alpha = {"0.2"};
                }
                else
                {   // no exx in scf, but will change to non-zero in postprocess like rpa
                    para.input.exx_fock_alpha = {};
                }
            }
        };
        sync_stringvec(input.exx_fock_alpha, para.input.exx_fock_alpha.size(), "");
        this->add_item(item);
    }
    {
        Input_Item item("exx_erfc_alpha");
        item.annotation = "fraction of exchange erfc(wr)/r in hybrid functionals";
        item.category = "Exact Exchange (Common)";
        item.type = "Real";
        item.description = R"(Fraction of short-ranged Fock exchange $\mathrm{erfc}(\omega r)/r$ in range-separated hybrid functionals.)";
        item.default_value = "see hybrid_func_params";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para)
        {
            para.input.exx_erfc_alpha = item.str_values;
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) 
        {
            if (para.input.exx_erfc_alpha.size()==1 &&  para.input.exx_erfc_alpha[0]=="default")
            {
                std::string& dft_functional = para.input.dft_functional;
                std::string dft_functional_lower = dft_functional;
                std::transform(dft_functional.begin(), dft_functional.end(), dft_functional_lower.begin(), tolower);
                if (dft_functional_lower == "hse")
                {
                    para.input.exx_erfc_alpha = {"0.25"};
                }
                else if (dft_functional_lower == "lrc_wpbeh")
                {
                    para.input.exx_erfc_alpha = {"-0.8"};
                }
                else if (dft_functional_lower == "cam_pbeh")
                {
                    para.input.exx_erfc_alpha = {"0.8"};
                }
                else if (dft_functional_lower == "cwp22")
                {
                    para.input.exx_erfc_alpha = {"1"};
                }
                else if (dft_functional_lower == "lc_pbe" || dft_functional_lower == "lc_wpbe" ||
                    dft_functional_lower == "lrc_wpbe" || dft_functional_lower == "wp22")
                {
                    para.input.exx_erfc_alpha = {"-1"};
                }
                else
                { // no exx in scf, but will change to non-zero in postprocess like rpa
                    para.input.exx_erfc_alpha = {};
                }
            }
        };
        sync_stringvec(input.exx_erfc_alpha, para.input.exx_erfc_alpha.size(), "");
        this->add_item(item);
    }
    {
        Input_Item item("exx_erfc_omega");
        item.annotation = "range-separation parameter erfc(wr)/r in hybrid functionals";
        item.category = "Exact Exchange (Common)";
        item.type = "Real";
        item.description = R"(Range-separation parameter $\omega$ in the short-ranged Fock term $\mathrm{erfc}(\omega r)/r$.)";
        item.default_value = "see hybrid_func_params";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para)
        {
            para.input.exx_erfc_omega = item.str_values;
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) 
        {
            if (para.input.exx_erfc_omega.size()==1 &&  para.input.exx_erfc_omega[0]=="default")
            {
                std::string& dft_functional = para.input.dft_functional;
                std::string dft_functional_lower = dft_functional;
                std::transform(dft_functional.begin(), dft_functional.end(), dft_functional_lower.begin(), tolower);
                if (dft_functional_lower == "hse" || dft_functional_lower == "cwp22" || dft_functional_lower == "wp22")
                {
                    para.input.exx_erfc_omega = {"0.11"};
                }
                else if (dft_functional_lower == "lc_pbe")
                {
                    para.input.exx_erfc_omega = {"0.33"};
                }
                else if (dft_functional_lower == "lc_wpbe")
                {
                    para.input.exx_erfc_omega = {"0.4"};
                }
                else if (dft_functional_lower == "lrc_wpbe")
                {
                    para.input.exx_erfc_omega = {"0.3"};
                }
                else if (dft_functional_lower == "lrc_wpbeh")
                {
                    para.input.exx_erfc_omega = {"0.2"};
                }
                else if (dft_functional_lower == "cam_pbeh")
                {
                    para.input.exx_erfc_omega = {"0.7"};
                }
                else
                {
                    para.input.exx_erfc_omega = {};
                }
            }
        };
        sync_stringvec(input.exx_erfc_omega, para.input.exx_erfc_omega.size(), "");
        this->add_item(item);
    }
    {
        Input_Item item("exx_separate_loop");
        item.annotation = "if 1, a two-step method is employed, else it will "
                          "start with a GGA-Loop, and then Hybrid-Loop";
        item.category = "Exact Exchange (Common)";
        item.type = "Boolean";
        item.description = R"(There are two types of iterative approaches provided by ABACUS to evaluate Fock exchange.
* False: Start with a GGA-Loop, and then Hybrid-Loop, in which EXX Hamiltonian is updated with electronic iterations.
* True: A two-step method is employed, i.e. in the inner iterations, density matrix is updated, while in the outer iterations, is calculated based on density matrix that converges in the inner iteration.)";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.esolver_type == "tddft" && para.input.exx_separate_loop)
            {
                GlobalV::ofs_running << "For RT-TDDFT with hybrid functionals, only exx_separate_loop = 0 is supported" << std::endl;
                para.input.exx_separate_loop = false;
            }
        };
        read_sync_bool(input.exx_separate_loop);
        this->add_item(item);
    }
    {
        Input_Item item("exx_hybrid_step");
        item.annotation = "the maximal electronic iteration number in the "
                          "evaluation of Fock exchange";
        item.category = "Exact Exchange (Common)";
        item.type = "Integer";
        item.description = "The maximal iteration number of the outer-loop, where the Fock exchange is calculated";
        item.default_value = "100";
        item.unit = "";
        item.availability = "exx_separate_loop==1";
        read_sync_int(input.exx_hybrid_step);
        item.check_value = [](const Input_Item& item, const Parameter& para) 
        {
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
        item.category = "Exact Exchange (Common)";
        item.type = "Real";
        item.description = "Mixing parameter for densty matrix in each iteration of the outer-loop";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "exx_separate_loop==1";
        read_sync_double(input.exx_mixing_beta);
        this->add_item(item);
    }
    {
        Input_Item item("exx_fock_lambda");
        item.annotation = "used to compensate for divergence points at G=0 in "
                          "the evaluation of Fock exchange using "
                          "lcao_in_pw method";
        item.category = "Exact Exchange (LCAO in PW)";
        item.type = "Real";
        item.description = "It is used to compensate for divergence points at G=0 in the evaluation of Fock exchange using lcao_in_pw method.";
        item.default_value = "0.3";
        item.unit = "";
        item.availability = "basis_type==lcao_in_pw";
        item.read_value = [](const Input_Item& item, Parameter& para)
        {
            para.input.exx_fock_lambda = item.str_values;
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) 
        {
            if (para.input.exx_fock_lambda.size()==1 &&  para.input.exx_fock_lambda[0]=="default")
            {
                para.input.exx_fock_lambda = std::vector<std::string>(para.input.exx_fock_alpha.size(), "0.3");
            }
        };
        sync_stringvec(input.exx_fock_lambda, para.input.exx_fock_lambda.size(), "");
        this->add_item(item);
    }
    {
        Input_Item item("exx_pca_threshold");
        item.annotation = "threshold to screen on-site ABFs in exx";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "To accelerate the evaluation of four-center integrals (), the product of atomic orbitals are expanded in the basis of auxiliary basis functions (ABF): . The size of the ABF (i.e. number of ) is reduced using principal component analysis. When a large PCA threshold is used, the number of ABF will be reduced, hence the calculation becomes faster. However, this comes at the cost of computational accuracy. A relatively safe choice of the value is 1e-4.";
        item.default_value = "1E-4";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.exx_pca_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_c_threshold");
        item.annotation = "threshold to screen C matrix in exx";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "See also the entry exx_pca_threshold. Smaller components (less than exx_c_threshold) of the matrix are neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.";
        item.default_value = "1E-4";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.exx_c_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_cs_inv_thr");
        item.annotation = "threshold to inverse Vq in abfs for generating Cs";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "By default, the Coulomb matrix inversion required for obtaining LRI coefficients is performed using LU decomposition. However, this approach may suffer from numerical instabilities when a large set of auxiliary basis functions (ABFs) is employed. When exx_cs_inv_thr > 0, the inversion is instead carried out via matrix diagonalization. Eigenvalues smaller than exx_cs_inv_thr are discarded to improve numerical stability. A relatively safe and commonly recommended value is 1e-5.";
        item.default_value = "-1";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.exx_cs_inv_thr);
        this->add_item(item);
    }
    {
        Input_Item item("exx_v_threshold");
        item.annotation = "threshold to screen C matrix in exx";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "See also the entry exx_pca_threshold. With the approximation , the four-center integral in Fock exchange is expressed as , where is a double-center integral. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.";
        item.default_value = "1E-1";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.exx_v_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_dm_threshold");
        item.annotation = "threshold to screen density matrix in exx";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "The Fock exchange can be expressed as where D is the density matrix. Smaller values of the density matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.";
        item.default_value = "1E-4";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.exx_dm_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_c_grad_threshold");
        item.annotation = "threshold to screen nabla C matrix in exx";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "See also the entry exx_pca_threshold. is used in force. Smaller components (less than exx_c_grad_threshold) of the matrix are neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.";
        item.default_value = "1E-4";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.exx_c_grad_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_v_grad_threshold");
        item.annotation = "threshold to screen nabla V matrix in exx";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "See also the entry exx_pca_threshold. With the approximation , the four-center integral in Fock exchange is expressed as , where is a double-center integral. is used in force. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.";
        item.default_value = "1E-1";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.exx_v_grad_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_c_grad_r_threshold");
        item.annotation = "threshold to screen nabla C * R matrix in exx";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "See also the entry exx_pca_threshold. is used in stress. Smaller components (less than exx_c_grad_r_threshold) of the matrix are neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.";
        item.default_value = "1E-4";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.exx_c_grad_r_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_v_grad_r_threshold");
        item.annotation = "threshold to screen nabla V * R matrix in exx";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "See also the entry exx_pca_threshold. With the approximation , the four-center integral in Fock exchange is expressed as , where is a double-center integral. is used in force and stress. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.";
        item.default_value = "1E-1";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.exx_v_grad_r_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("exx_ccp_rmesh_times");
        item.annotation = "how many times larger the radial mesh required for "
                          "calculating Columb potential is to that "
                          "of atomic orbitals";
        item.category = "Exact Exchange (LCAO)";
        item.type = "String";
        item.description = "This parameter determines how many times larger the radial mesh required for calculating Columb potential is to that of atomic orbitals. The value should be larger than 0. Reducing this value can effectively increase the speed of self-consistent calculations using hybrid functionals.";
        item.default_value = "";
        item.unit = "";
        item.availability = "";
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
                // added by jghan 2024-07-06
                else if (dft_functional_lower == "muller" || dft_functional_lower == "power")
                {
                    para.input.exx_ccp_rmesh_times = "5";
                }
                else if (dft_functional_lower == "wp22")
                {
                    para.input.exx_ccp_rmesh_times = "5";
                    // exx_ccp_rmesh_times = "1.5";
                }
                else if (dft_functional_lower == "cwp22")
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
            if (std::stod(para.input.exx_ccp_rmesh_times) <=0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exx_ccp_rmesh_times must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_opt_orb_lmax");
        item.annotation = "the maximum l of the spherical Bessel functions for opt ABFs";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Integer";
        item.description = "The maximum l of the spherical Bessel functions, when the radial part of opt-ABFs are generated as linear combinations of spherical Bessel functions. A reasonable choice is 2.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "calculation==gen_opt_abfs";
        read_sync_int(input.exx_opt_orb_lmax);
        this->add_item(item);
    }
    {
        Input_Item item("exx_opt_orb_ecut");
        item.annotation = "the cut-off of plane wave expansion for opt ABFs";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "The cut-off of plane wave expansion, when the plane wave basis is used to optimize the radial ABFs. A reasonable choice is 60.";
        item.default_value = "0";
        item.unit = "Ry";
        item.availability = "calculation==gen_opt_abfs";
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
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "The threshold when solving for the zeros of spherical Bessel functions. A reasonable choice is 1e-12.";
        item.default_value = "1E-12";
        item.unit = "";
        item.availability = "calculation==gen_opt_abfs";
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
        Input_Item item("exx_real_number");
        item.annotation = "exx calculated in real or complex";
        item.category = "Exact Exchange (LCAO)";
        item.type = "String";
        item.description = R"(* True: Enforce LibRI to use double data type.
* False: Enforce LibRI to use complex data type. Setting it to True can effectively improve the speed of self-consistent calculations with hybrid functionals.)";
        item.default_value = "depends on the gamma_only option";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("exx_singularity_correction");
        item.annotation = "set the scheme of Coulomb singularity correction";
        item.category = "Exact Exchange (LCAO)";
        item.type = "String";
        item.description = R"(* spencer: see Phys. Rev. B 77, 193110 (2008).
* revised_spencer: see Phys. Rev. Mater. 5, 013807 (2021). Set the scheme of Coulomb singularity correction.)";
        item.default_value = "default";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.exx_singularity_correction);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.exx_singularity_correction == "default")
            {  
                std::string& dft_functional = para.input.dft_functional;
                std::string dft_functional_lower = dft_functional;
                std::transform(dft_functional.begin(), dft_functional.end(), dft_functional_lower.begin(), tolower);
                if (dft_functional_lower == "hf"
                    || dft_functional_lower == "pbe0" || dft_functional_lower == "b3lyp"
                    || dft_functional_lower == "scan0"
                    || dft_functional_lower == "muller" || dft_functional_lower == "power"
                    || dft_functional_lower == "wp22" 
                    || dft_functional_lower == "lc_pbe"
                    || dft_functional_lower == "lc_wpbe" 
                    || dft_functional_lower == "lrc_wpbe"
                    || dft_functional_lower == "lrc_wpbeh"
                    || dft_functional_lower == "cam_pbeh")
                {
                    para.input.exx_singularity_correction = "spencer";
                }
                else if (dft_functional_lower == "hse" || dft_functional_lower == "cwp22")
                {
                    para.input.exx_singularity_correction = "limits";
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rpa_ccp_rmesh_times");
        item.annotation = "how many times larger the radial mesh required for "
                          "calculating Columb potential is to that "
                          "of atomic orbitals";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Real";
        item.description = "How many times larger the radial mesh required is to that of atomic orbitals in the postprocess calculation of the bare Coulomb matrix for RPA, GW, etc.";
        item.default_value = "10";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("exx_symmetry_realspace");
        item.annotation = "whether to reduce real-space sector in Hexx calculation";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Boolean";
        item.description = R"(* False: only rotate k-space density matrix D(k) from irreducible k-points to accelerate diagonalization
* True: rotate both D(k) and Hexx(R) to accelerate both diagonalization and EXX calculation)";
        item.default_value = "True";
        item.unit = "";
        item.availability = "symmetry==1 and exx calculation (dft_fuctional==hse/hf/pbe0/scan0 or rpa==True)";
        read_sync_bool(input.exx_symmetry_realspace);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.symmetry != "1") { para.input.exx_symmetry_realspace = false; }
            };
        this->add_item(item);
    }
    {
        Input_Item item("out_ri_cv");
        item.annotation = "Whether to output the coefficient tensor C and ABFs-representation Coulomb matrix V";
        item.category = "Exact Exchange (LCAO)";
        item.type = "Boolean";
        item.description = "Whether to output the coefficient tensor C(R) and ABFs-representation Coulomb matrix V(R) for each atom pair and cell in real space.";
        item.default_value = "false";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_ri_cv);
        this->add_item(item);
    }
}
void ReadInput::item_dftu()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    // dft+u
    {
        Input_Item item("dft_plus_u");
        item.annotation = "DFT+U correction method";
        item.category = "DFT+U correction";
        item.type = "Integer";
        item.description = R"(Determines whether to calculate the plus U correction, which is especially important for correlated electrons.
* 1: Calculate plus U correction with radius-adjustable localized projections (with parameter onsite_radius).
* 2: Calculate plus U correction using first zeta of NAOs as projections (this is old method for testing).
* 0: Do not calculate plus U correction.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
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
                if (input.basis_type == "pw" && input.nspin != 4)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "WRONG ARGUMENTS, only nspin2 with PW base is not supported now");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("dft_plus_dmft");
        item.annotation = "true:DFT+DMFT; false: standard DFT calcullation(default)";
        item.category = "DFT+U correction";
        item.type = "Boolean";
        item.description = "Whether to enable DFT+DMFT calculation. True: DFT+DMFT; False: standard DFT calculation.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "basis_type==lcao";
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
        Input_Item item("orbital_corr");
        item.annotation = "which correlated orbitals need corrected ; d:2 "
                          ",f:3, do not need correction:-1";
        item.category = "DFT+U correction";
        item.type = "Vector of Integer (n values where n is the number of atomic types)";
        item.description = R"(Specifies which orbits need plus U correction for each atom type ( for atom type 1, 2, 3, respectively).
* -1: The plus U correction will not be calculated for this atom.
* 1: For p-electron orbits, the plus U correction is needed.
* 2: For d-electron orbits, the plus U correction is needed.
* 3: For f-electron orbits, the plus U correction is needed.)";
        item.default_value = "-1";
        item.unit = "";
        item.availability = "";
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
    {
        Input_Item item("hubbard_u");
        item.annotation = "Hubbard Coulomb interaction parameter U(ev)";
        item.category = "DFT+U correction";
        item.type = "Vector of Real (n values where n is the number of atomic types)";
        item.description = R"(Specifies the Hubbard Coulomb interaction parameter U (eV) in plus U correction, which should be specified for each atom unless the Yukawa potential is used.

[NOTE] Note: Since only the simplified scheme by Duradev is implemented, the 'U' here is actually U-effective, which is given by Hubbard U minus Hund J.)";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("yukawa_potential");
        item.annotation = "default: false";
        item.category = "DFT+U correction";
        item.type = "Boolean";
        item.description = R"(Determines whether to use the local screen Coulomb potential method to calculate the values of U and J.
* True: hubbard_u does not need to be specified.
* False: hubbard_u does need to be specified.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.yukawa_potential);
        this->add_item(item);
    }
    {
        Input_Item item("yukawa_lambda");
        item.annotation = "default:0.0";
        item.category = "DFT+U correction";
        item.type = "Real";
        item.description = "The screen length of Yukawa potential. If left to default, the screen length will be calculated as an average of the entire system. It's better to stick to the default setting unless there is a very good reason.";
        item.default_value = "Calculated on the fly.";
        item.unit = "";
        item.availability = "DFT+U with yukawa_potential = True.";
        read_sync_double(input.yukawa_lambda);
        this->add_item(item);
    }
    {
        Input_Item item("uramping");
        item.annotation = "increasing U values during SCF";
        item.category = "DFT+U correction";
        item.type = "Real";
        item.description = "Once uramping > 0.15 eV. DFT+U calculations will start SCF with U = 0 eV, namely normal LDA/PBE calculations. Once SCF restarts when drho<mixing_restart, U value will increase by uramping eV. SCF will repeat above calcuations until U values reach target defined in hubbard_u. As for uramping=1.0 eV, the recommendations of mixing_restart is around 5e-4.";
        item.default_value = "-1.0.";
        item.unit = "eV";
        item.availability = "DFT+U calculations with mixing_restart > 0.";
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
        add_double_bcast(sys.uramping);
        this->add_item(item);
    }
    {
        Input_Item item("omc");
        item.annotation = "the mode of occupation matrix control";
        item.category = "DFT+U correction";
        item.type = "Integer";
        item.description = R"(The parameter controls the form of occupation matrix control used.
* 0: No occupation matrix control is performed, and the onsite density matrix will be calculated from wavefunctions in each SCF step.
* 1: The first SCF step will use an initial density matrix read from a file named initial_onsite.dm, but for later steps, the onsite density matrix will be updated.
* 2: The same onsite density matrix from initial_onsite.dm will be used throughout the entire calculation.

[NOTE] The easiest way to create initial_onsite.dm is to run a DFT+U calculation, look for a file named onsite.dm in the OUT.prefix directory, and make replacements there. The format of the file is rather straight-forward.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.omc);
        this->add_item(item);
    }
    {
        Input_Item item("onsite_radius");
        item.annotation = "radius of the sphere for onsite projection (Bohr)";
        item.category = "DFT+U correction";
        item.type = "Real";
        item.description = R"(* The onsite_radius parameter facilitates modulation of the single-zeta portion of numerical atomic orbitals used for DFT+U projections.
* The modulation algorithm applies a smooth truncation to the orbital tail followed by normalization. A representative profile is $f(r)=\frac{1}{2}\left[1+\operatorname{erf}\!\left(\frac{r_c-r}{\sigma}\right)\right]$, where $r_c$ is the cutoff radius and $\sigma=\gamma r_c$ controls smoothness.)";
        item.default_value = "3.0";
        item.unit = "Bohr";
        item.availability = "dft_plus_u is set to 1";
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
}
} // namespace ModuleIO
