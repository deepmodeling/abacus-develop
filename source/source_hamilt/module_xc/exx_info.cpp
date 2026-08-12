#include "exx_info.h"

#include "source_io/module_parameter/input_parameter.h"
#include "source_base/global_function.h"

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

//----------------------------------------------------------
// Initialize Exx_Info from input parameters.
// Extracted from Input_Conv::Convert() to allow ESolver to own
// its own Exx_Info object without depending on the global instance.
// Peize Lin add 2018-06-20, refactored 2026.
//----------------------------------------------------------
bool init_exx_info(Exx_Info& exx_info, const Input_para& inp)
{
    std::string dft_functional_lower = inp.dft_functional;
    std::transform(inp.dft_functional.begin(),
                   inp.dft_functional.end(),
                   dft_functional_lower.begin(),
                   tolower);
    bool generate_opt_orb = false;
    if (dft_functional_lower == "hf"
    || dft_functional_lower == "pbe0" || dft_functional_lower == "b3lyp" || dft_functional_lower == "hse"
    || dft_functional_lower == "scan0"
    || dft_functional_lower == "muller" || dft_functional_lower == "power"
    || dft_functional_lower == "cwp22" || dft_functional_lower == "wp22"
    || dft_functional_lower == "lc_pbe"
    || dft_functional_lower == "lc_wpbe"
    || dft_functional_lower == "lrc_wpbe"
    || dft_functional_lower == "lrc_wpbeh"
    || dft_functional_lower == "cam_pbeh")
    {
        exx_info.info_global.cal_exx = true;

        exx_info.info_global.hybrid_alpha = 0;
        std::vector<double> fock_alpha(inp.exx_fock_alpha.size());
        for(std::size_t i=0; i<fock_alpha.size(); ++i)
        {
            fock_alpha[i] = std::stod(inp.exx_fock_alpha[i]);
            exx_info.info_global.hybrid_alpha = std::max(std::abs(fock_alpha[i]), exx_info.info_global.hybrid_alpha);
        }
        std::vector<double> erfc_alpha(inp.exx_erfc_alpha.size());
        for(std::size_t i=0; i<erfc_alpha.size(); ++i)
        {
            erfc_alpha[i] = std::stod(inp.exx_erfc_alpha[i]);
            exx_info.info_global.hybrid_alpha = std::max(std::abs(erfc_alpha[i]), exx_info.info_global.hybrid_alpha);
        }
        assert(exx_info.info_global.hybrid_alpha>0);
        for(std::size_t i=0; i<fock_alpha.size(); ++i)
            { fock_alpha[i] /= exx_info.info_global.hybrid_alpha; }
        for(std::size_t i=0; i<erfc_alpha.size(); ++i)
            { erfc_alpha[i] /= exx_info.info_global.hybrid_alpha; }

        if(!fock_alpha.empty())
        {
            if(inp.basis_type == "lcao")
            {
                exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock].resize(fock_alpha.size());
                for(std::size_t i=0; i<fock_alpha.size(); ++i)
                {
                    exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock] = {{
                        {"alpha", ModuleBase::GlobalFunc::TO_STRING(fock_alpha[i])},
                        {"singularity_correction", inp.exx_singularity_correction} }};
                }
            }
            else if(inp.basis_type == "lcao_in_pw")
            {
                assert(fock_alpha.size() == inp.exx_fock_lambda.size());
                exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock].resize(fock_alpha.size());
                for(std::size_t i=0; i<fock_alpha.size(); ++i)
                {
                    exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock] = {{
                        {"alpha", ModuleBase::GlobalFunc::TO_STRING(fock_alpha[i])},
                        {"lambda", inp.exx_fock_lambda[i]} }};
                }
            }
            else if(inp.basis_type == "pw")
            {
                exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock].resize(fock_alpha.size());
                for(std::size_t i=0; i<fock_alpha.size(); ++i)
                {
                    exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock] = {{
                        {"alpha", ModuleBase::GlobalFunc::TO_STRING(fock_alpha[i])} }};
                }
            }
            else
            {
                throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
            }
        }
        if(!erfc_alpha.empty())
        {
            assert(erfc_alpha.size() == inp.exx_erfc_omega.size());
            if(inp.basis_type == "lcao")
            {
                exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc].resize(erfc_alpha.size());
                for(std::size_t i=0; i<erfc_alpha.size(); ++i)
                {
                    exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc] = {{
                        {"alpha", ModuleBase::GlobalFunc::TO_STRING(erfc_alpha[i])},
                        {"omega", ModuleBase::GlobalFunc::TO_STRING(inp.exx_erfc_omega[i])},
                        {"singularity_correction", inp.exx_singularity_correction} }};
                }
            }
            else if(inp.basis_type == "pw" || inp.basis_type == "lcao_in_pw")
            {
                exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc].resize(erfc_alpha.size());
                for(std::size_t i=0; i<erfc_alpha.size(); ++i)
                {
                    exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc] = {{
                        {"alpha", ModuleBase::GlobalFunc::TO_STRING(erfc_alpha[i])},
                        {"omega", ModuleBase::GlobalFunc::TO_STRING(inp.exx_erfc_omega[i])} }};
                }
            }
        }
    }
#ifdef __EXX
    else if (dft_functional_lower == "opt_orb")
    {
        exx_info.info_global.cal_exx = false;
        generate_opt_orb = true;
    }
#endif
    else
    {
        exx_info.info_global.cal_exx = false;
    }

    if (inp.rpa && exx_info.info_global.coulomb_param.empty())
    {
        if (inp.basis_type != "lcao")
        {
            throw std::invalid_argument("RPA currently expects basis_type=lcao when initializing RI Coulomb parameters.");
        }
        exx_info.info_global.hybrid_alpha = 1.0;
        exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hf;
        exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock] = {{
            {"alpha", "1"},
            {"singularity_correction", inp.exx_singularity_correction}
        }};
    }

    // info_global.ccp_type will be removed in the future. these codes for pw and lcao_in_pw temporarily
    if (dft_functional_lower == "hf"
     || dft_functional_lower == "pbe0" || dft_functional_lower == "b3lyp"
     || dft_functional_lower == "scan0"
     || dft_functional_lower == "muller" || dft_functional_lower == "power")
    {
        exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hf;
    }
    // use the error function erf(w|r-r'|), exx just has the short-range part
    else if (dft_functional_lower == "hse"
          || dft_functional_lower == "cwp22")
    {
        exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Erfc;
    }
    // use the error function erf(w|r-r'|), exx just has the long-range part
    else if ( dft_functional_lower == "wp22" )
    {
        exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Erf;
    }

    if (exx_info.info_global.cal_exx
#ifdef __EXX
        || generate_opt_orb
        || inp.rpa
#endif
        )
    {
        // EXX case, convert all EXX related variables
        if(!inp.exx_erfc_omega.empty())
            { exx_info.info_global.hse_omega = std::stod(inp.exx_erfc_omega[0]); }
        if(!inp.exx_fock_lambda.empty())
            { exx_info.info_lip.lambda = std::stod(inp.exx_fock_lambda[0]); }
        exx_info.info_global.separate_loop = inp.exx_separate_loop;
        exx_info.info_global.hybrid_step = inp.exx_hybrid_step;
        exx_info.info_global.mixing_beta_for_loop1 = inp.exx_mixing_beta;

        exx_info.info_ri.real_number = std::stoi(inp.exx_real_number);
        exx_info.info_ri.pca_threshold = inp.exx_pca_threshold;
        exx_info.info_ri.C_threshold = inp.exx_c_threshold;
        exx_info.info_ri.V_threshold = inp.exx_v_threshold;
        exx_info.info_ri.dm_threshold = inp.exx_dm_threshold;
        exx_info.info_ri.C_grad_threshold = inp.exx_c_grad_threshold;
        exx_info.info_ri.V_grad_threshold = inp.exx_v_grad_threshold;
        exx_info.info_ri.C_grad_R_threshold = inp.exx_c_grad_r_threshold;
        exx_info.info_ri.V_grad_R_threshold = inp.exx_v_grad_r_threshold;
        exx_info.info_ri.ccp_rmesh_times = std::stod(inp.exx_ccp_rmesh_times);
        exx_info.info_ri.exx_symmetry_realspace = inp.exx_symmetry_realspace;
        exx_info.info_ri.Cs_inv_thr = inp.exx_cs_inv_thr;
        exx_info.info_ri.shrink_abfs_pca_thr = inp.shrink_abfs_pca_thr;
        exx_info.info_ri.shrink_LU_inv_thr = inp.shrink_LU_inv_thr;
        exx_info.info_ri.coul_moment = inp.exx_coul_moment;
        exx_info.info_ri.rotate_abfs = inp.exx_rotate_abfs;
        exx_info.info_ri.multip_moments_threshold = inp.exx_multip_moments_threshold;
        exx_info.info_opt_abfs.pca_threshold = inp.exx_pca_threshold;
        exx_info.info_opt_abfs.abfs_Lmax = inp.exx_opt_orb_lmax;
        exx_info.info_opt_abfs.ecut_exx = inp.exx_opt_orb_ecut;
        exx_info.info_opt_abfs.tolerence = inp.exx_opt_orb_tolerence;

        // Space-group symmetry is supported for LCAO EXX (nspin=1,2 via restore_dm/restore_HR;
        // nspin=4/SOC via restore_dm + restore_HR_nspin4), so symmetry=1 is honored here.

        exx_info.sync_from_global();
    }

    return generate_opt_orb;
}
