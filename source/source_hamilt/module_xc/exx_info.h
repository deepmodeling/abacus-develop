#ifndef EXX_INFO_H
#define EXX_INFO_H

#include "source_lcao/module_ri/conv_coulomb_pot_k.h"

#include <vector>
#include <map>
#include <string>

struct Exx_Info
{
    struct Exx_Info_Global
    {
        bool cal_exx = false;

        std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string, std::string>>> coulomb_param;

        // Fock:
        //      "alpha":        "0"
        //      "singularity_correction":   "limits" / "spencer" / "revised_spencer" / "massidda" / "carrier"
        //      "lambda":       "0.3"
        //      "Rcut"
        // Erfc:
        //      "alpha":        "0"
        //      "omega":        "0.11"
        //      "singularity_correction":   "limits" / "spencer" / "revised_spencer"
        //      "Rcut"

        Conv_Coulomb_Pot_K::Ccp_Type ccp_type;
        double hybrid_alpha = 0.25;
        double hse_omega = 0.11;
        double mixing_beta_for_loop1 = 1.0;

        bool separate_loop = true;
        size_t hybrid_step = 1;
    };
    Exx_Info_Global info_global;

    // WARNING: Lifetime dependency
    // Exx_Info_Lip and Exx_Info_RI hold const references to members of Exx_Info_Global.
    // Their lifetimes MUST NOT exceed the Exx_Info_Global they were initialized from.
    // Currently this is safe because Exx_Info owns both info_global and info_lip/info_ri,
    // and GlobalC::exx_info is a global variable with program-lifetime scope.
    // However, if Exx_Info is ever copied or moved, these references will become dangling.
    // Do NOT add copy/move constructors to Exx_Info without addressing this.
    struct Exx_Info_Lip
    {
        const Conv_Coulomb_Pot_K::Ccp_Type& ccp_type;  // reference to info_global.ccp_type
        const double& hse_omega;                        // reference to info_global.hse_omega
        double lambda = 0.3;

        Exx_Info_Lip(const Exx_Info::Exx_Info_Global& info_global)
            : ccp_type(info_global.ccp_type),
              hse_omega(info_global.hse_omega) {}
    };
    Exx_Info_Lip info_lip;

    struct Exx_Info_RI
    {
        // reference to info_global.coulomb_param
        // Note: coulomb_param is populated AFTER Exx_Info construction (in input_conv.cpp),
        // so this reference is valid but initially points to an empty map.
        // The data is filled later into info_global, and this reference sees it via aliasing.
        const std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string, std::string>>> &coulomb_param;

        bool real_number = false;
        bool coul_moment = false;
        bool rotate_abfs = false;

        double pca_threshold = 0;
        std::vector<std::string> files_abfs;
        std::vector<std::string> files_shrink_abfs;
        double C_threshold = 0;
        double V_threshold = 0;
        double dm_threshold = 0;
        double C_grad_threshold = 0;
        double V_grad_threshold = 0;
        double C_grad_R_threshold = 0;
        double V_grad_R_threshold = 0;
        double ccp_rmesh_times = 10;
        bool exx_symmetry_realspace = true;
        double kmesh_times = 4;
        double Cs_inv_thr = -1;

        double shrink_abfs_pca_thr = -1;
        double shrink_LU_inv_thr = 1e-6;
        double multip_moments_threshold = 1e-10;
        double exx_cs_inv_thr = -1;

        int abfs_Lmax = 0; // tmp

        Exx_Info_RI(const Exx_Info::Exx_Info_Global& info_global)
            : coulomb_param(info_global.coulomb_param)
        {
        }
    };
    Exx_Info_RI info_ri;

    struct Exx_Info_Opt_ABFs
    {
        int abfs_Lmax = 0;
        double ecut_exx = 60;
        double tolerence = 1E-12;
        std::vector<std::string> files_jles;

        double pca_threshold = 0;
        std::vector<std::string> files_abfs;

        double kmesh_times = 4;
    };
    Exx_Info_Opt_ABFs info_opt_abfs;

    Exx_Info() : info_lip(this->info_global), info_ri(this->info_global)
    {
    }
};

namespace GlobalC
{
    extern Exx_Info exx_info;
}

#endif
