#ifndef EXX_INFO_H
#define EXX_INFO_H

#include "exx_info_global.h"
#include "exx_info_lip.h"
#include "exx_info_ri.h"
#include "exx_info_opt_abfs.h"

struct Exx_Info
{
    Exx_Info_Global info_global;
    Exx_Info_Lip info_lip;

    struct Exx_Info_RI
    {
        const std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>> &coulomb_param;

        bool real_number = false;
        bool coul_moment = false;
        bool rotate_abfs = false;

        double pca_threshold = 0;
        std::vector<std::string> files_abfs;
        std::vector<std::string> files_shrink_abfs;
        double C_threshold = 0;
        double V_threshold = 0;
        double V_threshold_long = 0;
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
    Exx_Info_Opt_ABFs info_opt_abfs;

    void sync_from_global()
    {
        info_lip.ccp_type = info_global.ccp_type;
        info_lip.hse_omega = info_global.hse_omega;
        info_ri.coulomb_param = info_global.coulomb_param;
    }
};

/// Forward declaration for Input_para (full definition in input_parameter.h)
struct Input_para;

/// Initialize an Exx_Info object from input parameters.
/// Returns true if opt_orb mode is requested (generate_opt_orb).
bool init_exx_info(Exx_Info& exx_info, const Input_para& inp);

#endif