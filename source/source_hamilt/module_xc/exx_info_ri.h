#ifndef EXX_INFO_RI_H
#define EXX_INFO_RI_H

#include "source_lcao/module_ri/conv_coulomb_pot_k.h"

#include <vector>
#include <string>
#include <map>

struct Exx_Info_RI
{
    std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string, std::string>>> coulomb_param;

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

    int abfs_Lmax = 0;
};

#endif