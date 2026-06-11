#ifndef EXX_INFO_OPT_ABFS_H
#define EXX_INFO_OPT_ABFS_H

#include <vector>
#include <string>

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

#endif