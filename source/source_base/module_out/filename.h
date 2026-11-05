#ifndef FILENAME_H
#define FILENAME_H
#include <vector>
#include <string>

namespace ModuleIO
{

/**
 * Generates the filename for the output files
 * @param directory: directory of the file 
 * @param property: wave function (wf), charge density (chg) or matrix (mat)
 * @param basis: nao or pw 
 * @param ik_local: index of the k-points within each pool, and starting from 0.
 * @param ik2iktot: map from ik to iktot
 * @param nspin: number of spin channels, 1,2 or 4 
 * @param nkstot: number of total k-points 
 * @param out_type: two types of output file format, 1 for .txt and 2 for .dat (binary)
 * @param out_app_flag: 0=fixed filename, overwritten at each write; 1=fixed filename, appended;
 *                      2=subdirectory layout with a g{step} suffix.
 * @param gamma_only: gamma_only algorithm or not.
 * @param istep: index of the ion step starting from 0. The step number is only embedded in the
 *               filename in mode 2 (out_app_flag == 2), and only when istep >= 0.
 * @param iter: index of the electronic step starting from 0. The step number is only embedded in
 *              the filename in mode 2 (out_app_flag == 2), and only when iter >= 0.
 * @return The generated filename.
 */
std::string filename_output(
            const std::string &directory,
            const std::string &property,
            const std::string &basis,
            const int ik_local,
            const std::vector<int> &ik2iktot,
            const int nspin,
            const int nkstot,
            const int out_type,
            const int out_app_flag,
            const bool gamma_only,
            const int istep=-1,
            const int iter=-1);

}
#endif
