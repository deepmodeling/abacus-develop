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
 * @param out_app_flag: whether to append to existing file.
 * @param gamma_only: gamma_only algorithm or not.
 * @param istep: index of the ion step starting from 0. If < 0, the step number is not included in the file name.
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
			const bool out_app_flag,
			const bool gamma_only,
			const int istep=-1);

}
#endif
