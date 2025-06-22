#ifndef READ_WF2RHO_PW_H
#define READ_WF2RHO_PW_H

#include "source_basis/module_pw/pw_basis_k.h"
#include "module_elecstate/module_charge/charge.h"

#include <string>
#include <fstream>

namespace ModuleIO
{
/**
 * @brief read wave functions and occupation numbers to charge density
 *
 * @param pw_wfc pw basis for wave functions
 * @param symm symmetry
 * @param nkstot total number of k points
 * @param isk k index to spin index
 * @param chg charge density
 */

void read_wf2rho_pw(
		const ModulePW::PW_Basis_K* pw_wfc,
		ModuleSymmetry::Symmetry& symm,
		Charge& chg,
        const std::string &readin_dir,
		const int kpar,
		const int my_pool,
		const int my_rank,
        const int nproc_in_pool,
        const int rank_in_pool,
		const int nbands,
		const int nspin,
		const int npol,
		const int nkstot,
		const std::vector<int> &ik2iktot,
		const std::vector<int> &isk,
		std::ofstream &ofs_running);

} // namespace ModuleIO

#endif
