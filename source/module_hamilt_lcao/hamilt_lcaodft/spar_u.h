#ifndef SPARSE_FORMAT_U_H 
#define SPARSE_FORMAT_U_H

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"


namespace sparse_format
{

	void cal_HR_dftu(
		const Parallel_Orbitals &pv,
		std::set<Abfs::Vector3_Order<int>> &all_R_coor,
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &SR_sparse,
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> *HR_sparse,
		const int &current_spin, 
		const double &sparse_thr);

	void cal_HR_dftu_soc(
		const Parallel_Orbitals &pv,
		std::set<Abfs::Vector3_Order<int>> &all_R_coor,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> &SR_soc_sparse,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> &HR_soc_sparse,
		const int &current_spin, 
		const double &sparse_thr);

}

#endif
