#ifndef SPARSE_FORMAT_U_H 
#define SPARSE_FORMAT_U_H

#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_pw/module_pwdft/global.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_dftu/dftu.h" // mohan add 20251107


namespace sparse_format
{

	void cal_HR_dftu(
        Plus_U &dftu, // mohan add 2025-11-07
		const Parallel_Orbitals &pv,
		std::set<Abfs::Vector3_Order<int>> &all_R_coor,
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &SR_sparse,
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> *HR_sparse,
		const int &current_spin, 
		const double &sparse_thr);

	void cal_HR_dftu_soc(
        Plus_U &dftu, // mohan add 2025-11-07
		const Parallel_Orbitals &pv,
		std::set<Abfs::Vector3_Order<int>> &all_R_coor,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> &SR_soc_sparse,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> &HR_soc_sparse,
		const int &current_spin, 
		const double &sparse_thr);

}

#endif
