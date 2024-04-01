#ifndef WRITE_DM_SPARSE_H
#define WRITE_DM_SPARSE_H

#include <string>
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_base/abfs-vector3_order.h"
#include "single_R_io.h"

namespace ModuleIO
{
	void write_dm1(const int &is, const int &istep, double** dm2d, const Parallel_Orbitals* ParaV,
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &DMR_sparse);

	void get_dm_sparse(const int &is, double** dm2d, const Parallel_Orbitals* ParaV,
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &DMR_sparse);

	void write_dm_sparse(const int &is, const int &istep, const Parallel_Orbitals* ParaV,
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &DMR_sparse);

	void destroy_dm_sparse(
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &DMR_sparse);
}

#endif

