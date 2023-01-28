#ifndef DENSITY_MATRIX
#define DENSITY_MATRIX

#include <string>
#include "module_orbital/parallel_orbitals.h"
#include "module_base/abfs-vector3_order.h"

namespace ModuleIO
{
	void read_dm(const int &is, const std::string &fn, double*** DM, double** DM_R);
	void write_dm(const int &is, const int &iter, const std::string &fn, const int &precision, const int &out_dm, double*** DM);
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

