#ifndef EXX_ABFS_IO_H
#define EXX_ABFS_IO_H

#include "exx_abfs.h"

#include <map>
#include <vector>
#include "../module_basis/module_ao/ORB_atomic_lm.h"
#include "../module_base/matrix.h"
#include "../module_base/element_basis_index.h"
#include "module_cell/klist.h"
#ifdef __MPI
#include "mpi.h"
#endif

	class LCAO_Orbitals;

namespace Exx_Abfs
{
namespace IO
{
	extern std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> construct_abfs(
		const LCAO_Orbitals &orbs,
		const std::vector<std::string> &files_abfs,
		const double kmesh_times=1 );				// close dK, keep Kcut

	extern std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> construct_abfs(
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs_pre,
		const LCAO_Orbitals &orbs,
		const std::vector<std::string> &files_abfs,
		const double kmesh_times=1 );				// close dK, keep Kcut

  // private:
	extern std::vector<std::vector<Numerical_Orbital_Lm>> construct_abfs_T(
		const std::string & file_name,
		const int &T,
		const int &nk,
		const double &dk,
		const double &dr_uniform);
}
}

#endif	// EXX_ABFS_IO_H
