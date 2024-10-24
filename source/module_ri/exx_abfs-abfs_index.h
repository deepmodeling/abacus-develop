#ifndef EXX_ABFS_ABFS_INDEX_H
#define EXX_ABFS_ABFS_INDEX_H

#include "exx_abfs.h"

#include <vector>
#include "../module_base/element_basis_index.h"
#include "../module_basis/module_ao/ORB_atomic_lm.h"

	class LCAO_Orbitals;

namespace Exx_Abfs
{
namespace Abfs_Index
{
	extern ModuleBase::Element_Basis_Index::Range construct_range( const LCAO_Orbitals &orb );
	extern ModuleBase::Element_Basis_Index::Range construct_range( const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb );
}
}

#endif	// EXX_ABFS_ABFS_INDEX_H