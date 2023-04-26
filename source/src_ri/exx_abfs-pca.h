#ifndef EXX_ABFS_PCA_H
#define EXX_ABFS_PCA_H

#include "exx_abfs.h"
#include "../module_basis/module_ao/ORB_atomic_lm.h"
#include "../module_base/matrix.h"
#include "../module_base/element_basis_index.h"
#include <vector>

//	training data: lcaos[i] * lcaos[j]
//	old basis:     abfs
//	new basis:     to be constructed
//	( all lcaos and abfs on same atom )

class Exx_Abfs::PCA
{
public:
	static std::vector<std::vector<std::pair<std::vector<double>,ModuleBase::matrix>>> cal_PCA( 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &lcaos, 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs,		// abfs must be orthonormal
		const double kmesh_times );
private:
	static ModuleBase::matrix get_sub_matrix( 
		const ModuleBase::matrix & m,
		const size_t & T,
		const size_t & L,
		const ModuleBase::Element_Basis_Index::Range & range,
		const ModuleBase::Element_Basis_Index::IndexLNM & index );
	static ModuleBase::matrix get_column_mean0_matrix( const ModuleBase::matrix & m );
};

#endif	// EXX_ABFS_PCA_H