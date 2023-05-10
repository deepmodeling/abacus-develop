#ifndef EXX_ABFS_PCA_H
#define EXX_ABFS_PCA_H

#include "exx_abfs.h"
#include "../module_basis/module_ao/ORB_atomic_lm.h"
#include "../module_base/matrix.h"
#include "../module_base/element_basis_index.h"
#include <vector>
#include <RI/global/Tensor.h>

//	training data: lcaos[i] * lcaos[j]
//	old basis:     abfs
//	new basis:     to be constructed
//	( all lcaos and abfs on same atom )

class Exx_Abfs::PCA
{
public:
	static std::vector<std::vector<std::pair<std::vector<double>,RI::Tensor<double>>>> cal_PCA( 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &lcaos, 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs,		// abfs must be orthonormal
		const double kmesh_times );
private:
	static RI::Tensor<double> get_sub_matrix( 
		const RI::Tensor<double> & m,
		const size_t & T,
		const size_t & L,
		const ModuleBase::Element_Basis_Index::Range & range,
		const ModuleBase::Element_Basis_Index::IndexLNM & index );
	static RI::Tensor<double> get_column_mean0_matrix( const RI::Tensor<double> & m );
};

#endif	// EXX_ABFS_PCA_H