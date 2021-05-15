#ifndef EXX_ABFS_PCA_H
#define EXX_ABFS_PCA_H

#include "exx_abfs.h"
#include "module_ORB/ORB_atomic_lm.h"
#include "src_global/matrix.h"
#include "src_global/element_basis_index.h"
#include <vector>

//	training data: lcaos[i] * lcaos[j]
//	old basis:     abfs
//	new basis:     to be constructed
//	( all lcaos and abfs on same atom )

class Exx_Abfs::PCA
{
public:
	static vector<vector<pair<vector<double>,matrix>>> cal_PCA( 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &lcaos, 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &abfs,		// abfs must be orthonormal
		const double kmesh_times );
private:
	static matrix get_sub_matrix( 
		const matrix & m,
		const size_t & T,
		const size_t & L,
		const Element_Basis_Index::Range & range,
		const Element_Basis_Index::IndexLNM & index );
	static matrix get_column_mean0_matrix( const matrix & m );
};

#endif	// EXX_ABFS_PCA_H