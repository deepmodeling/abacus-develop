#ifndef EXX_ABFS_ABFS_INDEX_H
#define EXX_ABFS_ABFS_INDEX_H

#include "exx_abfs.h"

#include<vector>
#include "src_global/element_basis_index.h"
#include "module_ORB/ORB_atomic_lm.h"

class LCAO_Orbitals;

class Exx_Abfs::Abfs_Index
{
public:
	static size_t get_index_index( 
		const Element_Basis_Index::IndexLNM &index_A, const size_t &TA, const size_t &LA, const size_t &NA, const size_t &MA, 
		const Element_Basis_Index::IndexLNM &index_B, const size_t &TB, const size_t &LB, const size_t &NB, const size_t &MB )
	{	return index_A[TA][LA][NA][MA] * index_B[TB].count_size + index_B[TB][LB][NB][MB];	}
	static Element_Basis_Index::Range construct_range( const LCAO_Orbitals &orb );	
	static Element_Basis_Index::Range construct_range( const vector<vector<vector<Numerical_Orbital_Lm>>> &orb );
};

#endif	// EXX_ABFS_ABFS_INDEX_H