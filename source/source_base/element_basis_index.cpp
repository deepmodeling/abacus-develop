//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-06-02
//==========================================================

#include "element_basis_index.h"
#include "../../source_basis/module_ao/ORB_read.h"
#include "../../source_basis/module_ao/ORB_atomic_lm.h"

namespace ModuleBase
{

Element_Basis_Index::IndexLNM
Element_Basis_Index::construct_index( const Range &range )
{
	IndexLNM index;
	index.resize( range.size() );
	for( size_t T=0; T!=range.size(); ++T )
	{
		size_t count=0;
		index[T].resize( range[T].size() );
		for( size_t L=0; L!=range[T].size(); ++L )
		{
			index[T][L].resize( range[T][L].N );
			for( size_t N=0; N!=range[T][L].N; ++N )
			{
				index[T][L][N].resize( range[T][L].M );
				for( size_t M=0; M!=range[T][L].M; ++M )
				{
					index[T][L][N][M] = count;
					++count;
				}
			}
			index[T][L].N = range[T][L].N;
			index[T][L].M = range[T][L].M;
		}
		index[T].count_size = count;
	}
	return index;
}


ModuleBase::Element_Basis_Index::Range
Element_Basis_Index::construct_range( const LCAO_Orbitals &orb )
{
	ModuleBase::Element_Basis_Index::Range range;
	range.resize( orb.get_ntype() );
	for( size_t T=0; T!=range.size(); ++T )
	{
		range[T].resize( orb.Phi[T].getLmax()+1 );
		for( size_t L=0; L!=range[T].size(); ++L )
		{
			range[T][L].N = orb.Phi[T].getNchi(L);
			range[T][L].M = 2*L+1;
		}
	}
	return range;
}


ModuleBase::Element_Basis_Index::Range
Element_Basis_Index::construct_range( const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb )
{
	ModuleBase::Element_Basis_Index::Range range;
	range.resize( orb.size() );
	for( size_t T=0; T!=range.size(); ++T )
	{
		range[T].resize( orb[T].size() );
		for( size_t L=0; L!=range[T].size(); ++L )
		{
			range[T][L].N = orb[T][L].size();	
			range[T][L].M = 2*L+1;	
		}			
	}
	return range;
}

}