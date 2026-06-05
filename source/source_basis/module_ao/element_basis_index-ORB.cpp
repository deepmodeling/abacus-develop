#include "element_basis_index-ORB.h"

#include "ORB_read.h"
#include "ORB_atomic_lm.h"

namespace ModuleBase
{

ModuleBase::Element_Basis_Index::Range
Element_Basis_Index::construct_range( const LCAO_Orbitals &orb )
{
	ModuleBase::Element_Basis_Index::Range range;
	range.resize( orb.get_ntype() );
	for( std::size_t T=0; T!=range.size(); ++T )
	{
		range[T].resize( orb.Phi[T].getLmax()+1 );
		for( std::size_t L=0; L!=range[T].size(); ++L )
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
	for( std::size_t T=0; T!=range.size(); ++T )
	{
		range[T].resize( orb[T].size() );
		for( std::size_t L=0; L!=range[T].size(); ++L )
		{
			range[T][L].N = orb[T][L].size();
			range[T][L].M = 2*L+1;
		}
	}
	return range;
}

}