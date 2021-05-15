#include "exx_abfs-abfs_index.h"
#include "module_ORB/ORB_read.h"

Element_Basis_Index::Range
	Exx_Abfs::Abfs_Index::construct_range( const LCAO_Orbitals &orb )
{
	Element_Basis_Index::Range range;
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


Element_Basis_Index::Range
	Exx_Abfs::Abfs_Index::construct_range( const vector<vector<vector<Numerical_Orbital_Lm>>> &orb )
{
	Element_Basis_Index::Range range;
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
