//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-06-02
//==========================================================

#include "element_basis_index.h"

#include <stdexcept>
#include <vector>

namespace ModuleBase
{

Element_Basis_Index::IndexLNM
Element_Basis_Index::construct_index( const Range &range )
{
	IndexLNM index;
	index.resize( range.size() );
	for( std::size_t T=0; T!=range.size(); ++T )
	{
		std::size_t count=0;
		index[T].resize( range[T].size() );
		for( std::size_t L=0; L!=range[T].size(); ++L )
		{
			index[T][L].resize( range[T][L].N );
			for( std::size_t N=0; N!=range[T][L].N; ++N )
			{
				index[T][L][N].resize( range[T][L].M );
				for( std::size_t M=0; M!=range[T][L].M; ++M )
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

Element_Basis_Index::IndexLNM
Element_Basis_Index::construct_index( const Range &range, const IndexPermutation &old_to_new )
{
	if( old_to_new.empty() )
	{
		return construct_index(range);
	}

	IndexLNM index;
	index.resize( range.size() );
	for( std::size_t T=0; T!=range.size(); ++T )
	{
		std::size_t count_size = 0;
		for( std::size_t L=0; L!=range[T].size(); ++L )
		{
			count_size += range[T][L].N * range[T][L].M;
		}

		const bool use_permutation = (T < old_to_new.size() && !old_to_new[T].empty());
		if( use_permutation )
		{
			if( old_to_new[T].size() != count_size )
			{
				throw std::invalid_argument("Element_Basis_Index::construct_index permutation size mismatch.");
			}

			std::vector<bool> seen(count_size, false);
			for( const std::size_t new_index : old_to_new[T] )
			{
				if( new_index >= count_size || seen[new_index] )
				{
					throw std::invalid_argument("Element_Basis_Index::construct_index permutation is not bijective.");
				}
				seen[new_index] = true;
			}
		}

		std::size_t count = 0;
		index[T].resize( range[T].size() );
		for( std::size_t L=0; L!=range[T].size(); ++L )
		{
			index[T][L].resize( range[T][L].N );
			for( std::size_t N=0; N!=range[T][L].N; ++N )
			{
				index[T][L][N].resize( range[T][L].M );
				for( std::size_t M=0; M!=range[T][L].M; ++M )
				{
					index[T][L][N][M] = use_permutation ? old_to_new[T][count] : count;
					++count;
				}
			}
			index[T][L].N = range[T][L].N;
			index[T][L].M = range[T][L].M;
		}
		index[T].count_size = count_size;
	}
	return index;
}

}
