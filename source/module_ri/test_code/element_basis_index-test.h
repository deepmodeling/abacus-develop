//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-06-02
//==========================================================

#ifndef ELEMENT_BASIS_INDEX_TEST_H
#define ELEMENT_BASIS_INDEX_TEST_H

#include <iostream>

static std::ostream & operator << (std::ostream & os, const ModuleBase::Element_Basis_Index::Range & range)
{
	for( size_t T=0; T!=range.size(); ++T )
	{	
		os<<T<<std::endl;
		for( size_t L=0; L!=range[T].size(); ++L )
		{
			os<<"\t"<<L<<"\t"<<range[T][L].N<<"\t"<<range[T][L].M<<std::endl;
		}
	}
	return os;
}

static std::ostream & operator << (std::ostream & os, const ModuleBase::Element_Basis_Index::IndexLNM & index)
{
	for( size_t T=0; T!=index.size(); ++T )
	{
		os<<T<<":\t"<<index[T].count_size<<std::endl;;
		for( size_t L=0; L!=index[T].size(); ++L )
		{
			os<<"\t"<<L<<std::endl;
			for( size_t N=0; N!=index[T][L].size(); ++N )
			{
				os<<"\t\t"<<N<<std::endl;
				for( size_t M=0; M!=index[T][L][N].size(); ++M )
				{
					os<<"\t\t\t"<<M<<"\t"<<index[T][L][N][M]<<std::endl;
				}
			}
		}
	}
	return os;
}

#endif