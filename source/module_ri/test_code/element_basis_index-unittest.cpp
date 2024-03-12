//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-06-02
//==========================================================

#include "module_base/element_basis_index.h"
#include "element_basis_index-test.h"

ModuleBase::Element_Basis_Index::Range construct_range()
{
	ModuleBase::Element_Basis_Index::Range range;
	range.resize(2);
	
	range[0].resize(3);
	range[0][0].N = range[0][1].N = 2;
	range[0][2].N = 1;
	
	range[1].resize(2);
	range[1][0].N = 3;
	range[1][1].N = 2;
	
	for( size_t T=0; T!=range.size(); ++T )
		for( size_t L=0; L!=range[T].size(); ++L )
			range[T][L].M=2*L+1;
	
	return range;
}

int main()
{
	const ModuleBase::Element_Basis_Index::Range range = construct_range();
	std::cout<<"range:"<<std::endl<<range<<std::endl;
	
	const ModuleBase::Element_Basis_Index::IndexLNM index = ModuleBase::Element_Basis_Index::construct_index( range );
	std::cout<<"index:"<<std::endl<<index<<std::endl;

	return 0;
}

/*
range:
0
	0	2	1
	1	2	3
	2	1	5
1
	0	3	1
	1	2	3

index:
0:	13
	0
		0
			0	0
		1
			0	1
	1
		0
			0	2
			1	3
			2	4
		1
			0	5
			1	6
			2	7
	2
		0
			0	8
			1	9
			2	10
			3	11
			4	12
1:	9
	0
		0
			0	0
		1
			0	1
		2
			0	2
	1
		0
			0	3
			1	4
			2	5
		1
			0	6
			1	7
			2	8
*/