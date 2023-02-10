#include "src_ri/abfs.h"
#include <iostream>

static void test_adjs( std::ostream & os )
{
	for( size_t it=0; it!=GlobalC::ucell.ntype; ++it )
		for( size_t ia=0; ia!=GlobalC::ucell.atoms[it].na; ++ia )
		{
			const int iat = GlobalC::ucell.itia2iat(it,ia);
			const std::map<size_t,std::vector<Abfs::Vector3_Order<int>>> adjs = Abfs::get_adjs(iat);
			os<<"@@@\t"<<iat<<std::endl;
			for( const auto & atom2 : adjs )
			{
				const size_t iat2 = atom2.first;
				for( const Abfs::Vector3_Order<int> &box2 : atom2.second )
					os<<iat2<<"\t"<<box2<<std::endl;
			}
		}
}
