#ifndef EXX_ABFS_SCREEN_TEST_H
#define EXX_ABFS_SCREEN_TEST_H

#include <map>
#include <string>
#include "src_ri/abfs.h"

static void test_screen( const std::string & file_name, const std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,double>>> & m )
{
	std::ofstream ofs(file_name);
	for( const auto m1 : m )
		for( const auto m2 : m1.second )
			for( const auto m3 : m2.second )
				ofs<<m1.first<<"\t"<<m2.first<<"\t"<<m3.first<<"\t"<<m3.second<<std::endl;
	ofs.close();
}

#endif
