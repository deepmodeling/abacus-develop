#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <fstream>
#include <cstdio>
#include "../memory.h"

namespace GlobalV
{
	std::ofstream ofs_running;
}

/************************************************
 *  unit test of class Memory
 ***********************************************/

/*
 * - Tested Functions:
 *   - print_all
 *     - print memory consumed (> MB) in a
 *     - std::ofstream file
*/

class MemoryTest : public testing::Test
{
protected:
	// for capturing stdout
	std::string line;
	std::string name;
	// for output in file
	std::ofstream ofs;
	std::ifstream ifs;
	double factor = 1024*1024;
	double value;
	void TearDown()
	{
		remove("tmp");
	}
};

TEST_F(MemoryTest, printall)
{
	ofs.open("tmp");
	// total memory is an internal parameter and added inside the class Memory
	ModuleBase::Memory::record("item1",10.0*factor);
	ModuleBase::Memory::record("item1",12.0*factor);
	ModuleBase::Memory::record("item1",8.0*factor);
	ModuleBase::Memory::record("item2",5.0*factor);
	ModuleBase::Memory::record("item3",15.0*factor);
	ModuleBase::Memory::record("item4",0.8*factor);

	ModuleBase::Memory::print_all(ofs);
	ofs.close();
	ifs.open("tmp");
	getline(ifs,line);
	getline(ifs,line);
	EXPECT_THAT(line,testing::HasSubstr("NAME---------------|MEMORY(MB)--------"));
	ifs >> name >> value;
	EXPECT_THAT(name,testing::HasSubstr("total"));
	EXPECT_NEAR(value,32.8,1e-8);
	
	ifs >> name >> value;
	EXPECT_THAT(name,testing::HasSubstr("item3"));
	EXPECT_NEAR(value,15.0,1e-8);
	
	ifs >> name >> value;
	EXPECT_THAT(name,testing::HasSubstr("item1"));
	EXPECT_NEAR(value,12.0,1e-8);
	
	ifs >> name >> value;
	EXPECT_THAT(name,testing::HasSubstr("item2"));
	EXPECT_NEAR(value, 5.0,1e-8);

	getline(ifs,line);
	getline(ifs,line);
	EXPECT_THAT(line,testing::HasSubstr("< 1.0 MB has been ignored"));
	ifs.close();
}