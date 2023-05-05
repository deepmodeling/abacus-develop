#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <fstream>
#include <cstdio>

namespace GlobalV
{
	std::ofstream ofs_running;
}

/************************************************
 *  unit test of class Memory
 ***********************************************/

/**
 * - Tested Functions:
 *   - print_all
 *     - print memory consumed (> MB) in a
 *     - std::ofstream file
 */

#define private public
#include "../memory.h"

class MemoryTest : public testing::Test
{
protected:
	// definition according to ../memory_psi.cpp
	double factor = 1.0 / 1024.0 / 1024.0; // MB
	double complex_matrix_mem = 2*sizeof(double) * factor; // byte to MB
	double double_mem = sizeof(double) * factor;
	double int_mem = sizeof(int) * factor;
	double bool_mem = sizeof(bool) * factor;
	double float_mem = sizeof(float) * factor;
	double short_mem = sizeof(short) * factor;
	int n = 1024;
	// for capturing stdout
	std::string output;
	// for output in file
	std::ofstream ofs;
	std::ifstream ifs;
	void TearDown()
	{
		remove("tmp");
	}
};

TEST_F(MemoryTest, Constructor)
{
	EXPECT_NO_THROW(ModuleBase::Memory mem);
}

TEST_F(MemoryTest, printall)
{
	ofs.open("tmp");
	// total memory is an internal parameter and added inside the class Memory
	ModuleBase::Memory::record("Rrho",1024*1024);
	ModuleBase::Memory::record("drho",1024*1024);
	ModuleBase::Memory::record("evc",1024*1024);
	ModuleBase::Memory::print_all(ofs);
	ofs.close();
	ifs.open("tmp");
	getline(ifs,output);
	getline(ifs,output);
	EXPECT_THAT(output,testing::HasSubstr("MEMORY(MB)"));
	ifs.close();
}

TEST_F(MemoryTest, finish)
{
	*ModuleBase::Memory::name = "tmp_name";
	*ModuleBase::Memory::consume = 100.0;
	ModuleBase::Memory::init_flag = true;
	ofs.open("tmp");
	// total memory is an internal parameter and added inside the class Memory
	ModuleBase::Memory::record("Rrho",1024*1024);
	EXPECT_NO_THROW(ModuleBase::Memory::finish(ofs));
	ofs.close();
	EXPECT_FALSE(ModuleBase::Memory::init_flag);
}
#undef private
