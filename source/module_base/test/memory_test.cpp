#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <fstream>
#include <cstdio>

namespace GlobalV
{
	std::ofstream ofs_running;
    int KPAR = 1;
} // namespace GlobalV

/************************************************
 *  unit test of class Memory
 ***********************************************/

/**
 * - Tested Functions:
 *   - calculate_mem
 *     - caculate memory consumed for various types
 *     - of data
 *   - record
 *     - record memory consumed during running
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

TEST_F(MemoryTest, CalculateMem)
{
	// three types of complex number
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"ModuleBase::ComplexMatrix"),
			n * complex_matrix_mem);
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"complexmatrix"),
			n * complex_matrix_mem);
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"cdouble"), 
			n * complex_matrix_mem);

	// three types of integral number
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"bool"), n * bool_mem);
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"int"), n * int_mem);
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"short"), n * short_mem);

	// three types of float point number
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"real"), n * double_mem);
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"double"), n * double_mem);
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"float"), n * float_mem);

	// vector with 3 double float point number
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"ModuleBase::Vector3<double>"),
			n * 3 * double_mem);

	// test a struct AtomLink defined in module_neighbor/sltk_grid.h
	// AtomLink defined as class FAtom (module_neighbor/sltk_atom.h)
	// which includes 3 double and 2 int numbers
	EXPECT_EQ(ModuleBase::Memory::calculate_mem(n,"AtomLink"), 
			n * (int_mem * 2 + 3 * double_mem));

	// test types of data not defined
	testing::internal::CaptureStdout();
	ModuleBase::Memory::calculate_mem(n,"Exception");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("not this type in memory storage"));
}

TEST_F(MemoryTest,Record)
{
	// abacus record mem in MB
	// for double
	double mem = ModuleBase::Memory::record("wavefunc","evc",1024*1024,"double");
	EXPECT_EQ(mem,double_mem/factor);
	// for cdouble
	mem = ModuleBase::Memory::record("wavefunc","evc",1024*1024,"cdouble");
	EXPECT_EQ(mem,complex_matrix_mem/factor);
	// for int
	mem = ModuleBase::Memory::record("wavefunc","evc",1024*1024,"int");
	EXPECT_EQ(mem,int_mem/factor);
	// for bool
	mem = ModuleBase::Memory::record("wavefunc","evc",1024*1024,"bool");
	EXPECT_EQ(mem,bool_mem/factor);
	// for float
	mem = ModuleBase::Memory::record("wavefunc","evc",1024*1024,"float");
	EXPECT_EQ(mem,float_mem/factor);
	// for short
	mem = ModuleBase::Memory::record("wavefunc","evc",1024*1024,"short");
	EXPECT_EQ(mem,short_mem/factor);
	// for Vector3
	mem = ModuleBase::Memory::record("wavefunc","evc",1024*1024,"ModuleBase::Vector3<double>");
	EXPECT_EQ(mem,double_mem/factor*3);
	// for AtomLink
	mem = ModuleBase::Memory::record("wavefunc","evc",1024*1024,"AtomLink");
	EXPECT_EQ(mem,double_mem/factor*3+int_mem/factor*2);
}

TEST_F(MemoryTest, printall)
{
	ofs.open("tmp");
	// total memory is an internal parameter and added inside the class Memory
	ModuleBase::Memory::record("Charge_Mixing","Rrho",1024*1024,"ModuleBase::Vector3<double>");
	ModuleBase::Memory::record("Charge_Mixing","drho",1024*1024,"AtomLink");
	ModuleBase::Memory::record("wavefunc","evc",1024*1024,"float");
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
	*ModuleBase::Memory::class_name = "tmp_class_name";
	*ModuleBase::Memory::consume = 100.0;
	ModuleBase::Memory::init_flag = true;
	ofs.open("tmp");
	// total memory is an internal parameter and added inside the class Memory
	ModuleBase::Memory::record("Charge_Mixing","Rrho",1024*1024,"ModuleBase::Vector3<double>");
	EXPECT_NO_THROW(ModuleBase::Memory::finish(ofs));
	ofs.close();
	EXPECT_FALSE(ModuleBase::Memory::init_flag);
}
#undef private
