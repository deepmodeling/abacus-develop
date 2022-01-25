//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-18
//==========================================================
#ifndef MEMORY_H
#define MEMORY_H

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
namespace ModuleBase
{

//    8 bit  = 1 Byte
// 1024 Byte = 1 KB
// 1024 KB   = 1 MB
// 1024 MB   = 1 GB

class Memory
{
public:

	Memory();
	~Memory();

	static double record
	(
	 	const std::string &class_name,
		const std::string &name,
		const long &n,
		const std::string &type,
		const bool accumulate = false
	);

	static double& get_total(void){return total;}
	static void finish(std::ofstream &ofs);
	static void print_all(std::ofstream &ofs);
	static void print(const int find_in);
	static double calculate_mem(const long &n,const std::string &type);

private:

	static double total;
	static std::string *name;
	static std::string *class_name;
	static double *consume;
	static int n_memory;
	static int n_now;
	static bool init_flag;

	static double complex_matrix_memory;//(16 Byte)
	static double double_memory;//(8 Byte)
	static double int_memory;//(4 Byte)
	static double bool_memory;
	static double short_memory;//(2 Byte)
	static double float_memory;//(4 Byte)
};

}

#endif
