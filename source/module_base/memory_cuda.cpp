//==========================================================
// AUTHOR : Zhang Xiaoyang
// DATE : 2024-2-20
//==========================================================

// This is absolutely a copy from memory.cpp

#include "memory_cuda.h"
#include "global_variable.h"
#include "module_base/parallel_reduce.h"

namespace ModuleBase
{
//    8 bit  = 1 Byte
// 1024 Byte = 1 KB
// 1024 KB   = 1 MB
// 1024 MB   = 1 GB
double Memory_CUDA::total = 0.0;

int Memory_CUDA::n_memory = 1000;
int Memory_CUDA::n_now = 0;
bool Memory_CUDA::init_flag =  false;

std::string *Memory_CUDA::name;
std::string *Memory_CUDA::class_name;
double *Memory_CUDA::consume;

Memory_CUDA::Memory_CUDA()
{
}

Memory_CUDA::~Memory_CUDA()
{
}
	

double Memory_CUDA::record
(
 	const std::string &class_name_in,
	const std::string &name_in,
	const long &size,
	const bool accumulate
)
{
	if(!Memory_CUDA::init_flag)
	{
		name = new std::string[n_memory];
		class_name = new std::string[n_memory];
		consume = new double[n_memory];
		for(int i=0;i<n_memory;i++)
		{
			consume[i] = 0.0;
		}
		Memory_CUDA::init_flag = true;
	}

	int find = 0;
	for(find = 0; find < n_now; find++)
	{
		if( name_in == name[find] )
		{
			break;
		}
	}

	if(find == n_now)
	{
		n_now++;
		name[find] = name_in;
		class_name[find] = class_name_in;
	}
	if(n_now >= n_memory)
	{
		std::cout<<" Error! Too many memories required.";
		return 0.0;
	}

	consume[find] += size / 1024.0 / 1024.0;
	Memory_CUDA::total += size / 1024.0 / 1024.0;

	//if(consume[find] > 5)
	//{
	//	print(find);
	//}
	return consume[find];
}

void Memory_CUDA::print(const int find)
{
	GlobalV::ofs_running <<"\n Warning_GPU_Memory_Consuming allocated: "
	<<" "<<name[find]<<" "<<consume[find]<<" MB" << std::endl;
	return;
}


void Memory_CUDA::finish(std::ofstream &ofs)
{
	print_all(ofs);
	if(init_flag)
	{
		delete[] name;
		delete[] class_name;
		delete[] consume;
		init_flag = false;
	}
	return;
}

void Memory_CUDA::print_all(std::ofstream &ofs)
{
	if(!init_flag) return;

	const double small = 1.0; 
#ifdef __MPI
		Parallel_Reduce::reduce_all(Memory_CUDA::total);
#endif
    ofs <<"\n NAME---------------|GPU_MEMORY(MB)----" << std::endl;
	ofs <<std::setw(20)<< "total" << std::setw(15) <<std::setprecision(4)<< Memory_CUDA::total << std::endl;
    
	bool *print_flag = new bool[n_memory];
	for(int i=0; i<n_memory; i++) print_flag[i] = false;
	

	for (int i=0; i<n_memory; i++)
    {
#ifdef __MPI
		Parallel_Reduce::reduce_all(consume[i]);
#endif
	}

	for (int i=0; i<n_memory; i++)
	{
		int k = 0;
		double tmp = -1.0;
		for(int j=0; j<n_memory; j++)
		{
			if(print_flag[j])
			{
				continue;
			}
			else if(tmp < consume[j])
			{
				k = j;
				tmp = consume[j];
			}
		}
		print_flag[k] = true;
		if ( consume[k] < small ){
			continue;
		}
		else
		{
			ofs << std::setw(20) << name[k]
            << std::setw(15) << consume[k] << std::endl;
		}

	}
	ofs<<" -------------   < 1.0 MB has been ignored ----------------"<<std::endl;
    ofs<<" ----------------------------------------------------------"<<std::endl;
	delete[] print_flag;
	return;
}

}
