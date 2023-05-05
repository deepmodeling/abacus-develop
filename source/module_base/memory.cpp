//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-18
//==========================================================
#include "memory.h"
#include "global_variable.h"
#include "module_base/parallel_reduce.h"

namespace ModuleBase
{
//    8 bit  = 1 Byte
// 1024 Byte = 1 KB
// 1024 KB   = 1 MB
// 1024 MB   = 1 GB
double Memory::total = 0.0;
int Memory::n_memory = 1000;
int Memory::n_now = 0;
bool Memory::init_flag =  false;

std::string *Memory::name;
double *Memory::consume;

Memory::Memory()
{
}

Memory::~Memory()
{
}

void Memory::record
(
	const std::string &name_in,
	const size_t &n_in
)
{
	if(!Memory::init_flag)
	{
		name = new std::string[n_memory];
		consume = new double[n_memory];
		for(int i=0;i<n_memory;i++)
		{
			consume[i] = 0.0;
		}
		Memory::init_flag = true;
	}

	int find = 0;
	for(find = 0; find < n_now; find++)
	{
		if( name_in == name[find] )
		{
			break;
		}
	}

	// find == n_now : found a new record.	
	if(find == n_now)
	{
		n_now++;
		name[find] = name_in;
	}
	if(n_now >= n_memory)
	{
		std::cout<<" Error! Too many memories has been recorded.";
		return;
	}

	const double factor = 1.0/1024.0/1024.0;
	double size_mb = n_in * factor;

	if(consume[find] < size_mb)
	{
		Memory::total += size_mb - consume[find];
		consume[find] = size_mb;
		if(consume[find] > 5)
		{
			print(find);
		}
	}

	return;
}

void Memory::print(const int find)
{
	GlobalV::ofs_running <<"\n Warning_Memory_Consuming allocated: "
	<<" "<<name[find]<<" "<<consume[find]<<" MB" << std::endl;
	return;
}


void Memory::finish(std::ofstream &ofs)
{
	print_all(ofs);
	if(init_flag)
	{
		delete[] name;
		delete[] consume;
		init_flag = false;
	}
	return;
}

void Memory::print_all(std::ofstream &ofs)
{
	if(!init_flag) return;

	const double small = 1.0; 
#ifdef __MPI
		Parallel_Reduce::reduce_double_all(Memory::total);
#endif
    ofs <<"\n NAME---------------|MEMORY(MB)--------" << std::endl;
	ofs <<std::setw(20)<< "total" << std::setw(15) <<std::setprecision(4)<< Memory::total << std::endl;
    
	bool *print_flag = new bool[n_memory];
	for(int i=0; i<n_memory; i++) print_flag[i] = false;
	
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
#ifdef __MPI
		Parallel_Reduce::reduce_double_all(consume[k]);
#endif
	    if ( consume[k] < small ) 
        {
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
	delete[] print_flag; //mohan fix by valgrind at 2012-04-02
	return;
}

}
