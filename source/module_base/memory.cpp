//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-18
//==========================================================

//updated 2023/05/05 Wenfei

#include "memory.h"
#include "global_variable.h"
#include "module_base/parallel_reduce.h"
#include <numeric>   //std::iota
#include <algorithm> //std::sort

namespace ModuleBase
{
//    8 bit  = 1 Byte
// 1024 Byte = 1 KB
// 1024 KB   = 1 MB
// 1024 MB   = 1 GB
double Memory::total = 0.0;

std::vector<std::string> Memory::name;
std::vector<double> Memory::consume;

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
	const double factor = 1.0/1024.0/1024.0;
	double size_mb = n_in * factor;

	int find = -1;
	for(int i=0;i<name.size();i++)
	{
		if(name[i] == name_in) find = i;
	}

	if(find < 0)
	{
		name.push_back(name_in);
		consume.push_back(size_mb);
		Memory::total += size_mb;
		find = name.size()-1;
		if(consume[find] > 5)
		{
			print(find);
		}
	}
	else if(consume[find] < size_mb)
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
	GlobalV::ofs_running << std::setprecision(4);
	GlobalV::ofs_running <<"\n Warning_Memory_Consuming allocated: "
	<<" "<<name[find]<<" "<<consume[find]<<" MB" << std::endl;
	return;
}

void Memory::print_all(std::ofstream &ofs)
{
	const double small = 1.0; 
#ifdef __MPI
	Parallel_Reduce::reduce_double_all(total);
	for(int i=0;i<consume.size();i++)
	{
		Parallel_Reduce::reduce_double_all(consume[i]);
	}
#endif

    ofs <<"\n NAME---------------|MEMORY(MB)--------" << std::endl;
	ofs <<std::setw(20)<< "total" << std::setw(15) <<std::setprecision(4)<< Memory::total << std::endl;

	//sort 'consume' and keep track of the index
	std::vector<int> index;
	index.resize(consume.size());
	std::iota(index.begin(),index.end(),0);
	std::sort(index.begin(),index.end(), [&](int i,int j){return consume[i]>consume[j];} );

	for(int i=0;i<consume.size();i++)
	{
		if(consume[index[i]] < small) break;
		ofs << std::setw(20) << name[index[i]] << std::setw(15) << consume[index[i]] << std::endl;
	}

	ofs<<" -------------   < 1.0 MB has been ignored ----------------"<<std::endl;
    ofs<<" ----------------------------------------------------------"<<std::endl;

	return;
}

}
