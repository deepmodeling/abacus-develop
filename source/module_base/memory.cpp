//==========================================================
// AUTHOR : Zhang Xiaoyang
// DATE : 2024-2-26
// A reconfiguration version of original memory.cpp
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
double Memory::total_gpu = 0.0;
int Memory::complex_matrix_memory = 2*sizeof(double); // 16 byte
int Memory::double_memory = sizeof(double); // 8 byte
int Memory::int_memory = sizeof(int); // 4.0 Byte
int Memory::bool_memory = sizeof(bool); // 1.0 Byte
int Memory::float_memory = sizeof(float); // 4.0 Byte
int Memory::short_memory = sizeof(short); // 2.0 Byte

//std::string *Memory::name;
//std::string *Memory::class_name;
//double *Memory::consume;
std::map<std::string,double> Memory::name_mem_map;
std::map<std::string,double> Memory::name_mem_gpu_map;
std::map<std::string,bool> Memory::name_print_flag_map;
std::map<std::string,bool> Memory::name_print_flag_gpu_map;

bool Memory::init_flag =  false;
bool Memory::init_flag_gpu =  false;


Memory::Memory()
{
}

Memory::~Memory()
{
}

double Memory::calculate_mem(const long &n_in,const std::string &type)
{
	double n = static_cast<double>(n_in);
	double mem = 0.0;
	
	double factor = 1.0/1024.0/1024.0;
	double complex_matrix_mem = complex_matrix_memory * factor;
	double double_mem = double_memory * factor;
	double int_mem = int_memory * factor;
	double bool_mem = bool_memory * factor;
	double float_mem = float_memory * factor;
	double short_mem = short_memory * factor;

	if(type=="ModuleBase::ComplexMatrix" || type=="complexmatrix" || type=="cdouble")
	{
		mem = complex_matrix_mem;
	}
	else if(type=="real" || type=="double")
	{
		mem = double_mem;
	}
	else if(type=="int")
	{
		mem = int_mem;
	}
	else if(type=="bool")
	{
		mem = bool_mem;
	}
	else if(type=="short")
	{
		mem = short_mem;
	}
	else if(type=="float")
	{
		mem = float_mem;
	}
	else if(type=="AtomLink")
	{
		mem =  int_mem * 2 + double_mem * 3;
	}
	else if(type=="ModuleBase::Vector3<double>")
	{
		mem = 3 * double_mem;
	}
	else
	{
		std::cout<<"not this type in memory storage : "<<type << std::endl;
	}
	total += n * mem;	
	return n*mem;
}
	

double Memory::record
(
 	const std::string &class_name_in,
	const std::string &name_in,
	const long &n_in,
	const std::string &type,
	const bool accumulate
)
{
	std::map<std::string,double>::iterator iter;
	iter = name_mem_map.find(name_in);
	if(iter != name_mem_map.end()){
		name_mem_map[name_in] += Memory::calculate_mem(n_in,type);
	}else{
		name_mem_map[name_in] = Memory::calculate_mem(n_in,type);
		init_flag = true;
	}
	total += Memory::calculate_mem(n_in,type);
	name_print_flag_map[name_in] = false;
	return name_mem_map[name_in];
}

void Memory::record
(
	const std::string &name_in,
	const size_t &n_in,
	const bool accumulate
)
{
	const double factor = 1.0/1024.0/1024.0;
	std::map<std::string,double>::iterator iter;
	iter = name_mem_map.find(name_in);
	if(iter != name_mem_map.end()){
		name_mem_map[name_in] += n_in * factor;
	}else{
		name_mem_map[name_in] = n_in * factor;
		init_flag = true;
	}
	name_print_flag_map[name_in] = false;
	total += n_in * factor;
	return;
}

double Memory::record_gpu
(
	const std::string &class_name,
	const std::string &name_in,
	const long &size,
	const bool accumulate
)
{
	const double factor = 1.0/1024.0/1024.0;
	std::map<std::string,double>::iterator iter;
	iter = name_mem_gpu_map.find(name_in);
	if(iter != name_mem_gpu_map.end()){
		name_mem_gpu_map[name_in] += size * factor;
	}else{
		name_mem_gpu_map[name_in] = size * factor;
		init_flag_gpu = true;
	}
	total_gpu += size * factor;
	name_print_flag_gpu_map[name_in] = false;
	return name_mem_gpu_map[name_in];
}

void Memory::print(const int find)
{
	//GlobalV::ofs_running <<"\n Warning_Memory_Consuming allocated: "
	//<<" "<<name[find]<<" "<<consume[find]<<" MB" << std::endl;
	return;
}


void Memory::finish(std::ofstream &ofs)
{
	//print_all(ofs);
	//if(init_flag)
	//{
	//	delete[] name;
	//	delete[] class_name;
	//	delete[] consume;
	//	init_flag = false;
	//}
	return;
}


void Memory::print_all(std::ofstream &ofs){
	const double small = 1.0; 
	if (init_flag){
#ifdef __MPI
		Parallel_Reduce::reduce_all(Memory::total);
#endif
		ofs <<"\n NAME---------------|MEMORY(MB)--------" << std::endl;
	//	std::cout<<"\n"<<std::setw(41)<< " " <<std::setprecision(4)<<total;
		ofs <<std::setw(20)<< "total" << std::setw(15) <<std::setprecision(4)<< Memory::total << std::endl;

		std::map<std::string,double>::iterator iter;
		std::map<std::string,double>::iterator sub_iter;
		for(iter = name_mem_map.begin(); iter != name_mem_map.end(); iter++){
#ifdef __MPI
			Parallel_Reduce::reduce_all(iter->second);
#endif
		}
		for(iter = name_mem_map.begin(); iter != name_mem_map.end(); iter++){
			double tmp = -1.0;
			std::string current = "";
			for (sub_iter = name_mem_map.begin(); sub_iter != name_mem_map.end(); sub_iter++){
				if (name_print_flag_map[sub_iter->first]){
					continue;
				}
				else if(tmp < sub_iter->second){
					tmp = sub_iter->second;
					current = sub_iter->first;
				}
			}
			name_print_flag_map[current] = true;
			if ( name_mem_map[current] < small ){
				continue;
			}
			else
			{
				ofs << std::setw(20) << current
				<< std::setw(15) << name_mem_map[current] << std::endl;
			}
		}
		ofs<<" -------------   < 1.0 MB has been ignored ----------------"<<std::endl;
		ofs<<" ----------------------------------------------------------"<<std::endl;
	}


	if (init_flag_gpu){
#ifdef __MPI
		Parallel_Reduce::reduce_all(Memory::total_gpu);
#endif
		ofs <<"\n NAME---------------|MEMORY_GPU(MB)----" << std::endl;
	//	std::cout<<"\n"<<std::setw(41)<< " " <<std::setprecision(4)<<total;
		ofs <<std::setw(20)<< "total" << std::setw(15) <<std::setprecision(4)<< Memory::total_gpu << std::endl;

		std::map<std::string,double>::iterator iter;
		std::map<std::string,double>::iterator sub_iter;
		for(iter = name_mem_gpu_map.begin(); iter != name_mem_gpu_map.end(); iter++){
#ifdef __MPI
			Parallel_Reduce::reduce_all(iter->second);
#endif
		}
		for(iter = name_mem_gpu_map.begin(); iter != name_mem_gpu_map.end(); iter++){
			double tmp = -1.0;
			std::string current = "";
			for (sub_iter = name_mem_gpu_map.begin(); sub_iter != name_mem_gpu_map.end(); sub_iter++){
				if (name_print_flag_gpu_map[sub_iter->first]){
					continue;
				}
				else if(tmp < sub_iter->second){
					tmp = sub_iter->second;
					current = sub_iter->first;
				}
			}
			name_print_flag_gpu_map[current] = true;
			if ( name_mem_gpu_map[current] < small ){
				continue;
			}
			else
			{
				ofs << std::setw(20) << current
				<< std::setw(15) << name_mem_gpu_map[current] << std::endl;
			}
		}
		ofs<<" -------------   < 1.0 MB has been ignored ----------------"<<std::endl;
		ofs<<" ----------------------------------------------------------"<<std::endl;
	}

}

}
