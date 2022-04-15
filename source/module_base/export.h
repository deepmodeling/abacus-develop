//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-07
//==========================================================
#ifndef EXPORT_H
#define EXPORT_H
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "global_variable.h"

#ifdef __MPI
#include "mpi.h"
#endif
//these two function is not used yet!!!!
/*namespace ModuleBase
{

template <class T>
void ModuleBase::GlobalFunc::AUTO_SET(std::ofstream &ofs,const std::string &name,const T &a)
{
	ofs<<" AUTO_SET "<<name<<" : "<<a << std::endl;
	return;
}

//===================
// OUT
//===================

template <class T>
void IF_MATCH(const T &a,const T &b)
{
	if(a!=b)
	{
		if(GlobalV::MY_RANK == 0)
		{
			std::cout<<"\n Can not match : "<<a<<"  "<<b<<std::endl;
		}
#ifdef __MPI
		MPI_Finalize();
#endif
		exit(0);
	}
	//std::cout<<std::setw(12)<<a<<std::endl;
	return;
}

void IF_MATCH(const std::string &name,const std::string &name2);

}*/

#endif 
