//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-07
//==========================================================
#ifndef EXPORT_H
#define EXPORT_H
using namespace std;
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "../module_base/global_variable.h"

#ifdef __MPI
#include <mpi.h>
#endif

template <class T>
void AUTO_SET(ofstream &ofs,const string &name,const T &a)
{
	ofs<<" AUTO_SET "<<name<<" : "<<a << endl;
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
		if(MY_RANK == 0)
		{
			cout<<"\n Can not match : "<<a<<"  "<<b<<endl;
		}
#ifdef __MPI
		MPI_Finalize();
#endif
		exit(0);
	}
	//cout<<setw(12)<<a<<endl;
	return;
}

void IF_MATCH(const string &name,const string &name2);


#endif 
