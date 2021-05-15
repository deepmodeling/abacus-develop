//==========================================================
// AUTHOR : mohan
//==========================================================
#include "global_function.h"
#include "global_file.h"
#include "../src_parallel/parallel_common.h"

//==========================================================
// USE FILE timer.h
// ONLY :  output time after quit.
//==========================================================
#include "timer.h"
#include "memory.h"

#include <string>
#include <fstream>
#include <iostream>
using namespace std;

void NOTE(const string &words)
{
	return;
	if(ofs_running)
	{
		//ofs_running << " *********************************************************************************" << endl;
		ofs_running << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		ofs_running << " " << words << endl;
		ofs_running << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	}
}

void NEW_PART(const string &words)
{
	ofs_running << "\n ><><><><><><><><><><><><><><><><><><><><><><" << endl;
	ofs_running << "\n " << words << endl;
	ofs_running << "\n ><><><><><><><><><><><><><><><><><><><><><><\n" << endl;
	return;
}


//==========================================================
// GLOBAL FUNCTION :
// NAME : OUT( output date for checking )
//==========================================================
void OUT(ofstream &ofs,const string &name)
{
    ofs<<"\n"<<setw(18)<<name<<endl;
    return;
}

//==========================================================
// GLOBAL FUNCTION :
// NAME : MAKE_DIR( make dir ,using system function)
//==========================================================
void MAKE_DIR(const string &fn)
{
//	TITLE("global_function","MAKE_DIR");
    if (MY_RANK==0)
    {
        stringstream ss;
        ss << " test -d " << fn << " || mkdir " << fn ;
//----------------------------------------------------------
// EXPLAIN : 'system' function return '0' if success
//----------------------------------------------------------
        if ( system( ss.str().c_str() ) )
        {
            WARNING_QUIT( "MAKE_DIR", fn );
        }
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    return;
}


void DONE(ofstream &ofs,const string &description, const bool only_rank0)
{
    if (only_rank0)
    {
        if (MY_RANK==0)
        {
     //       ofs << " ---------------------------------------------------------------------------------\n";
            ofs << " DONE : " << description;
            ofs << " Time : "<< timer::print_until_now() << " (SEC)" ;
			ofs << endl << endl;
     //       ofs << "\n ---------------------------------------------------------------------------------\n";
        }
    }
    else
    {
     //   ofs << " ---------------------------------------------------------------------------------\n";
        ofs << " DONE : " << description;
        ofs << " Time : "<< timer::print_until_now() << " (SEC)" ;
		ofs << endl << endl;
     //   ofs << "\n ---------------------------------------------------------------------------------\n";
    }
//   	cout << "\n---------------------------------------------------------------------------------\n";
    cout << " DONE(" << setw(10) << timer::print_until_now() <<" SEC) : "<< description << endl;
//   	cout << "\n---------------------------------------------------------------------------------\n";
    return;
}

//==========================================================
// GLOBAL FUNCTION :
// NAME : TEST_LEVEL
// control the test_level
//==========================================================
void TEST_LEVEL(const string &name)
{
    bool disable = true;
    if (disable) return;

    if (name == "none")
    {
        test_wf = 0;
        test_potential = 0;
        test_charge = 0;
    }
    else if (name == "init_potential")
    {
        test_wf = 1;
        test_potential = 1;
        test_charge = 1;
    }
    else if (name == "init_read")
    {
        test_input = 1;
        test_winput = 1;
        test_kpoint = 1;
        test_atom = 1;
        test_unitcell = 1;
#ifndef __EPM
        test_pseudo_cell = 1;
#else
        test_epm_unitcell = 1;
#endif
    }
    else if (name == "pw_init")
    {
        test_pw = 1;

    }

    return;
}


bool SCAN_BEGIN(ifstream &ifs, const string &TargetName, const bool restart)
{
    string SearchName;
    bool find = false;
    if (restart)
    {
        ifs.clear();
        ifs.seekg(0);
    }
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> SearchName;
        if ( SearchName == TargetName)
        {
            find = true;
            break;
        }
    }
    if (!find)
    {
        ofs_warning<<" In SCAN_BEGIN, can't find: "<<TargetName<<" block."<<endl;
    }
    return find;
}


void SCAN_END(ifstream &ifs, const string &TargetName)
{
    string SearchName;
    ifs >> SearchName;
    if ( SearchName != TargetName)
    {
        ofs_warning<<" In SCAN_END, can't find: "<<TargetName<<" block."<<endl;
    }
    return;
}

void BLOCK_HERE( const string &description)
{
//	return;
	cout << "\n********************************************";
    cout << "\n Here is a Block, 1: go on 0: quit";
    cout << "\n " << description;
	cout << "\n********************************************" << endl;
    bool go_on = false;
	if(MY_RANK==0)
	{
    	cin >> go_on;
	}

#ifdef __MPI
	int swap = go_on;
	if(MY_RANK == 0)swap = go_on;
	MPI_Bcast(&swap, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(MY_RANK != 0)go_on = static_cast<bool>( swap );
#endif
	if(go_on)
	{
		return;
	}
	else
	{
		QUIT();
	}
}


void OUT_TIME(const string &name, time_t &start, time_t &end)
{
	double mini = difftime(end, start)/60.0;
	if(mini>0.1)
	{
		ofs_warning << setprecision(2);
		ofs_warning << " -------------------------------------------------------" << endl;
		ofs_warning << " NAME < " << name << " > = " << endl;
		ofs_warning << " -> " << ctime(&start) << " -> " << ctime(&end);	
		ofs_warning << " TIME = " << mini << " [Minutes]" << endl;
		ofs_warning << " -------------------------------------------------------" << endl;
		ofs_warning << setprecision(6);
	}
}

size_t MemAvailable()
{
	size_t mem_sum = 0;
	int i=0;
	ifstream ifs("/proc/meminfo");
	while(ifs.good())
	{
		string label, size, kB;
		ifs>>label>>size>>kB;
		if(label=="MemAvailable:")
			return std::stol(size);
		else if(label=="MemFree:" || label=="Buffers:" || label=="Cached:")
		{
			mem_sum += std::stol(size);
			++i;
		}
		if(i==3)
			return mem_sum;
	}
	throw runtime_error("read /proc/meminfo error in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
}
