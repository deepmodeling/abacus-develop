//==========================================================
// AUTHOR : mohan
// DATE : 2021-01-30
//==========================================================
#ifndef PRINT_INFO_H
#define PRINT_INFO_H

#include "module_base/timer.h"
#include "module_cell/unitcell.h"
#include "module_cell/klist.h"

class Print_Info
{
	public:
	
	Print_Info();
	~Print_Info();

	// print out to screen about the readin parameters
	static void setup_parameters(UnitCell &ucell, K_Vectors &kv);

	static void print_time(time_t &time_start, time_t &time_finish);
    static void print_scf(const int &istep, const int &iter);
    static void print_screen(const int &stress_step, const int &force_step, const int &istep);

};

#endif
