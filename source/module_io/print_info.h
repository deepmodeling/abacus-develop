//==========================================================
// AUTHOR : mohan
// DATE : 2021-01-30
//==========================================================
#ifndef PRINT_INFO
#define PRINT_INFO

#include "module_base/timer.h"
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"

namespace ModuleIO
{

// print out to screen about the readin parameters
void setup_parameters(UnitCell& ucell, K_Vectors& kv);
void print_time(time_t& time_start, time_t& time_finish);
void print_scf(const int& istep, const int& iter);
void print_screen(const int& stress_step, const int& force_step, const int& istep);

} // namespace ModuleIO

#endif
