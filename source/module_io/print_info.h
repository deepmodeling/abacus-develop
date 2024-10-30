//==========================================================
// AUTHOR : mohan
// DATE : 2021-01-30
//==========================================================
#ifndef PRINT_INFO
#define PRINT_INFO

#include "module_base/timer.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_parameter/input_parameter.h"

namespace ModuleIO
{

// print out to screen about the readin parameters
void setup_parameters(UnitCell& ucell, K_Vectors& kv);
void print_time(time_t& time_start, time_t& time_finish);
void print_screen(const int& stress_step, const int& force_step, const int& istep);
//! Print charge density using FFT
void print_rhofft(ModulePW::PW_Basis* pw_rhod,
                  ModulePW::PW_Basis* pw_rho,
                  ModulePW::PW_Basis_Big* pw_big,
                  std::ofstream& ofs);
void print_wfcfft(const Input_para& inp, ModulePW::PW_Basis_K& pw_wfc, std::ofstream& ofs);

} // namespace ModuleIO

#endif
