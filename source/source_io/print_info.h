#ifndef PRINT_INFO_H
#define PRINT_INFO_H

#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"

namespace ModuleIO
{
// print out to screen about the readin parameters
void print_parameters(
	const UnitCell& ucell, 
	K_Vectors& kv,
    const Input_para& inp);

void print_time(time_t& time_start, time_t& time_finish);

void print_screen(const int& stress_step, const int& force_step, const int& istep);

//! Print charge density using FFT
void print_rhofft(ModulePW::PW_Basis* pw_rhod,
                  ModulePW::PW_Basis* pw_rho,
                  ModulePW::PW_Basis_Big* pw_big,
                  std::ofstream& ofs);

void print_wfcfft(const Input_para& inp, ModulePW::PW_Basis_K& pw_wfc, std::ofstream& ofs);

void print_kpar(const int &nks, const int &kpar_lcao);

} // namespace ModuleIO

#endif
