#ifndef WAVEFUNC_IN_PW_H
#define WAVEFUNC_IN_PW_H

#include "source_base/complexmatrix.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/realarray.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_pw/hamilt_pwdft/structure_factor.h"

//---------------------------------------------------
// FUNCTION: expand the local basis sets into plane
// wave basis sets
//---------------------------------------------------
namespace Wavefunc_in_pw
{

	void make_table_q(
		const UnitCell &ucell,
		std::vector<std::string> &orbital_files, 
		ModuleBase::realArray &table_local);

	void integral(
		const UnitCell& ucell,
		const int meshr, // number of mesh points 
		const double *psir,
		const double *r,
		const double *rab, 
		const int &l, 
		double* table);
	
	//mohan add 2010-04-20
	double smearing(
		const double &energy_x,
		const double &ecut,
		const double &beta);

    void produce_local_basis_in_pw(const UnitCell& ucell,
								   const int& ik,
                                   const ModulePW::PW_Basis_K* wfc_basis,
                                   const Structure_Factor& sf,
                                   ModuleBase::ComplexMatrix& psi,
                                   const ModuleBase::realArray& table_local);

}
#endif
