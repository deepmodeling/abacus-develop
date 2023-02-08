#ifndef WRITE_WFC_NAO_H
#define WRITE_WFC_NAO_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_elecstate/elecstate.h"

// mohan add 2010-09-09
namespace ModuleIO
{
	void write_wfc_nao(const std::string &name, double** ctot, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg);
	void write_wfc_nao_complex(const std::string &name, std::complex<double>** ctot, const int &ik, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg);
}

#endif
