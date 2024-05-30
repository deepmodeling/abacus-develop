#ifndef WRITE_WFC_NAO_H
#define WRITE_WFC_NAO_H
#include <complex>
#include <vector>

#include "module_base/matrix.h"
#include "module_base/vector3.h"

// mohan add 2010-09-09
namespace ModuleIO
{
	void write_wfc_nao(const std::string &name, std::vector<std::vector<double>>& ctot, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, bool writeBinary = false);
	void write_wfc_nao_complex(const std::string &name, std::vector<std::vector<std::complex<double>>>& ctot, const int &ik, const ModuleBase::Vector3<double> &kvec_c, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, bool writeBinary = false);
}

#endif
