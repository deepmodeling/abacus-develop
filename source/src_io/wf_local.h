#ifndef WF_LOCAL_H
#define WF_LOCAL_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "src_lcao/local_orbital_charge.h"

// mohan add 2010-09-09
namespace WF_Local
{
	void write_lowf(const std::string &name, double** ctot);
	void write_lowf_complex(const std::string &name, std::complex<double>** ctot, const int &ik);

	void distri_lowf(double** ctot, double **c);
	void distri_lowf_complex(std::complex<double>** ctot, std::complex<double> **cc);

    void distri_lowf_new(double** ctot, const int& is,
        Local_Orbital_wfc &lowf);
    void distri_lowf_complex_new(std::complex<double>** ctot, const int& ik,
        Local_Orbital_wfc &lowf);

    int read_lowf(double** ctot, const int& is,
        Local_Orbital_wfc &lowf);

    int read_lowf_complex(std::complex<double>** ctot, const int& ik,
        Local_Orbital_wfc &lowf);
}

#endif
