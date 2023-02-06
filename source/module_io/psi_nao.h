#ifndef PSI_NAO_H
#define PSI_NAO_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_elecstate/elecstate.h"

// mohan add 2010-09-09
namespace ModuleIO
{
	void write_psi_nao(const std::string &name, double** ctot, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg);
	void write_psi_nao_complex(const std::string &name, std::complex<double>** ctot, const int &ik, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg);

	void distri_psi_nao_complex(std::complex<double>** ctot, std::complex<double> **cc);

    void distri_psi_nao_new(double** ctot, const int& is,
        const Parallel_Orbitals* ParaV, psi::Psi<double>* psid);
    void distri_psi_nao_complex_new(std::complex<double>** ctot, const int& ik,
        const Parallel_Orbitals* ParaV, psi::Psi<std::complex<double>>* psi);

    int read_psi_nao(
        double** ctot, 
        const int& is,
        const Parallel_Orbitals* ParaV, 
        psi::Psi<double>* psid,
        elecstate::ElecState* pelec);

    int read_psi_nao_complex(
        std::complex<double>** ctot, 
        const int& ik,
        const Parallel_Orbitals* ParaV, 
        psi::Psi<std::complex<double>>* psi,
        elecstate::ElecState* pelec);
}

#endif
