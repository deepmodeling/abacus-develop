#ifndef READ_WFC_NAO_H
#define READ_WFC_NAO_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_elecstate/elecstate.h"

// mohan add 2010-09-09
namespace ModuleIO
{
    void distri_wfc_nao(double** ctot, const int& is,
        const Parallel_Orbitals* ParaV, psi::Psi<double>* psid);
    void distri_wfc_nao_complex(std::complex<double>** ctot, const int& ik,
        const Parallel_Orbitals* ParaV, psi::Psi<std::complex<double>>* psi);

    int read_wfc_nao(
        double** ctot, 
        const int& is,
        const Parallel_Orbitals* ParaV, 
        psi::Psi<double>* psid,
        elecstate::ElecState* pelec);

    int read_wfc_nao_complex(
        std::complex<double>** ctot, 
        const int& ik,
        const ModuleBase::Vector3<double> kvec_c,
        const Parallel_Orbitals* ParaV, 
        psi::Psi<std::complex<double>>* psi,
        elecstate::ElecState* pelec);
}

#endif
