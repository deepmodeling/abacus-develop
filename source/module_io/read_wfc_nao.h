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
    void distri_wfc_nao(double** ctot, const int& is, const int& nb2d,const int& nbands_g,
                        const int& nlocal_g, const Parallel_Orbitals* ParaV, psi::Psi<double>* psid);
    void distri_wfc_nao_complex(std::complex<double>** ctot, const int& ik,const int& nb2d, const int& nbands_g,
                                const Parallel_Orbitals* ParaV, psi::Psi<std::complex<double>>* psi);

    int read_wfc_nao(
        double** ctot, 
        const int& is,
        const bool& gamma_only_local,
        const int& nb2d,
        const int& nbands_g,
        const int& nlocal_g,
        const std::string& global_readin_dir,
        const Parallel_Orbitals* ParaV, 
        psi::Psi<double>* psid,
        elecstate::ElecState* pelec);

    int read_wfc_nao_complex(
        std::complex<double>** ctot, 
        const int& ik,
        const int& nb2d,
        const int& nbands_g,
        const int& nlocal_g,
        const std::string& global_readin_dir,
        const ModuleBase::Vector3<double> kvec_c,
        const Parallel_Orbitals* ParaV, 
        psi::Psi<std::complex<double>>* psi,
        elecstate::ElecState* pelec);
}

#endif
