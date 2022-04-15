#ifndef DIAGOELPA_H
#define DIAGOELPA_H

#include "diagh.h"
#include "module_orbital/parallel_orbitals.h"

namespace ModuleHSolver
{

class DiagoElpa : public DiagH
{

public:
    void diag(
        ModuleHamilt::Hamilt* phm_in,
        ModulePsi::Psi<double> &psi,
        double *eigenvalue_in)override;

    void diag(
        ModuleHamilt::Hamilt* phm_in,
        ModulePsi::Psi<std::complex<double>> &psi,
        double *eigenvalue_in)override;

private:

    bool ifElpaHandle(const bool& newIteration, const bool& ifNSCF);

    static bool is_already_decomposed;


};

}

#endif