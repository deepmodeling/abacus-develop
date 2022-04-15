#ifndef DIAGH_H
#define DIAGH_H

#include "module_base/complexmatrix.h"
#include "module_hamilt/hamilt.h"
#include "module_psi/psi.h"

namespace ModuleHSolver
{

class DiagH
{
    public:

    //virtual void init()=0;

    virtual void diag(
        ModuleHamilt::Hamilt* phm_in,
        ModulePsi::Psi<std::complex<double>> &psi,
        double *eigenvalue_in)=0;
    
    virtual void diag(
        ModuleHamilt::Hamilt* phm_in,
        ModulePsi::Psi<double> &psi,
        double *eigenvalue_in){return;}

};

}

#endif