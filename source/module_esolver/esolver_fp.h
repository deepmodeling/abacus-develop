#ifndef ESOLVER_FP_H
#define ESOLVER_FP_H
#include "./esolver.h"
#include "../module_pw/pw_basis.h"
#include <fstream>
// #include "hamilt.h"
namespace ModuleESolver
{
    class ESolver_FP : public ESolver
    {
    public:
        ModulePW::PW_Basis* pw_rho;
        ESolver_FP();
        virtual ~ESolver_FP();
        virtual void Init(Input& inp, UnitCell& cell) override;
        // Hamilt* phamilt;
    private:
        void print_rhofft(Input& inp, ofstream &ofs);
    };
}

#endif
