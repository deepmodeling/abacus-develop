#ifndef ESOLVER_FP_H
#define ESOLVER_FP_H
#include <fstream>

#include "esolver.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_elecstate/elecstate.h"
#include "module_psi/psi.h"
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
        virtual void init_after_vc(Input& inp, UnitCell& cell);    // liuyu add 2023-03-09
        // Hamilt* phamilt;

        elecstate::ElecState* pelec = nullptr;
        Charge chr;
        //--------------temporary----------------------------
        // this is the interface of non-self-consistant calculation
        virtual void nscf(){};

        // wavefunction coefficients
        psi::Psi<std::complex<double>>* psi = nullptr;
        psi::Psi<double>* psid = nullptr;

      private:
        void print_rhofft(Input& inp, ofstream &ofs);
    };
}

#endif
