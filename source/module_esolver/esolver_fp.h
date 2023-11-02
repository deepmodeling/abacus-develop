#ifndef ESOLVER_FP_H
#define ESOLVER_FP_H
#include <fstream>

#include "esolver.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
// #include "hamilt.h"
namespace ModuleESolver
{
    class ESolver_FP : public ESolver
    {
    public:
        ModulePW::PW_Basis* pw_rho;
        ModulePW::PW_Basis* pw_rhod;    ///< dense grid for dense part of the charge density (for uspp)
        ModulePW::PW_Basis_Big* pw_big; ///< [temp] pw_basis_big class
        ESolver_FP();
        virtual ~ESolver_FP();
        virtual void Init(Input& inp, UnitCell& cell) override;
        virtual void init_after_vc(Input& inp, UnitCell& cell);    // liuyu add 2023-03-09
        // Hamilt* phamilt;

        elecstate::ElecState* pelec = nullptr;
        Charge chr;
        ModuleSymmetry::Symmetry symm;
        //--------------temporary----------------------------
        // this is the interface of non-self-consistant calculation
        virtual void nscf(){};

        Structure_Factor sf;
        K_Vectors kv;

      private:
        void print_rhofft(Input& inp, std::ofstream &ofs);
    };
}

#endif
