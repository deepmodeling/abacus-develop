#ifndef ESOLVER_LJ_H
#define ESOLVER_LJ_H

#include "esolver.h"

namespace ModuleESolver
{

    class ESolver_LJ : public ESolver
    {
    public:
        ESolver_LJ()
        {
            classname = "ESolver_LJ";
        }

        void Init(Input& inp, UnitCell& cell) override;
        void Run(const int istep, UnitCell& cell) override;
        double cal_Energy() override;
        void cal_Force(ModuleBase::matrix& force) override;
        void cal_Stress(ModuleBase::matrix& stress) override;
        void postprocess() override;


        double LJ_energy(const double d);
        ModuleBase::Vector3<double> LJ_force(const double d,
            const ModuleBase::Vector3<double> dr);
        void LJ_virial(const ModuleBase::Vector3<double>& force,
            const ModuleBase::Vector3<double>& dtau);

        //--------------temporary----------------------------
        double lj_rcut;
        double lj_sigma;
        double lj_epsilon;
        double lj_potential;
        ModuleBase::matrix lj_force;
        ModuleBase::matrix lj_virial;
        //---------------------------------------------------
    };
}
#endif
