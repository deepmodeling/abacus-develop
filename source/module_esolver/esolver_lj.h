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

        void before_all_runners(Input& inp, UnitCell& cell) override;

        void runner(const int istep, UnitCell& cell) override;

        double cal_energy() override;

        void cal_force(ModuleBase::matrix& force) override;

        void cal_stress(ModuleBase::matrix& stress) override;

        void after_all_runners() override;

      private:

        double LJ_energy(const double d, const int i, const int j);

        ModuleBase::Vector3<double> LJ_force(const ModuleBase::Vector3<double> dr, const int i, const int j);

        void LJ_virial(const ModuleBase::Vector3<double>& force,
            const ModuleBase::Vector3<double>& dtau);

        //--------------temporary----------------------------
        int lj_rule;
        double lj_rcut;
        ModuleBase::matrix lj_c12;
        ModuleBase::matrix lj_c6;
        ModuleBase::matrix en_shift;

        double lj_potential;
        ModuleBase::matrix lj_force;
        ModuleBase::matrix lj_virial;
        UnitCell* ucell_; ///< pointer to the unitcell information
        //---------------------------------------------------
    };
}
#endif
