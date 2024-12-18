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

        void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

        void runner(UnitCell& cell, const int istep) override;

        double cal_energy() override;

        void cal_force(UnitCell& ucell, ModuleBase::matrix& force) override;

        void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) override;

        void after_all_runners(UnitCell& ucell) override;

      private:
        double LJ_energy(const double& d, const int& i, const int& j);

        ModuleBase::Vector3<double> LJ_force(const ModuleBase::Vector3<double>& dr, const int& i, const int& j);

        void LJ_virial(const ModuleBase::Vector3<double>& force, const ModuleBase::Vector3<double>& dtau);

        void rcut_search_radius(const int& ntype, const std::vector<double>& rcut);

        void set_c6_c12(const int& ntype,
                        const int& rule,
                        const std::vector<double>& epsilon,
                        const std::vector<double>& sigma);

        void cal_en_shift(const int& ntype, const bool& is_shift);

        //--------------temporary----------------------------
        double search_radius=-1.0;
        ModuleBase::matrix lj_rcut;
        ModuleBase::matrix lj_c12;
        ModuleBase::matrix lj_c6;
        ModuleBase::matrix en_shift;

        double lj_potential;
        ModuleBase::matrix lj_force;
        ModuleBase::matrix lj_virial;
        //---------------------------------------------------
    };
}
#endif
