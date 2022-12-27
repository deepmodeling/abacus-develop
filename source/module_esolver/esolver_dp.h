#ifndef ESOLVER_DP_H
#define ESOLVER_DP_H

#include "./esolver.h"
#ifdef __DPMDC
#define __DPMD 1
#include "deepmd/deepmd.hpp"
using deepmd::hpp::DeepPot;
#else
#ifdef __DPMD
#include "deepmd/DeepPot.h"
using deepmd::DeepPot;
#endif
#endif

namespace ModuleESolver
{

    class ESolver_DP : public ESolver
    {
    public:
#ifdef __DPMD
        ESolver_DP(std::string pot_file) : dp(pot_file)
        {
            classname = "ESolver_DP";
        }
#else
        ESolver_DP(std::string pot_file)
        {
            classname = "ESolver_DP";
        }
#endif

        void Init(Input& inp, UnitCell& cell) override;
        void Run(const int istep, UnitCell& cell) override;
        void cal_Energy(double& etot) override;
        void cal_Force(ModuleBase::matrix& force) override;
        void cal_Stress(ModuleBase::matrix& stress) override;
        void postprocess() override;

        //--------------temporary----------------------------
#ifdef __DPMD
        DeepPot dp;
#endif
        std::vector<double> cell;
        std::vector<int> atype;
        std::vector<double> coord;
        double dp_potential;
        ModuleBase::matrix dp_force;
        ModuleBase::matrix dp_virial;
        //---------------------------------------------------
    };
}
#endif
