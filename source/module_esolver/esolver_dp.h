#ifndef ESOLVER_DP_H
#define ESOLVER_DP_H

#include "./esolver.h"

namespace ModuleESolver
{

class ESolver_DP: public ESolver
{
public:
    ESolver_DP()
    {
        tag = "ESolver_DP";
    }
    
    void Init(Input &inp, UnitCell_pseudo &cell) override;
    void Run(int istep, UnitCell_pseudo& cell) override;
    void cal_Energy(energy& en) override;
    void cal_Force(ModuleBase::matrix &force) override;
    void cal_Stress(ModuleBase::matrix &stress) override;

//--------------temporary----------------------------
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