#ifndef ESOLVER_LJ_H
#define ESOLVER_LJ_H

#include "./esolver.h"
//--------------temporary----------------------------
#include "../module_neighbor/sltk_grid_driver.h"
//---------------------------------------------------

namespace ModuleESolver
{

class ESolver_LJ: public ESolver
{
public:
    ESolver_LJ() : grid_neigh(GlobalV::test_deconstructor, GlobalV::test_grid_driver, GlobalV::test_grid)
    {
        classname = "ESolver_LJ";
    }
    
    void Init(Input &inp, UnitCell_pseudo &cell) override;
    void Run(const int istep, UnitCell_pseudo& cell) override;
    void cal_Energy(energy& en) override;
    void cal_Force(ModuleBase::matrix &force) override;
    void cal_Stress(ModuleBase::matrix& stress) override;


    double LJ_energy(const double d);
    ModuleBase::Vector3<double> LJ_force(const double d, 
        const ModuleBase::Vector3<double> dr);
    void LJ_virial(const ModuleBase::Vector3<double> &force, 
        const ModuleBase::Vector3<double> &dtau);

//--------------temporary----------------------------
    Grid_Driver grid_neigh;
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