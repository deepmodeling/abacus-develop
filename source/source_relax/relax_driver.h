#ifndef RELAX_DRIVER_H
#define RELAX_DRIVER_H

#include "source_cell/unitcell.h"
#include "source_esolver/esolver.h"
#include "relax_sync.h"
#include "relax_nsync.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_base/matrix.h"

class Relax_Driver
{

  public:
    Relax_Driver(){};
    ~Relax_Driver(){};

    void relax_driver(ModuleESolver::ESolver* p_esolver, 
            UnitCell& ucell, 
            const Input_para& inp);

  private:
    int force_step = 1;
    int stress_step = 1;

    Relax rl;
    IonCellOptimizer rl_old;

    void init_relax(const int nat, const Input_para& inp);
    void iter_info(const int istep, const Input_para& inp);
    void esolve(const int istep, ModuleESolver::ESolver* p_esolver, UnitCell& ucell, ModuleBase::matrix& force, ModuleBase::matrix& stress, double& etot);
    bool relax_step(const int istep, ModuleESolver::ESolver* p_esolver, UnitCell& ucell, const Input_para& inp, const ModuleBase::matrix& force, const ModuleBase::matrix& stress, const double etot);
    void stru_out(const int istep, UnitCell& ucell, const Input_para& inp);
    void json_out(ModuleESolver::ESolver* p_esolver, UnitCell& ucell, const Input_para& inp, const ModuleBase::matrix& force, const ModuleBase::matrix& stress);
    void final_out(const int istep, UnitCell& ucell, const Input_para& inp);
};

#endif
