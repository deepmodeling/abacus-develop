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
    int istep = 0;
    double etot = 0;
    int force_step = 1;
    int stress_step = 1;

    ModuleBase::matrix force_;
    ModuleBase::matrix stress_;

    Relax rl;
    IonCellOptimizer rl_old;

    void init_relax(const int nat, const Input_para& inp);
    void iter_info(const Input_para& inp);
    void esolve(ModuleESolver::ESolver* p_esolver, UnitCell& ucell);
    bool relax_step(ModuleESolver::ESolver* p_esolver, UnitCell& ucell, const Input_para& inp);
    void stru_out(UnitCell& ucell, const Input_para& inp);
    void json_out(ModuleESolver::ESolver* p_esolver, UnitCell& ucell, const Input_para& inp);
    bool stop_cond(bool stop);
    void final_out(UnitCell& ucell, const Input_para& inp);
};

#endif
