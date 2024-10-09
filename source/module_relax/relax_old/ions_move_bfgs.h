#ifndef IONS_MOVE_BFGS_H
#define IONS_MOVE_BFGS_H

#include "bfgs_basic.h"
#include "module_base/matrix.h"
#include "module_cell/unitcell.h"
class Ions_Move_BFGS : public BFGS_Basic
{
  public:
    Ions_Move_BFGS();
    ~Ions_Move_BFGS();

    void allocate(void);
    void start(UnitCell& ucell, const ModuleBase::matrix& force, const double& energy_in);

  private:
    bool init_done;
    void bfgs_routine(const double& lat0);
    void restart_bfgs(const double& lat0);
    bool first_step=true;   // If it is the first step of the relaxation. The pos is only generated from ucell in the first step, and in the following steps, the pos is generated from the previous step.
};

#endif
