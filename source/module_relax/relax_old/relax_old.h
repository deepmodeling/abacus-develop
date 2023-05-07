#ifndef RELAX_OLD_H
#define RELAX_OLD_H

#include "ions_move_methods.h"
#include "lattice_change_methods.h"
#include "module_cell/unitcell.h"

class Relax_old
{
  public:
    void init_relax(const int& natom);
    bool relax_step(const int& istep,
                    const double& energy,
                    UnitCell& ucell,
                    ModuleBase::matrix force,
                    ModuleBase::matrix stress,
                    int& force_step,
                    int& stress_step);

  private:
    Ions_Move_Methods IMM;
    Lattice_Change_Methods LCM;

    // seperate force_stress function first
    bool if_do_relax(const UnitCell& ucell);
    bool if_do_cellrelax(const UnitCell& ucell);
    bool do_relax(const int& istep,
                  const ModuleBase::matrix& ionic_force,
                  const double& total_energy,
                  UnitCell& ucell,
                  int& jstep);
    bool do_cellrelax(const int& istep,
                      const int& stress_step,
                      const ModuleBase::matrix& stress,
                      const double& total_energy,
                      UnitCell& ucell);
};

#endif