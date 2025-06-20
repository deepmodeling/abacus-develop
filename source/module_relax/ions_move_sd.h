#ifndef IONS_MOVE_SD_H
#define IONS_MOVE_SD_H

#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
class Ions_Move_SD
{
  public:
    Ions_Move_SD();
    ~Ions_Move_SD();

    void allocate(void);
    void start(UnitCell& ucell, const ModuleBase::matrix& force, const double& etot);

  private:
    double energy_saved;
    double* pos_saved;
    double* grad_saved;

    void cal_tradius_sd(void) const;
};

#endif
