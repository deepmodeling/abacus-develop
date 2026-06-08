#ifndef IONS_MOVE_SD_H
#define IONS_MOVE_SD_H

#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include <vector>

class Ions_Move_SD
{
  public:
    Ions_Move_SD();
    ~Ions_Move_SD() = default;

    void allocate(void);
    void start(UnitCell& ucell, const ModuleBase::matrix& force, const double& etot);

  private:
    double energy_saved;
    std::vector<double> pos_saved;
    std::vector<double> grad_saved;

    void cal_tradius_sd(void) const;
};

#endif
