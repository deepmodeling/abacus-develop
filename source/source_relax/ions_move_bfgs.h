#ifndef IONS_MOVE_BFGS_H
#define IONS_MOVE_BFGS_H

#include <fstream>
#include <iostream>
#include <vector>
#include "bfgs_basic.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
class Ions_Move_BFGS : public BFGS_Basic
{
  public:
    Ions_Move_BFGS();
    ~Ions_Move_BFGS();

    void allocate(void);
    void reset(void);
    bool start(UnitCell& ucell, const ModuleBase::matrix& force, const double& energy_in, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info);

  private:
    bool init_done;
    void bfgs_routine(const double& lat0, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info);
    void restart_bfgs(const double& lat0, int& update_iter, std::ofstream& ofs);
    bool first_step=true;   // If it is the first step of the relaxation. The pos is only generated from ucell in the first step, and in the following steps, the pos is generated from the previous step.
};

#endif
