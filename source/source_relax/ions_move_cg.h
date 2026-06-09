#ifndef IONS_MOVE_CG_H
#define IONS_MOVE_CG_H

#include <fstream>
#include <iostream>
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include "cg_base.h"
#include <vector>

class Ions_Move_CG : public CG_Base
{
  public:
    Ions_Move_CG();
    ~Ions_Move_CG() = default;

    void allocate(void);
    void start(UnitCell &ucell, const ModuleBase::matrix &force, const double &etot, std::ofstream& ofs);

    static double RELAX_CG_THR;
    int sd_step = 0;
    int cg_step = 0;

  private:
    std::vector<double> pos0;
    std::vector<double> grad0;
    std::vector<double> cg_grad0;
    std::vector<double> move0;
    double e0 = 0.0;
};

#endif