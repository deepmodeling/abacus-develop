#ifndef IONS_MOVE_SD_H
#define IONS_MOVE_SD_H

#include <fstream>
#include <iostream>
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include <vector>

class Ions_Move_SD
{
  public:
    Ions_Move_SD();
    ~Ions_Move_SD() = default;

    void allocate(void);
    bool start(UnitCell& ucell, const ModuleBase::matrix& force, const double& etot, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info);

  private:
    double energy_saved;
    std::vector<double> pos_saved;
    std::vector<double> grad_saved;

    void cal_tradius_sd(const int istep, std::vector<double>& etot_info) const;
};

#endif
