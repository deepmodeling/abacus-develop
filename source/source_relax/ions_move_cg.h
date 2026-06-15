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

    void allocate(const int dim);
    bool start(UnitCell &ucell, const ModuleBase::matrix &force, const double &etot, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info, std::vector<std::string>& relax_method);

    static double RELAX_CG_THR;

  private:
    std::vector<double> pos0;
    std::vector<double> grad0;
    std::vector<double> cg_grad0;
    std::vector<double> move0;
    double e0 = 0.0;

    bool sd = false;
    bool trial = false;
    int ncggrad = 0;
    int nbrent = 0;
    double fa = 0.0;
    double fb = 0.0;
    double fc = 0.0;
    double xa = 0.0;
    double xb = 0.0;
    double xc = 0.0;
    double xpt = 0.0;
    double steplength = 0.0;
    double fmax = 0.0;
};

#endif