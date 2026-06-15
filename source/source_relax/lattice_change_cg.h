#ifndef LATTICE_CHANGE_CG_H
#define LATTICE_CHANGE_CG_H

#include <fstream>
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include "cg_base.h"
#include <vector>

class Lattice_Change_CG : public CG_Base
{

  public:
    Lattice_Change_CG();
    ~Lattice_Change_CG() = default;

    void allocate(void);
    bool start(UnitCell &ucell, const ModuleBase::matrix &stress_in, const double &etot, std::ofstream& ofs, std::vector<double>& etot_info);

  private:
    std::vector<double> lat0;
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