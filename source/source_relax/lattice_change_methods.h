#ifndef LATTICE_CHANGE_METHODS_H
#define LATTICE_CHANGE_METHODS_H

#include <fstream>
#include <vector>
#include "lattice_change_basic.h"
#include "lattice_change_cg.h"

class Lattice_Change_Methods
{
  public:
    Lattice_Change_Methods();

    ~Lattice_Change_Methods();

    void allocate(void);

    void cal_lattice_change(const int &istep,
                            const int &stress_step,
                            const ModuleBase::matrix &stress,
                            const double &etot,
                            UnitCell &ucell,
                            std::ofstream& ofs);

    bool get_converged(void) const
    {
        return converged_;
    }

    double get_ediff(void) const
    {
        return etot_info_[0] - etot_info_[1];
    }

    double get_largest_grad(void) const
    {
        return Lattice_Change_Basic::largest_grad;
    }

  private:
    Lattice_Change_CG lccg;
    bool converged_ = false;
    std::vector<double> etot_info_{0.0, 0.0}; // [etot, etot_p]
};
#endif
