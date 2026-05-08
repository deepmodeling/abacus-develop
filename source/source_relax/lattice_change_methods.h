#ifndef LATTICE_CHANGE_METHODS_H
#define LATTICE_CHANGE_METHODS_H

#include "lattice_change_basic.h"
#include "lattice_change_cg.h"

class Lattice_Change_Methods
{
  public:
    Lattice_Change_Methods ();

    ~Lattice_Change_Methods ();

    void allocate ();

    void cal_lattice_change (const int& istep,
                             const int& stress_step,
                             const ModuleBase::matrix& stress,
                             const double& etot,
                             UnitCell& ucell);

    bool
        get_converged () const
    {
        return Lattice_Change_Basic::converged;
    }

    double
        get_ediff () const
    {
        return Lattice_Change_Basic::ediff;
    }

    double
        get_largest_grad () const
    {
        return Lattice_Change_Basic::largest_grad;
    }

  private:
    Lattice_Change_CG lccg;
};
#endif
