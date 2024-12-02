#ifndef IONS_MOVE_METHODS_H
#define IONS_MOVE_METHODS_H

#include "ions_move_basic.h"
#include "ions_move_bfgs.h"
#include "ions_move_cg.h"
#include "ions_move_sd.h"
#include "bfgs.h"

class Ions_Move_Methods
{
  public:
    Ions_Move_Methods();
    ~Ions_Move_Methods();

    void allocate(const int &natom);
    // void cal_movement(const int &istep, const ModuleBase::matrix &f, const double &etot);
    void cal_movement(const int &istep,
                      const int &force_step,
                      const ModuleBase::matrix &f,
                      const double &etot,
                      UnitCell &ucell);

    bool get_converged() const
    {
        return Ions_Move_Basic::converged;
    }
    double get_ediff() const
    {
        return Ions_Move_Basic::ediff;
    }
    double get_largest_grad() const
    {
        return Ions_Move_Basic::largest_grad;
    }
    double get_trust_radius() const
    {
        return Ions_Move_Basic::trust_radius;
    }
    double get_update_iter() const
    {
        return Ions_Move_Basic::update_iter;
    }

  private:
    Ions_Move_BFGS bfgs;
    Ions_Move_CG cg;
    Ions_Move_SD sd;
    BFGS bfgs_trad;
};
#endif
