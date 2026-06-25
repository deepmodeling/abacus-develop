#ifndef IONS_MOVE_METHODS_H
#define IONS_MOVE_METHODS_H

#include <fstream>
#include <iostream>
#include <vector>
#include "ions_move_basic.h"
#include "ions_move_bfgs.h"
#include "ions_move_cg.h"
#include "ions_move_sd.h"
#include "ions_move_bfgs2.h"
#include "ions_move_lbfgs.h"

class Ions_Move_Methods
{
  public:
    Ions_Move_Methods();
    ~Ions_Move_Methods();

    void allocate(const int &natom, const std::string& relax_method_0, const std::string& relax_method_1);
    void cal_movement(const int &istep,
                      const int &force_step,
                      const ModuleBase::matrix &f,
                      const double &etot,
                      UnitCell &ucell,
                      std::ofstream& ofs,
                      std::vector<std::string>& relax_method);
    void reset_after_cell_change(const std::vector<std::string>& relax_method, std::ofstream& ofs);

    bool get_converged() const
    {
        return converged_;
    }

    double get_ediff() const
    {
        return etot_info_[0] - etot_info_[1];
    }
    double get_largest_grad() const
    {
        return Ions_Move_Basic::largest_grad;
    }
    double get_trust_radius() const
    {
        return Ions_Move_Basic::trust_radius;
    }
    int get_update_iter() const
    {
        return update_iter_;
    }

  private:
    Ions_Move_BFGS bfgs;
    Ions_Move_CG cg;
    Ions_Move_SD sd;
    Ions_Move_BFGS2 bfgs_trad;
    Ions_Move_LBFGS lbfgs;
    bool converged_ = false;
    int update_iter_ = 0;
    std::vector<double> etot_info_{0.0, 0.0}; // [etot, etot_p]
};
#endif
