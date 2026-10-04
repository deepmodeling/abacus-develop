#ifndef IONS_MOVE_METHODS_H
#define IONS_MOVE_METHODS_H

#include <fstream>
#include <iostream>
#include <vector>
#include "ions_move_basic.h"
#include "relax_criteria.h"
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
                      std::vector<std::string>& relax_method,
                      const Relax_Criteria& criteria);
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

    //====================================================================
    // Test seam.
    //
    // reset_after_cell_change() is expected to clear the per-run state of
    // whichever method is active, so the test seeds that state through the
    // sub-optimiser it belongs to and checks it was cleared. get_converged()
    // and get_update_iter() above already cover the reads.
    //
    // Production code must keep using the private names directly.
    //====================================================================

    void set_converged(const bool value)
    {
        converged_ = value;
    }
    void set_update_iter(const int value)
    {
        update_iter_ = value;
    }
    /// @brief {etot, etot_p} of the current and previous step
    std::vector<double>& get_etot_info()
    {
        return etot_info_;
    }
    /// @brief the BFGS optimiser used by relax_method "bfgs"
    Ions_Move_BFGS& get_bfgs()
    {
        return bfgs;
    }
    /// @brief the traditional BFGS optimiser used by relax_method "bfgs_trad"
    Ions_Move_BFGS2& get_bfgs_trad()
    {
        return bfgs_trad;
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
