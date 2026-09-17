#ifndef IONS_MOVE_BFGS_H
#define IONS_MOVE_BFGS_H

#include <fstream>
#include <iostream>
#include <vector>
#include "bfgs_basic.h"
#include "relax_criteria.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
class Ions_Move_BFGS : public BFGS_Basic
{
  public:
    Ions_Move_BFGS();
    ~Ions_Move_BFGS();

    void allocate(void);
    void reset(void);
    bool start(UnitCell& ucell, const ModuleBase::matrix& force, const double& energy_in, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info, const Relax_Criteria& criteria);

    //====================================================================
    // Test seam; see the equivalent block in BFGS_Basic. Production code
    // must keep using the private names directly.
    //====================================================================

    /// @brief whether allocate() has already run
    bool& get_init_done()
    {
        return init_done;
    }
    /// @brief whether this is the first step of the relaxation
    bool& get_first_step()
    {
        return first_step;
    }

    void bfgs_routine_for_testing(const double& lat0,
                                  const int istep,
                                  int& update_iter,
                                  std::ofstream& ofs,
                                  std::vector<double>& etot_info,
                                  const std::string& out_level,
                                  const int test_relax_method)
    {
        bfgs_routine(lat0, istep, update_iter, ofs, etot_info, out_level, test_relax_method);
    }
    void restart_bfgs_for_testing(const double& lat0,
                                  int& update_iter,
                                  std::ofstream& ofs,
                                  const int test_relax_method)
    {
        restart_bfgs(lat0, update_iter, ofs, test_relax_method);
    }

  private:
    bool init_done;
    void bfgs_routine(const double& lat0, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info, const std::string& out_level, const int test_relax_method);
    void restart_bfgs(const double& lat0, int& update_iter, std::ofstream& ofs, const int test_relax_method);
    bool first_step=true;   // If it is the first step of the relaxation. The pos is only generated from ucell in the first step, and in the following steps, the pos is generated from the previous step.
};

#endif
