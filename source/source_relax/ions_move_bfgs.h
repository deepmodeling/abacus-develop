#ifndef IONS_MOVE_BFGS_H
#define IONS_MOVE_BFGS_H

#include "bfgs_basic.h"
#include "relax_criteria.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"

#include <fstream>
#include <iostream>
#include <vector>
class Ions_Move_BFGS : public BFGS_Basic
{
  public:
    Ions_Move_BFGS();
    ~Ions_Move_BFGS();

    void allocate(void);
    void reset(void);
    bool start(
        UnitCell &ucell,
        const ModuleBase::matrix &force,
        const double &energy_in,
        const int istep,
        int &update_iter,
        std::ofstream &ofs,
        std::vector<double> &etot_info,
        const Relax_Criteria &criteria);

    //====================================================================
    // Test seam; see the equivalent block in BFGS_Basic. Production code
    // must keep using the private names directly.
    //====================================================================

    /// @brief whether allocate() has already run
    bool &get_init_done() { return init_done; }

    void bfgs_routine_for_testing(
        const int istep,
        int &update_iter,
        std::ofstream &ofs,
        std::vector<double> &etot_info,
        const std::string &out_level,
        const int test_relax_method)
    {
        bfgs_routine(istep, update_iter, ofs, etot_info, out_level, test_relax_method);
    }
    void restart_bfgs_for_testing(
        int &update_iter,
        std::ofstream &ofs,
        const int test_relax_method)
    {
        restart_bfgs(update_iter, ofs, test_relax_method);
    }

  private:
    bool init_done;
    void bfgs_routine(
        const int istep,
        int &update_iter,
        std::ofstream &ofs,
        std::vector<double> &etot_info,
        const std::string &out_level,
        const int test_relax_method);
    void restart_bfgs(
        int &update_iter,
        std::ofstream &ofs,
        const int test_relax_method);
};

#endif
