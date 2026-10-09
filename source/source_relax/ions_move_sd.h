#ifndef IONS_MOVE_SD_H
#define IONS_MOVE_SD_H

#include "relax_criteria.h"
#include <fstream>
#include <iostream>
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include <string>
#include <vector>

namespace ions_move_sd
{
/**
 * @brief Update the steepest-descent trust radius from the energy history.
 *
 * Depends only on its arguments and on the shared Ions_Move_Basic state, so it
 * is a free function rather than a member.
 *
 * @param istep Current ionic step index (must be >= 1)
 * @param etot_info Energy information array [etot, etot_p]
 * @param out_level Output verbosity level ("ie" prints the radius to stdout)
 */
void cal_tradius_sd(const int istep, std::vector<double>& etot_info, const std::string& out_level);
} // namespace ions_move_sd

class Ions_Move_SD
{
  public:
    Ions_Move_SD();
    ~Ions_Move_SD() = default;

    void allocate(void);
    bool start(UnitCell& ucell, const ModuleBase::matrix& force, const double& etot, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info, const Relax_Criteria& criteria);

    /// @brief Energy of the last accepted step.
    double get_energy_saved() const
    {
        return energy_saved;
    }
    /// @brief Atomic positions of the last accepted step (dimension: Ions_Move_Basic::dim).
    const std::vector<double>& get_pos_saved() const
    {
        return pos_saved;
    }
    /// @brief Normalized gradient of the last accepted step (dimension: Ions_Move_Basic::dim).
    const std::vector<double>& get_grad_saved() const
    {
        return grad_saved;
    }

  private:
    double energy_saved;
    std::vector<double> pos_saved;
    std::vector<double> grad_saved;
};

#endif
