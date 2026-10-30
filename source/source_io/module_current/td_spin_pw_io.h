#ifndef TD_SPIN_PW_IO_H
#define TD_SPIN_PW_IO_H
#include "source_pw/module_rt/td_obs_pw.h"

#include <iosfwd>
#include <string>
namespace ModuleIO
{
/**
 * @brief Write globally reduced endpoint data on the designated output rank only.
 * @param step Zero-based internal electronic step, written as step + 1.
 * @param time_fs Physical time in fs, written with six decimal places; the initial state is at zero.
 */
void write_td_spin_pw(const std::string& directory,
                      int step,
                      double time_fs,
                      const std::vector<std::string>& orbitals,
                      const std::string& mag_partition,
                      const std::string& orbital_dir,
                      double onsite_radius,
                      bool calc_torque,
                      bool spin_current,
                      const pw::TDSpinResultPW& spin,
                      const pw::TDVelResultPW& current,
                      std::ostream& log);
} // namespace ModuleIO
#endif
