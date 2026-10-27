#ifndef CTRL_OUTPUT_TD_PW_H
#define CTRL_OUTPUT_TD_PW_H

#include "source_pw/module_rt/td_obs_pw.h"

#include <ostream>

namespace ModuleIO
{
/** @brief Write precomputed endpoint observables on the designated output rank. */
void ctrl_output_td_pw(const pw::TDVelResultPW& currents,
                       const pw::TDSpinResultPW& moments,
                       const pw::TDVelRequestPW& request,
                       const pw::TDSpinOptsPW& spin_opts,
                       const std::vector<std::string>& orbital_fn,
                       int step,
                       double dt_ha,
                       int nspin,
                       int nkstot,
                       int rank,
                       const std::string& directory,
                       std::ostream& log);
} // namespace ModuleIO
#endif
