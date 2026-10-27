#ifndef TD_CURRENT_PW_H
#define TD_CURRENT_PW_H

#include "source_pw/module_rt/td_obs_pw.h"

#include <string>

namespace ModuleIO
{
/** @brief Write reduced charge currents on the output rank; performs no collectives. */
void write_current_pw(const pw::TDVelResultPW& current, int istep, int nk_per_spin, bool out_current_k, const std::string& out_dir);
} // namespace ModuleIO
#endif
