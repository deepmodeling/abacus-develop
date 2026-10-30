#include "source_io/module_ctrl/ctrl_output_td_pw.h"

#include "source_base/constants.h"
#include "source_io/module_current/td_current_pw.h"
#include "source_io/module_current/td_spin_pw_io.h"

namespace ModuleIO
{
void ctrl_output_td_pw(const pw::TDVelResultPW& currents,
                       const pw::TDSpinResultPW& moments,
                       const pw::CurrentRequest& request,
                       const std::string& mag_partition,
                       const std::string& orbital_dir,
                       double onsite_radius,
                       bool calc_torque,
                       const std::vector<std::string>& orbital_fn,
                       int step,
                       double dt_ha,
                       int nspin,
                       int nkstot,
                       int rank,
                       const std::string& directory,
                       std::ostream& log)
{
    if (rank != 0)
    {
        return;
    }
    if (request.charge_total)
    {
        const int spin_mult = nspin == 2 ? 2 : 1;
        const int nk = nkstot / spin_mult;
        write_current_pw(currents, step, nk, request.charge_per_k, directory);
    }
    if (nspin == 4)
    {
        const double time = step * dt_ha * ModuleBase::AU_to_FS;
        write_td_spin_pw(directory,
                         step,
                         time,
                         orbital_fn,
                         mag_partition,
                         orbital_dir,
                         onsite_radius,
                         calc_torque,
                         request.spin_current,
                         moments,
                         currents,
                         log);
    }
}
} // namespace ModuleIO
