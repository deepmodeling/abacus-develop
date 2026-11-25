#include "source_io/module_current/td_current_pw.h"

#include "source_base/tool_quit.h"

#include <cstdint>
#include <fstream>
#include <iomanip>

namespace ModuleIO
{

void write_current_pw(const pw::TDVelResultPW& current,
                      const int istep,
                      const int nk_per_spin,
                      const bool out_current_k,
                      const std::string& out_dir)
{
    // The initial endpoint starts a new series; later endpoints append.
    const std::ios::openmode mode = istep == 0 ? std::ios::out | std::ios::trunc : std::ios::app;
    // Always write the total current when this routine is enabled.
    std::string filename_tot = out_dir + "current_tot.txt";
    std::ofstream fout_tot;
    fout_tot.open(filename_tot, mode);
    if (!fout_tot)
    {
        ModuleBase::WARNING_QUIT("write_current_pw", "Cannot open " + filename_tot);
    }
    fout_tot << std::right << std::setw(6) << istep + 1 << std::scientific << std::setprecision(12);
    for (double value: current.total)
    {
        fout_tot << "  " << std::setw(19) << value;
    }
    fout_tot << std::endl;
    fout_tot.close();
    if (!fout_tot)
    {
        ModuleBase::WARNING_QUIT("write_current_pw", "Failed to write " + filename_tot);
    }

    // Optionally write the contribution from each k point.
    if (out_current_k)
    {
        for (int ik = 0; ik < static_cast<int>(current.per_k.size() / 3); ++ik)
        {
            // Use a one-based spin index in the output filename.
            int is = ik / nk_per_spin + 1;
            int k_idx = ik % nk_per_spin + 1;

            std::string filename_k = out_dir + "current_s" + std::to_string(is) + "k" + std::to_string(k_idx) + ".txt";
            std::ofstream fout_k;
            fout_k.open(filename_k, mode);
            if (!fout_k)
            {
                ModuleBase::WARNING_QUIT("write_current_pw", "Cannot open " + filename_k);
            }
            const std::int64_t offset = static_cast<std::int64_t>(ik) * 3;
            fout_k << std::right << std::setw(6) << istep + 1 << std::scientific << std::setprecision(12);
            for (int direction = 0; direction < 3; ++direction)
            {
                fout_k << "  " << std::setw(19) << current.per_k[offset + direction];
            }
            fout_k << std::endl;
            fout_k.close();
            if (!fout_k)
            {
                ModuleBase::WARNING_QUIT("write_current_pw", "Failed to write " + filename_k);
            }
        }
    }
}

} // namespace ModuleIO
