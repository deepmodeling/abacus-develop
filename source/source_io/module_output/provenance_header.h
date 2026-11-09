#ifndef MODULEIO_PROVENANCE_HEADER_H
#define MODULEIO_PROVENANCE_HEADER_H

#include "source_base/formatter.h"
#include "source_main/version.h"

#include <ctime>
#include <string>

namespace ModuleIO
{

/// Local wall-clock timestamp "YYYY-MM-DD HH:MM:SS" for the provenance
/// header of output files. Call once per frame and pass the result to
/// provenance_header() so all files of the same frame carry one timestamp.
inline std::string local_timestamp()
{
    const std::time_t now = std::time(nullptr);
    std::tm tm_local;
    localtime_r(&now, &tm_local);
    char time_buf[64];
    std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", &tm_local);
    return std::string(time_buf);
}

/// First two comment lines shared by per-step output files
/// (band*.txt, wfk*_nao.txt, eband_*.txt, ...):
/// "# ABACUS version: ..." and "# Written at ...".
/// Callers append their own step/energy lines after this block.
inline std::string provenance_header(const std::string& time_stamp)
{
    return FmtCore::format("# ABACUS version: %s\n# Written at %s\n",
                           VERSION,
                           time_stamp.c_str());
}

/// Two provenance lines followed by the unified "# Ionic_Step N" marker.
inline std::string prov_step(const int istep_one_based, const std::string& time_stamp)
{
    return FmtCore::format("# ABACUS version: %s\n# Written at %s\n# Ionic_Step %d\n",
                           VERSION,
                           time_stamp.c_str(),
                           istep_one_based);
}

/// Inline provenance for formats such as Cube that reserve a fixed number of
/// header lines. Callers append this to an existing comment line.
inline std::string prov_inline(const std::string& time_stamp)
{
    return FmtCore::format("  # ABACUS version: %s  # Written at %s",
                           VERSION,
                           time_stamp.c_str());
}

} // namespace ModuleIO

#endif // MODULEIO_PROVENANCE_HEADER_H
