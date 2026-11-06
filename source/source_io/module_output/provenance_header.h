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

} // namespace ModuleIO

#endif // MODULEIO_PROVENANCE_HEADER_H
