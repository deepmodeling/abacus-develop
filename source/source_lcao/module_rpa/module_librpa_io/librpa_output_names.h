#ifndef LIBRPA_OUTPUT_NAMES_H
#define LIBRPA_OUTPUT_NAMES_H

#include <string>

namespace RpaLriDetail
{

/// Prefix carried by every file the LibRPA producer writes. LibRPA reads
/// these files by name, and users see many other OUT.* files from the same
/// run, so the prefix keeps the producer's subset identifiable at a glance.
inline const char* librpa_output_prefix()
{
    return "librpa_";
}

/// Build one producer file name: <rpa_outdir>/librpa_<name>.
inline std::string librpa_output_path(const std::string& rpa_outdir, const std::string& name)
{
    return rpa_outdir + librpa_output_prefix() + name;
}

} // namespace RpaLriDetail

#endif
