#ifndef MODULEIO_SPIN_TAG_H
#define MODULEIO_SPIN_TAG_H

#include <string>

namespace ModuleIO
{

/// Running-log suffix naming the collinear spin channel of an output.
/// Returns " (spin up  )" for nspin==2 and is==0, " (spin down)" for nspin==2
/// and is==1, and an empty string otherwise (nspin==1, nspin==4, or is<0).
inline std::string make_spin_tag(const int is, const int nspin)
{
    if (nspin != 2 || is < 0)
    {
        return "";
    }
    if (is == 0)
    {
        return " (spin up  )";
    }
    return " (spin down)";
}

/// Compact spin-channel label used in per-spin output files
/// (bands1.txt, wfk1s1_nao.txt, ...). Unlike make_spin_tag, which pads "up"
/// so both channels share one width in the running log, this form is
/// unpadded because it heads a file that holds a single channel.
/// Returns " (spin up)", " (spin down)" or an empty string.
inline std::string make_compact_spin_tag(const int is, const int nspin)
{
    if (nspin != 2 || is < 0)
    {
        return "";
    }
    if (is == 0)
    {
        return " (spin up)";
    }
    return " (spin down)";
}

} // namespace ModuleIO

#endif
