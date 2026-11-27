#ifndef MDCELL_CIF_IO_H
#define MDCELL_CIF_IO_H

#include <string>

class MDCell;

namespace mdcell
{
/// Write one distributed MD geometry in CIF format, ordered by atom identity.
void print_cif_file(const MDCell& cell, const std::string& filename, const std::string& frame_comment);
}

#endif
