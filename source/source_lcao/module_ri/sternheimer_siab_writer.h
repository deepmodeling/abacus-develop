#ifndef STERNHEIMER_SIAB_WRITER_H
#define STERNHEIMER_SIAB_WRITER_H

#include "sternheimer_siab_data.h"

#include <complex>
#include <string>
#include <vector>

namespace module_ri
{
namespace sternheimer_siab
{

void write_v1(const std::string& path,
              double grid_volume_bohr3,
              const std::vector<PrimitiveBlock>& blocks,
              const std::vector<ReferenceRow>& rows,
              const std::vector<std::complex<double>>& overlap_s,
              const Provenance& provenance);

} // namespace sternheimer_siab
} // namespace module_ri

#endif
