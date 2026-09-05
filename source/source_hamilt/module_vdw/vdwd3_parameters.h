#ifndef ABACUS_D3_PARAMETERS_H
#define ABACUS_D3_PARAMETERS_H

#include "vdwd3_types.h"

#include <string>

namespace vdw
{
namespace d3
{

std::string canonicalize_method_name(const std::string& input);

bool lookup_parameters(const std::string& method,
                       Damping damping,
                       Parameters& parameters,
                       std::string& canonical_method);

} // namespace d3
} // namespace vdw

#endif // ABACUS_D3_PARAMETERS_H
