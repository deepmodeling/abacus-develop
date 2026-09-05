#ifndef ABACUS_D3_EVALUATOR_H
#define ABACUS_D3_EVALUATOR_H

#include "vdwd3_types.h"

#include <string>

namespace vdw
{
namespace d3
{

bool evaluate(const Structure& structure,
              const Parameters& parameters,
              const Cutoffs& cutoffs,
              bool derivatives,
              Result& result,
              std::string& error);

} // namespace d3
} // namespace vdw

#endif // ABACUS_D3_EVALUATOR_H
