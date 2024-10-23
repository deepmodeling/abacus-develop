#pragma once
#include "module_base/matrix.h"
// this namespace used to store global function for some stress operation
namespace StressTools
{
// set upper matrix to whole matrix
void stress_fill(const double& lat0_, const double& omega_, ModuleBase::matrix& stress_matrix);
} // namespace StressTools