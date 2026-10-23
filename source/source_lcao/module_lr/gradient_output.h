#ifndef ABACUS_LR_GRADIENT_OUTPUT_H
#define ABACUS_LR_GRADIENT_OUTPUT_H
#include "source_base/matrix.h"
#include <iosfwd>
#include <vector>
class UnitCell;
namespace LR
{
void print_lr_force(const std::vector<ModuleBase::matrix>& force, std::ostream& ofs, int istate_begin);
ModuleBase::matrix average_forces(const std::vector<ModuleBase::matrix>& forces);
void print_grad_matrix(const std::vector<std::vector<ModuleBase::matrix>>& gradient,
    const std::vector<int>& group, const UnitCell& ucell, std::ofstream& ofs);
}
#endif
