#ifndef SOURCE_MD_POTENTIAL_H
#define SOURCE_MD_POTENTIAL_H

#include "source_base/matrix.h"

class BaseCell;
struct Input_para;

namespace ModuleMD
{
// The potential supplies Ry, Ry/Bohr, and Ry/Bohr^3 values. For MDCell it
// writes owned and ghost forces to the atoms; for UnitCell it fills force.
class Potential
{
  public:
    virtual ~Potential() {}
    virtual void initialize(BaseCell& cell, const Input_para& inp) = 0;
    virtual void compute(BaseCell& cell,
                         double& energy,
                         ModuleBase::matrix& force,
                         ModuleBase::matrix& stress) = 0;
};
} // namespace ModuleMD

#endif
