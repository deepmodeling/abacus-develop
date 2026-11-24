#ifndef SOURCE_MD_FORCEFIELD_H
#define SOURCE_MD_FORCEFIELD_H

class BaseCell;
struct Input_para;
namespace ModuleBase
{
class matrix;
}

namespace ModuleMD
{
// The force field supplies Ry, Ry/Bohr, and Ry/Bohr^3 values. For MDCell it
// writes owned and ghost forces to the atoms; for UnitCell it fills force.
class ForceField
{
  public:
    virtual ~ForceField() {}
    virtual void initialize(BaseCell &cell, const Input_para &inp) = 0;
    virtual void compute(
        BaseCell &cell,
        double &energy,
        ModuleBase::matrix &force,
        ModuleBase::matrix &stress) = 0;
};
} // namespace ModuleMD

#endif
