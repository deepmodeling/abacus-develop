#ifndef ESOLVER_FORCEFIELDS_H
#define ESOLVER_FORCEFIELDS_H

#include "esolver.h"
#include <iosfwd>
#include <memory>

namespace ModuleMD
{
class ForceField;
}

namespace ModuleESolver
{
class esolver_forcefields : public ESolver
{
  public:
    esolver_forcefields(std::unique_ptr<ModuleMD::ForceField> potential, std::ofstream &running_log);
    ~esolver_forcefields() override;

    void before_all_runners(BaseCell &cell, const Input_para &inp) override;
    void runner(BaseCell &cell, int istep) override;
    double cal_energy() override;
    void cal_force(BaseCell &cell, ModuleBase::matrix &force) override;
    void cal_stress(BaseCell &cell, ModuleBase::matrix &stress) override;
    void after_all_runners(BaseCell &cell) override;

  private:
    std::unique_ptr<ModuleMD::ForceField> potential_;
    std::ofstream &running_log_;
    double energy_ = 0.0;
    ModuleBase::matrix force_;
    ModuleBase::matrix stress_;
};
} // namespace ModuleESolver

#endif
