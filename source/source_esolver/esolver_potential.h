#ifndef ESOLVER_POTENTIAL_H
#define ESOLVER_POTENTIAL_H

#include "esolver.h"
#include <iosfwd>
#include <memory>

namespace ModuleMD
{
class Potential;
}

namespace ModuleESolver
{
class ESolver_Potential : public ESolver
{
  public:
    ESolver_Potential(std::unique_ptr<ModuleMD::Potential> potential, std::ofstream& running_log);
    ~ESolver_Potential() override;

    void before_all_runners(BaseCell& cell, const Input_para& inp) override;
    void runner(BaseCell& cell, int istep) override;
    double cal_energy() override;
    void cal_force(BaseCell& cell, ModuleBase::matrix& force) override;
    void cal_stress(BaseCell& cell, ModuleBase::matrix& stress) override;
    void after_all_runners(BaseCell& cell) override;

  private:
    std::unique_ptr<ModuleMD::Potential> potential_;
    std::ofstream& running_log_;
    double energy_ = 0.0;
    ModuleBase::matrix force_;
    ModuleBase::matrix stress_;
};
} // namespace ModuleESolver

#endif
