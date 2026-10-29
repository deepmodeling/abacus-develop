#ifndef ESOLVER_SNAP_H
#define ESOLVER_SNAP_H

#include "esolver.h"
#include "source_md/snap/snap_calculator.h"

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace ModuleESolver
{

class ESolver_SNAP : public ESolver
{
  public:
    ESolver_SNAP(const std::string& coefficient_file,
                 const std::string& parameter_file,
                 std::ofstream& running_log);

    void before_all_runners(BaseCell& cell, const Input_para& inp) override;
    void runner(BaseCell& cell, int istep) override;
    double cal_energy() override;
    void cal_force(BaseCell& cell, ModuleBase::matrix& force) override;
    void cal_stress(BaseCell& cell, ModuleBase::matrix& stress) override;
    void after_all_runners(BaseCell& cell) override;

  private:
    void initialize(const std::vector<std::string>& labels);

    std::string coefficient_file_;
    std::string parameter_file_;
    std::ofstream& running_log_;
    ModuleMD::SNAP::ModelData model_;
    std::unique_ptr<ModuleMD::SNAP::Calculator> calculator_;
    double snap_potential_ = 0.0;
    ModuleBase::matrix snap_force_;
    ModuleBase::matrix snap_virial_;
};

} // namespace ModuleESolver

#endif // ESOLVER_SNAP_H
