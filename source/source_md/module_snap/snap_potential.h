#ifndef MODULE_MD_SNAP_POTENTIAL_H
#define MODULE_MD_SNAP_POTENTIAL_H

#include "source_md/potential.h"
#include "source_md/module_snap/snap_calculator.h"

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace ModuleMD
{

class SnapPotential : public Potential
{
  public:
    SnapPotential(const std::string& coefficient_file,
                 const std::string& parameter_file,
                 std::ofstream& running_log,
                 int rank);

    void initialize(BaseCell& cell, const Input_para& inp) override;
    void compute(BaseCell& cell,
                 double& energy,
                 ModuleBase::matrix& force,
                 ModuleBase::matrix& stress) override;

  private:
    void initialize(const std::vector<std::string>& labels);

    std::string coefficient_file_;
    std::string parameter_file_;
    std::ofstream& running_log_;
    const int rank_;
    ModuleMD::SNAP::ModelData model_;
    std::unique_ptr<ModuleMD::SNAP::Calculator> calculator_;
};

} // namespace ModuleMD

#endif // MODULE_MD_SNAP_POTENTIAL_H
