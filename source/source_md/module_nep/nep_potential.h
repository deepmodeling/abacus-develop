#ifndef SOURCE_MD_MODULE_NEP_POTENTIAL_H
#define SOURCE_MD_MODULE_NEP_POTENTIAL_H

#include "source_md/potential.h"

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

class NEP;
class UnitCell;

namespace ModuleMD
{
class NEPPotential : public Potential
{
  public:
    NEPPotential(const std::string& pot_file, const std::string& output_dir, std::ofstream& running_log);
    ~NEPPotential() override;

    void initialize(BaseCell& cell, const Input_para& inp) override;
    void compute(BaseCell& cell,
                 double& energy,
                 ModuleBase::matrix& force,
                 ModuleBase::matrix& stress) override;

  private:
    void type_map(const UnitCell& cell);
    void initialize_type_map_(const std::vector<std::string>& type_labels);

#ifdef __NEP
    std::unique_ptr<::NEP> nep;
#endif
    std::string nep_file;
    std::string output_dir_;
    std::ofstream& running_log_;
    std::vector<int> atype;
    std::vector<int> md_type_to_nep_type_;
    std::vector<double> _e;
    std::vector<double> _f;
    std::vector<double> _v;
};
} // namespace ModuleMD

#endif
