#ifndef SOURCE_MD_MODULE_DP_DP_FORCEFIELD_H
#define SOURCE_MD_MODULE_DP_DP_FORCEFIELD_H

#include "source_md/forcefield.h"

#ifdef __DPMD
#ifdef __DPMDC
namespace deepmd
{
namespace hpp
{
class DeepPot;
}
} // namespace deepmd
#else
namespace deepmd
{
class DeepPot;
}
#endif
#endif

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

class UnitCell;
class MDCell;
class DomainDecomposition;

namespace ModuleMD
{
class DP : public ForceField
{
  public:
    DP(const std::string &pot_file, const std::string &output_dir, std::ofstream &running_log);
    ~DP() override;

    void initialize(BaseCell &cell, const Input_para &inp) override;
    void compute(
        BaseCell &cell,
        double &energy,
        ModuleBase::matrix &force,
        ModuleBase::matrix &stress) override;

  private:
    void initialize_type_map_(const std::vector<std::string> &type_labels);
    void compute_mdcell_(MDCell &cell, double &energy, ModuleBase::matrix &stress);

#ifdef __DPMD
#ifdef __DPMDC
    std::unique_ptr<deepmd::hpp::DeepPot> dp;
#else
    std::unique_ptr<deepmd::DeepPot> dp;
#endif
#endif
    std::string dp_file;
    std::string output_dir_;
    std::ofstream &running_log_;
    std::vector<int> md_type_to_dp_type_;
    std::vector<double> fparam;
    std::vector<double> aparam;
    std::unique_ptr<MDCell> distributed_cell_;
    std::unique_ptr<DomainDecomposition> decomposition_;
    bool aparam_nall_ = false;
    double rescaling = 1.0;
};
} // namespace ModuleMD

#endif
