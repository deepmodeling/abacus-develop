#ifndef MODULE_MD_LJ_H
#define MODULE_MD_LJ_H

#include "source_md/potential.h"
#include "source_base/vector3.h"

#include <iosfwd>
#include <string>
#include <vector>

namespace ModuleMD
{

class LJ : public Potential
{
  public:
    LJ(const std::string& output_dir, std::ofstream& running_log);

    void initialize(BaseCell& cell, const Input_para& inp) override;
    void compute(BaseCell& cell,
                 double& energy,
                 ModuleBase::matrix& force,
                 ModuleBase::matrix& stress) override;

    //====================================================================
    // Test seam.
    //
    // initialize() derives the LJ tables in three steps. The unit tests
    // drive each step on its own and check the table it produced, so the steps
    // and the tables are reachable here rather than by reinterpreting the access
    // specifiers.
    //
    // Production code must keep going through initialize(); nothing
    // outside the tests should call the *_for_testing() wrappers.
    //====================================================================

    /// @brief neighbour search radius derived from the largest cutoff
    double get_search_radius() const
    {
        return search_radius;
    }
    /// @brief per-type-pair cutoff radii
    const ModuleBase::matrix& get_lj_rcut() const
    {
        return lj_rcut;
    }
    /// @brief per-type-pair c6 coefficients
    const ModuleBase::matrix& get_lj_c6() const
    {
        return lj_c6;
    }
    /// @brief per-type-pair c12 coefficients
    const ModuleBase::matrix& get_lj_c12() const
    {
        return lj_c12;
    }
    /// @brief per-type-pair energy shift at the cutoff
    const ModuleBase::matrix& get_en_shift() const
    {
        return en_shift;
    }
    void rcut_search_radius_for_testing(const int& ntype, const std::vector<double>& rcut)
    {
        rcut_search_radius(ntype, rcut);
    }
    void set_c6_c12_for_testing(const int& ntype,
                                const int& rule,
                                const std::vector<double>& epsilon,
                                const std::vector<double>& sigma)
    {
        set_c6_c12(ntype, rule, epsilon, sigma);
    }
    void cal_en_shift_for_testing(const int& ntype, const bool& is_shift)
    {
        cal_en_shift(ntype, is_shift);
    }

  private:
    std::string output_dir_;
    std::ofstream& running_log_;
    double LJ_energy(const double& d, const int& i, const int& j) const;

    ModuleBase::Vector3<double> LJ_force(const ModuleBase::Vector3<double>& dr, const int& i, const int& j) const;

    void rcut_search_radius(const int& ntype, const std::vector<double>& rcut);

    void set_c6_c12(const int& ntype,
                    const int& rule,
                    const std::vector<double>& epsilon,
                    const std::vector<double>& sigma);

    void cal_en_shift(const int& ntype, const bool& is_shift);

    double search_radius = -1.0;
    ModuleBase::matrix lj_rcut;
    ModuleBase::matrix lj_c12;
    ModuleBase::matrix lj_c6;
    ModuleBase::matrix en_shift;

};
} // namespace ModuleMD
#endif
