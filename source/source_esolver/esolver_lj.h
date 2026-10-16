#ifndef ESOLVER_LJ_H
#define ESOLVER_LJ_H

#include "esolver.h"

namespace ModuleESolver
{
class ESolver_LJ;
}

class MDCell;

namespace ModuleESolver
{

class ESolver_LJ : public ESolver
{
  public:
    ESolver_LJ()
    {
        classname = "ESolver_LJ";
    }

    void before_all_runners(BaseCell& cell, const Input_para& inp) override;

    void runner(BaseCell& cell, const int istep) override;

    double cal_energy() override;

    void cal_force(BaseCell& cell, ModuleBase::matrix& force) override;

    void cal_stress(BaseCell& cell, ModuleBase::matrix& stress) override;

        void after_all_runners(BaseCell& cell) override;

        void others(BaseCell& cell, const int istep) override;

    //====================================================================
    // Test seam.
    //
    // before_all_runners() derives the LJ tables in three steps. The unit tests
    // drive each step on its own and check the table it produced, so the steps
    // and the tables are reachable here rather than by reinterpreting the access
    // specifiers.
    //
    // Production code must keep going through before_all_runners(); nothing
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
    /// @brief computed lattice virials
    const ModuleBase::matrix& get_lj_virial() const
    {
        return lj_virial;
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
    double LJ_energy(const double& d, const int& i, const int& j) const;

    ModuleBase::Vector3<double> LJ_force(const ModuleBase::Vector3<double>& dr, const int& i, const int& j) const;

    void rcut_search_radius(const int& ntype, const std::vector<double>& rcut);

    void set_c6_c12(const int& ntype,
                    const int& rule,
                    const std::vector<double>& epsilon,
                    const std::vector<double>& sigma);

    void cal_en_shift(const int& ntype, const bool& is_shift);

    //--------------temporary----------------------------
    double search_radius = -1.0;
    ModuleBase::matrix lj_rcut;
    ModuleBase::matrix lj_c12;
    ModuleBase::matrix lj_c6;
    ModuleBase::matrix en_shift;

    double lj_potential = 0.0;
    ModuleBase::matrix lj_force;
    ModuleBase::matrix lj_virial;
    //---------------------------------------------------
};
} // namespace ModuleESolver
#endif
