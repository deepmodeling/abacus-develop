#ifndef H_TDDFT_PW_H
#define H_TDDFT_PW_H

#include "pot_base.h"
#include "source_base/vector3.h"

#include <memory>
#include <vector>

namespace elecstate
{

class TDFieldManager;

/**
 * @brief Length-gauge RT-TDDFT potential backed by TDFieldManager.
 *
 * The class consumes the manager's current field samples to construct the
 * periodic real-space potential. Static members publish only gauge and sampled electric-field data
 * for LCAO and force consumers. Time advancement belongs to the ESolver.
 */
class H_TDDFT_pw : public PotBase
{
  public:
    /**
     * @brief Construct the time-dependent potential component.
     *
     * @param rho_basis_in Real-space grid on which the potential is evaluated.
     * @param ucell_in Unit cell defining Cartesian coordinates and lattice size.
     * @param field_manager Shared source of RT-TDDFT field state.
     */
    H_TDDFT_pw(const ModulePW::PW_Basis* rho_basis_in, const UnitCell* ucell_in, const std::shared_ptr<TDFieldManager>& field_manager);

    /** @brief Destroy the time-dependent potential component. */
    ~H_TDDFT_pw()
    {
    }

    /**
     * @brief Add the current length-gauge external potential to the fixed part.
     *
     * @param vl_pseudo Real-space fixed potential updated in place.
     */
    void cal_fixed_v(double* vl_pseudo) override;

    /**
     * @brief Publish sampled fields for LCAO and ionic-force consumers.
     *
     * @param manager Active RT-TDDFT field manager.
     */
    static void set_field_state(const TDFieldManager& manager);

    /**
     * @brief Compute ionic force of electric field.
     *
     * @param[in] cell Information of cell.
     * @param[out] fe Electric-field force, F = qE.
     */
    static void compute_force(const UnitCell& cell, ModuleBase::matrix& fe);

    /** @brief Spatial-gauge selector for the remaining LCAO consumers. */
    static int stype;

    /** @brief Sampled electric field in Hartree atomic units for LCAO consumers. */
    static ModuleBase::Vector3<double> efield_ha;

  private:
    static std::vector<double> global_vext_time;

    void cal_v_space_length(std::vector<double>& vext_space, int direction);
    double cal_v_space_length_potential(double coordinate) const;

    const UnitCell* ucell_ = nullptr;
    std::shared_ptr<TDFieldManager> field_manager_;
};

} // namespace elecstate

#endif
