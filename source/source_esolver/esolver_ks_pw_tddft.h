#ifndef ESOLVER_KS_PW_TDDFT_H
#define ESOLVER_KS_PW_TDDFT_H

#include "source_esolver/esolver_ks_pw.h"
#include "source_estate/module_pot/td_field_manager.h"
#include "source_hsolver/hsolver_pw_tddft.h"
#include "source_io/module_current/td_current_pw.h"
#include "source_pw/module_pwdft/td_history_pw.h"

#include <memory>

namespace ModuleESolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class ESolver_KS_PW_TDDFT : public ESolver_KS_PW<T, Device>
{
  public:
    ESolver_KS_PW_TDDFT();
    ~ESolver_KS_PW_TDDFT() override = default;

    /** @brief Initialize the field, propagation solver and time-step history. */
    void before_all_runners(BaseCell& basecell, const Input_para& inp) override;

  protected:
    /** @brief Use ground-state diagonalization for the initial electronic step. */
    std::string diag_policy(const int istep) const override;

    /** @brief Require an unmixed endpoint prediction before accepting a corrected step. */
    ESolver_KS::DensityStage density_stage(const int istep, const int iter) const override;

    void before_scf(UnitCell& ucell, const int istep) override;

    /** @brief Print the TDDFT iteration header and initialize the PW iteration. */
    void iter_init(UnitCell& ucell, const int istep, const int iter) override;

    /** @brief Switch between ground-state diagonalization and real-time propagation. */
    void hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr) override;

    /** @brief Apply the convergence policy appropriate to the current stage. */
    void iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver) override;

    /** @brief Save the converged state as the origin of the next propagation step. */
    void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;

  private:
    using Real = typename GetTypeReal<T>::type;
    double q_unshifted_ = 0.0; // Global geometric momentum bound in inverse Bohr.
    std::unique_ptr<hsolver::HSolverPWTDDFT<T, Device>> td_solver_;

    std::shared_ptr<elecstate::TDFieldManager> td_field_manager_;
    ModuleIO::CurrentPW<Real, Device> current_output_;

    pw::TDHistoryPW<T, Device> history_;
    void prepare_td_step(const int istep);
};

} // namespace ModuleESolver
#endif
