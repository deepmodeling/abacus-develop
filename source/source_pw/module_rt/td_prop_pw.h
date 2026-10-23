#ifndef PW_TD_PROP_PW_H
#define PW_TD_PROP_PW_H
#include "source_estate/module_pot/td_field_manager.h"
#include "source_hsolver/module_rt_pw/cn_solver_pw.h"
#include "source_hsolver/module_rt_pw/exp_action.h"
#include "source_hsolver/module_rt_pw/td_orbital_ops.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_pw/module_rt/td_hamilt_pw.h"
#include "source_pw/module_rt/td_history_pw.h"

namespace pw
{
enum class PropStage
{
    predictor,
    corrector
};
enum class PropKind
{
    cn,
    etrs
};

/** @brief Resolved input controls, passed explicitly to numerical services. */
struct TDPropOpts
{
    PropKind prop_kind = PropKind::cn;
    hsolver::TDOrbOpts orb;
    hsolver::PWLinearOptions cn;
    hsolver::ExpOpts exp;
};
TDPropOpts td_prop_options(const Input_para& input, int npol, const std::vector<int>& global_k);

/** @brief Borrowed state for a single predictor or endpoint iteration. */
template <typename T, typename Device>
struct PropArgs
{
    TDHamiltPW<T, Device>& hamilt;
    const TDHistoryPW<T, Device>& history;
    elecstate::Potential& potential;
    const elecstate::TDFieldManager& field;
    int elec_step;
    int scf_iter;
    bool log_detail;
};

/** @brief Method-specific wavefunction map; density and acceptance stay above it. */
template <typename T, typename Device>
class PWProp
{
  public:
    virtual ~PWProp() = default;
    virtual void begin_trial() = 0;
    virtual void discard_trial() = 0;
    virtual void invalidate_basis() = 0;
    virtual void evolve(PropArgs<T, Device>& args, psi::Psi<T, Device>* current) = 0;
};

/** @brief Sole PW propagation entry point, independent of the selected method. */
template <typename T, typename Device>
class TDPropPW
{
  private:
    const ModulePW::PW_Basis_K& basis_;
    hsolver::TDOrbOps<T, Device> orb_ops_;
    std::unique_ptr<PWProp<T, Device>> prop_;

  public:
    TDPropPW(const ModulePW::PW_Basis_K& basis, const TDPropOpts& prop_opts, const hsolver::diag_comm_info& comm, std::ostream& log);
    PropStage stage(int scf_iter) const;
    void begin_trial();
    void finish_trial();
    void discard_trial();
    void invalidate_basis();
    hsolver::TDOrbOps<T, Device>& orbitals()
    {
        return orb_ops_;
    }
    /** @brief Restore the endpoint sample before returning to density construction. */
    void evolve(hamilt::HamiltPW<T, Device>& hamilt,
                elecstate::Potential& potential,
                const TDHistoryPW<T, Device>& history,
                const elecstate::TDFieldManager& field,
                psi::Psi<T, Device>* current,
                int elec_step,
                int scf_iter,
                bool log_detail);
};
} // namespace pw
#endif
