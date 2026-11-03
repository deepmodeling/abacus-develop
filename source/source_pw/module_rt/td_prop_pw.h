#ifndef PW_TD_PROP_PW_H
#define PW_TD_PROP_PW_H
#include "source_estate/module_pot/td_field_manager.h"
#include "source_pw/module_rt/td_orbital_ops.h"
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

/** @brief Method-specific wavefunction map; density and acceptance stay above it. */
template <typename T, typename Device>
class PWProp
{
  public:
    virtual ~PWProp() = default;
    virtual void begin_trial() = 0;
    virtual void discard_trial() = 0;
    virtual void invalidate_basis() = 0;
    virtual void evolve(TDHamiltPW<T, Device>& hamilt,
                        const TDHistoryPW<T, Device>& history,
                        elecstate::Potential& potential,
                        const elecstate::TDFieldManager& field,
                        psi::Psi<T, Device>* current,
                        int elec_step,
                        int scf_iter,
                        bool log_detail) = 0;
};

/** @brief Sole PW propagation entry point, independent of the selected method. */
template <typename T, typename Device>
class TDPropPW
{
  private:
    const ModulePW::PW_Basis_K& basis_;
    std::unique_ptr<hsolver::TDOrbOps<T, Device>> orb_ops_;
    std::unique_ptr<PWProp<T, Device>> prop_;

  public:
    /** @brief Take ownership of the initialized orbital service and propagator. */
    TDPropPW(const ModulePW::PW_Basis_K& basis,
             std::unique_ptr<hsolver::TDOrbOps<T, Device>> orb_ops,
             std::unique_ptr<PWProp<T, Device>> prop);
    PropStage stage(int scf_iter) const;
    void begin_trial();
    void finish_trial();
    void discard_trial();
    void invalidate_basis();
    hsolver::TDOrbOps<T, Device>& orbitals()
    {
        return *orb_ops_;
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

/** @brief Resolve input controls and construct only the selected propagation method. */
template <typename T, typename Device>
std::unique_ptr<TDPropPW<T, Device>> make_td_prop(const ModulePW::PW_Basis_K& basis,
                                                  const Input_para& input,
                                                  int npol,
                                                  const std::vector<int>& global_k,
                                                  const hsolver::diag_comm_info& comm,
                                                  std::ostream& log);
} // namespace pw
#endif
