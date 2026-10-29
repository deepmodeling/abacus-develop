#ifndef TD_SPIN_PW_H
#define TD_SPIN_PW_H

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_estate/elecstate.h"
#include "source_pw/module_proj/onsite_proj.h"
#include "source_pw/module_pwdft/nonlocal_workspace.h"
#include "source_pw/module_rt/td_obs_pw.h"

#include <complex>
#include <memory>
#include <string>
#include <vector>

namespace pw
{
/** @brief Periodic nearest-atom labels on the local charge grid, in global atom order. */
std::vector<int> td_voronoi_pw(const UnitCell& cell, const ModulePW::PW_Basis& grid);

/** @brief Endpoint magnetism with reusable device workspaces and position-aware caches. */
template <typename Real, typename Device>
class TDSpinPW
{
  private:
    using Complex = std::complex<Real>;
    const UnitCell* cell_ = nullptr;
    const ModulePW::PW_Basis* grid_ = nullptr;
    const ModulePW::PW_Basis_K* basis_ = nullptr;
    const K_Vectors* kv_ = nullptr;
    Structure_Factor* sf_ = nullptr;
    pseudopot_cell_vnl* pp_ = nullptr;
    std::string mag_partition_;
    std::string orbital_dir_;
    double onsite_radius_ = 0.0;
    bool calc_torque_ = false;
    std::vector<ModuleBase::Vector3<double>> pos_cache_;
    std::vector<int> labels_;
    hamilt::NonlocalWorkspace<Complex, Device> nonlocal_;
    ct::Tensor beta_;
    std::vector<Complex> projections_;
    std::vector<Complex> contracted_;
    std::unique_ptr<psi::Psi<std::complex<double>, Device>> onsite_psi_;
    // Destroy the projector before the wavefunction it borrows.
    std::unique_ptr<projectors::OnsiteProjector<double, Device>> onsite_projector_;
    double q_bound_ = 0.0;

    void sync_positions();
    void density(const elecstate::ElecState& state, TDSpinResultPW* result);
    void soc(const psi::Psi<Complex, Device>& psi,
             const ModuleBase::matrix& weights,
             const ModuleBase::Vector3<double>& a,
             TDSpinResultPW* result);
    void onsite(const psi::Psi<Complex, Device>& psi,
                const elecstate::ElecState& state,
                const ModuleBase::Vector3<double>& a,
                TDSpinResultPW* result);

  public:
    /** @brief Bind once; dependencies outlive this object and retain their cell and layout.
     * The caller must refresh structure factors and electronic data before sampling moved ions.
     * @param mag_partition Atomic moment partition: none, voronoi, or onsite.
     * @param orbital_dir Directory containing the onsite projector orbitals.
     * @param onsite_radius Onsite projector cutoff in Bohr.
     * @param calc_torque Whether to calculate XC and SOC spin torques.
     */
    void initialize(const UnitCell& cell,
                    const ModulePW::PW_Basis& grid,
                    const ModulePW::PW_Basis_K& basis,
                    const K_Vectors& kv,
                    Structure_Factor& sf,
                    pseudopot_cell_vnl& pp,
                    const std::string& mag_partition,
                    const std::string& orbital_dir,
                    double onsite_radius,
                    bool calc_torque);
    /** @brief Sample all ranks and restore the k cursor; coordinate changes invalidate geometric caches. */
    TDSpinResultPW calculate(const psi::Psi<Complex, Device>& psi, const elecstate::ElecState& state, const ModuleBase::Vector3<double>& a);
};
} // namespace pw
#endif
