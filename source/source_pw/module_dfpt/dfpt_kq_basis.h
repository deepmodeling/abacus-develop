// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in design phase and has not been
// put into production yet.
// It may change in the future.
// Please use this code with caution.
// Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_KQ_BASIS_H
#define DFPT_KQ_BASIS_H

#include "source_base/vector3.h"
#include <vector>

namespace ModulePW {
class PW_Basis;
class PW_Basis_K;
}

namespace ModuleDFPT {

/**
 * @brief Plane-wave basis at the perturbation wavevector k+q.
 *
 * C0: for every (ik, q) pair the first-order response (the Sternheimer
 * solution dpsi) lives in the k+q plane-wave basis. Rather than building a
 * full PW_Basis_K (new FFT grids, MP redistribution) per q, this helper
 * enumerates the G vectors of two already-initialized ground-state bases at
 * the shifted center k+q: every G with |G + (k+q)|^2 <= gk_ecut is kept.
 *
 * Two reservoirs are needed because the wavefunction G grid is distributed
 * with the ball radius (sqrt(gk_ecut) + max|k_mesh|)^2, which does not cover
 * the q-shifted balls (the k-mesh knows nothing about q). The dense charge
 * grid, built with ecutrho >= 4*ecutwfc, has radius 2*sqrt(gk_ecut) and
 * covers sqrt(gk_ecut) + |k+q| for every q inside the first Brillouin
 * zone, so it completes the enumeration.
 *
 * Preconditions:
 * - The ground-state k-basis must be complex (gamma_only=false): DFPT
 *   couples k and k+q symmetrically and the q-perturbation breaks the
 *   gamma-only half-space reduction.
 * - Both bases must share the same FFT grid dimensions (the k+q G vectors
 *   are exchanged between them through the shared FFT cell position).
 */
class DFPT_KQ_Basis {
public:
    DFPT_KQ_Basis();
    ~DFPT_KQ_Basis();

    /**
     * @brief Enumerate the local (per-processor) k+q plane-wave basis.
     * @param pw_wfc ground-state k-dependent plane-wave basis (complex)
     * @param pw_rho ground-state charge-density plane-wave basis
     * @param q_cart perturbation wavevector in Cartesian coordinates
     * @param ik index of the ground-state k point
     */
    void init(const ModulePW::PW_Basis_K* pw_wfc,
              const ModulePW::PW_Basis* pw_rho,
              const ModuleBase::Vector3<double>& q_cart,
              int ik);

    void clear();

    bool is_valid() const { return pw_wfc_ != nullptr; }
    ///< number of k+q plane waves on this processor
    int get_npwk() const { return npwk_; }
    ///< index of the G vector in the charge-density basis (-1 if the shared
    ///< FFT cell position carries no local rho-grid G)
    int get_ig_rho(int igl) const { return ig_rho_[igl]; }
    ///< G in Cartesian coordinates
    ModuleBase::Vector3<double> get_gcar(int igl) const { return gcar_[igl]; }
    ///< G + (k+q) in Cartesian coordinates
    ModuleBase::Vector3<double> get_gpluskq(int igl) const { return gcar_[igl] + kplusq_c_; }
    ///< |G + (k+q)|^2 in units of 1/lat0^2
    double get_gk2(int igl) const { return gk2_[igl]; }
    ///< k+q wavevector in Cartesian coordinates
    ModuleBase::Vector3<double> get_kplusq() const { return kplusq_c_; }
    const std::vector<double>& get_gk2_all() const { return gk2_; }
    const std::vector<ModuleBase::Vector3<double>>& get_gcar_all() const { return gcar_; }

private:
    const ModulePW::PW_Basis_K* pw_wfc_ = nullptr; ///< ground-state k-basis
    ModuleBase::Vector3<double> kplusq_c_;         ///< k+q in Cartesian coordinates
    int npwk_ = 0;                                 ///< number of k+q plane waves
    std::vector<int> ig_rho_;                      ///< k+q index -> charge-grid G index
    std::vector<double> gk2_;                      ///< |G + (k+q)|^2
    std::vector<ModuleBase::Vector3<double>> gcar_;///< G in Cartesian coordinates
};

} // namespace ModuleDFPT

#endif // DFPT_KQ_BASIS_H
