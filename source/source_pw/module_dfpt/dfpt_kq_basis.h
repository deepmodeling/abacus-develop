// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_KQ_BASIS_H
#define DFPT_KQ_BASIS_H

#include "source_base/vector3.h"
#include <vector>

namespace ModulePW {
class PW_Basis_K;
}

namespace ModuleDFPT {

/**
 * @brief Plane-wave basis at the perturbation wavevector k+q.
 *
 * C0: for every (ik, q) pair the first-order response (the Sternheimer
 * solution dpsi) lives in the k+q plane-wave basis. Rather than building a
 * full PW_Basis_K (new FFT grids, MP redistribution) per q, this helper
 * re-filters the G vectors of an already-initialized ground-state k-basis
 * PW_Basis_K at the shifted center k+q. Each G with |G + (k+q)|^2 <= gk_ecut
 * satisfies |G| <= sqrt(gk_ecut) + |k+q|, which is within the FFT-grid ball
 * (ggecut) the ground-state basis already distributed, so no G vector needed
 * by k+q is missing and only the shifted-center selection is performed.
 *
 * Preconditions:
 * - The ground-state basis must be a complex (gamma_only=false) k-basis:
 *   DFPT couples k and k+q symmetrically and the q-perturbation breaks the
 *   gamma-only half-space reduction.
 * - Every G vector needed by the largest k+q ball must lie inside the FFT
 *   grid of the ground-state basis, i.e. the FFT-grid cutoff (gridecut_lat)
 *   must satisfy gridecut_lat >= (sqrt(gk_ecut) + max|k| + max|q|)^2.
 *   In practice ecutrho >= 4*ecutwfc covers every q inside the first
 *   Brillouin zone, which is the DFPT q range.
 */
class DFPT_KQ_Basis {
public:
    DFPT_KQ_Basis();
    ~DFPT_KQ_Basis();

    /**
     * @brief Enumerate the local (per-processor) k+q plane-wave basis.
     * @param pw_wfc ground-state k-dependent plane-wave basis (complex)
     * @param q_cart perturbation wavevector in Cartesian coordinates
     * @param ik index of the ground-state k point
     */
    void init(const ModulePW::PW_Basis_K* pw_wfc,
              const ModuleBase::Vector3<double>& q_cart,
              int ik);

    void clear();

    bool is_valid() const { return pw_wfc_ != nullptr; }
    ///< number of k+q plane waves on this processor
    int get_npwk() const { return npwk_; }
    ///< index of the underlying ground-state G vector
    int get_ig(int igl) const { return igl2ig_[igl]; }
    ///< slab index (ig2isz) of the underlying ground-state G vector
    int get_ig2isz(int igl) const;
    ///< G in Cartesian coordinates
    ModuleBase::Vector3<double> get_gcar(int igl) const { return gcar_[igl]; }
    ///< G + (k+q) in Cartesian coordinates
    ModuleBase::Vector3<double> get_gpluskq(int igl) const { return gcar_[igl] + kplusq_c_; }
    ///< |G + (k+q)|^2 in units of 1/lat0^2
    double get_gk2(int igl) const { return gk2_[igl]; }
    ///< k+q wavevector in Cartesian coordinates
    ModuleBase::Vector3<double> get_kplusq() const { return kplusq_c_; }
    const std::vector<int>& get_igl2ig() const { return igl2ig_; }
    const std::vector<double>& get_gk2_all() const { return gk2_; }
    const std::vector<ModuleBase::Vector3<double>>& get_gcar_all() const { return gcar_; }

private:
    const ModulePW::PW_Basis_K* pw_wfc_ = nullptr; ///< ground-state k-basis
    ModuleBase::Vector3<double> kplusq_c_;         ///< k+q in Cartesian coordinates
    int npwk_ = 0;                                 ///< number of k+q plane waves
    std::vector<int> igl2ig_;                       ///< local k+q index -> base G index
    std::vector<double> gk2_;                       ///< |G + (k+q)|^2
    std::vector<ModuleBase::Vector3<double>> gcar_; ///< G in Cartesian coordinates
};

} // namespace ModuleDFPT

#endif // DFPT_KQ_BASIS_H