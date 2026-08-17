// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_HAMILT_SHIFT_H
#define DFPT_HAMILT_SHIFT_H

#include "dfpt_kq_basis.h"
#include "dfpt_stern.h"
#include "source_base/vector3.h"
#include <complex>
#include <vector>

namespace ModulePW {
class PW_Basis;
class PW_Basis_K;
}

class UnitCell;

namespace ModuleDFPT {

class DFPT_Pert;

/**
 * @brief Production adapter of the shifted Sternheimer operator (C7).
 *
 * Applies y = (H(k+q) - eps_n) x on the k+q plane-wave basis, with the
 * ground-state Hamiltonian assembled from module_dfpt primitives instead
 * of the GS operator chain (which is index-bound to the GS k list):
 *   - kinetic        tpiba^2 |G+k+q|^2 (diagonal, DFPT_KQ_Basis)
 *   - local V_eff    FFT convolution with the real-space effective
 *                    potential injected by the esolver after the GS SCF
 *                    (veff_smooth row, shared FFT grid)
 *   - nonlocal V_nl  separable projectors from DFPT_Pert::build_vkb at
 *                    the k+q center, cached per (k,q) context
 *   - shift          -eps_n fixed before each solve
 * All transforms are phase-free on the shared FFT grid (the same
 * convention as DFPT_Pert::real_space_dv), and the k+q scatter/gather
 * goes through the (ix,iy,iz) FFT-cell reverse map (C1 finding: the rho
 * and wfc stick encodings are not interchangeable).
 */
class DFPT_HamiltShift : public DFPT_Stern::LinearOperator {
public:
    DFPT_HamiltShift(const UnitCell& ucell,
                     ModulePW::PW_Basis* pw_rho,
                     ModulePW::PW_Basis_K* pw_wfc,
                     const std::vector<double>& veff_r,
                     const DFPT_Pert* pert);
    ~DFPT_HamiltShift();

    DFPT_HamiltShift(const DFPT_HamiltShift&) = delete;
    DFPT_HamiltShift& operator=(const DFPT_HamiltShift&) = delete;

    /// Fix the (k, q) context and cache the k+q projector set; the
    /// eigenvalue shift (Ry) is set per solve through set_shift.
    void set_context(const ModuleBase::Vector3<double>& q_cart, int k_idx);
    void set_shift(double shift);

    int dimension() const override;
    void apply(const std::complex<double>* x, std::complex<double>* y) const override;

private:
    const UnitCell* ucell_ = nullptr;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
    const DFPT_Pert* pert_ = nullptr;
    std::vector<double> veff_r_;
    double tpiba2_ = 0.0;
    int nrxx_ = 0;

    DFPT_KQ_Basis kq_;
    double shift_ = 0.0;
    ///< k+q list index -> rho-grid ig (-1 when the cell position carries no
    /// rho G; cannot happen with ecutrho >= 4 ecutwfc but kept defensive)
    std::vector<int> kq2rho_;

    ///< projector bookkeeping per type (mu -> (beta index, m channel))
    std::vector<std::vector<int>> mu_ib_;
    std::vector<std::vector<int>> mu_m_;
    ///< cached beta projectors per atom on the k+q list
    std::vector<std::vector<std::vector<std::complex<double>>>> vkb_cache_;

    ///< apply-time scratch (apply is const)
    mutable std::vector<std::complex<double>> x_recip_;
    mutable std::vector<std::complex<double>> x_r_;
    mutable std::vector<std::complex<double>> becp_;
    mutable std::vector<std::complex<double>> dbecp_;
};

} // namespace ModuleDFPT

#endif // DFPT_HAMILT_SHIFT_H
