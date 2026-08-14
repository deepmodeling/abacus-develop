/**
 * @file dfpt_irrep_data.h
 * @brief Interface-adaptation layer exposing irrep-indexed DFPT data access.
 * @author Mohan Chen (added on 2026-05-18)
 * @note Phase 4 (DFPT wiring) interim layer: it simulates the option-2
 *       signatures (an irrep dimension added to the DFPT_PW_Data storage)
 *       on top of the current per-q storage, so that DFPT_PW::run can drive
 *       a per-irrep SCF loop and be verified first. Once verified, this
 *       adapter is sunk into DFPT_PW_Data as its official data layer.
 */
#ifndef DFPT_IRREP_DATA_H
#define DFPT_IRREP_DATA_H

#include "dfpt_pw_data.h"
#include <complex>
#include <map>
#include <utility>
#include <vector>

namespace ModuleDFPT {

/**
 * @brief Wrapper exposing irrep-indexed DFPT data access (option-2 signature).
 *
 * The underlying DFPT_PW_Data storage is indexed per q-point only; this
 * wrapper adds the little-group irrep dimension on top, delegating the
 * first-order quantities to the per-q storage and keeping a per-(q, irrep)
 * convergence state that the official data layer will absorb later.
 */
class DFPT_IrrepData {
public:
    /**
     * @brief Construct the adapter over a DFPT_PW_Data.
     * @param data underlying per-q data
     */
    explicit DFPT_IrrepData(DFPT_PW_Data& data);

    /**
     * @brief Number of q-points.
     */
    int get_nq() const;

    /**
     * @brief Number of irreps at given q-point.
     * @param q_idx q-point index
     */
    int get_nirr(int q_idx) const;

    /**
     * @brief Representative modes of a given irrep.
     * @param q_idx q-point index
     * @param irrep irrep index
     */
    std::vector<int> get_irrep_modes(int q_idx, int irrep) const;

    /**
     * @brief Irrep-indexed first-order wave-function coefficients.
     */
    void set_dpsi(int q_idx, int irrep, int k_idx, int band_idx,
                  const std::vector<std::complex<double>>& psi);
    std::vector<std::complex<double>> get_dpsi(int q_idx, int irrep, int k_idx, int band_idx) const;

    /**
     * @brief Irrep-indexed first-order charge density (real space).
     */
    void set_drho_r(int q_idx, int irrep, int spin, const std::vector<double>& rho);
    std::vector<double> get_drho_r(int q_idx, int irrep, int spin) const;

    /**
     * @brief Irrep-indexed first-order charge density (G space).
     */
    void set_drho_g(int q_idx, int irrep, int spin, const std::vector<std::complex<double>>& rho);
    std::vector<std::complex<double>> get_drho_g(int q_idx, int irrep, int spin) const;

    /**
     * @brief Irrep-indexed first-order potential (real space).
     */
    void set_dv_r(int q_idx, int irrep, int spin, const std::vector<double>& v);
    std::vector<double> get_dv_r(int q_idx, int irrep, int spin) const;

    /**
     * @brief Per-(q, irrep) SCF convergence bookkeeping.
     */
    void set_converged(int q_idx, int irrep, bool flag);
    bool get_converged(int q_idx, int irrep) const;
    void add_residual(int q_idx, int irrep, double r);
    std::vector<double> get_residuals(int q_idx, int irrep) const;
    void set_current_iter(int q_idx, int irrep, int iter);
    int get_current_iter(int q_idx, int irrep) const;

private:
    DFPT_PW_Data& data_;

    // per-(q, irrep) SCF state; the official data layer will store these
    // keyed by irrep once the adapter is sunk into DFPT_PW_Data.
    std::map<std::pair<int, int>, bool> converged_;
    std::map<std::pair<int, int>, std::vector<double>> residuals_;
    std::map<std::pair<int, int>, int> current_iter_;
};

} // namespace ModuleDFPT

#endif // DFPT_IRREP_DATA_H
