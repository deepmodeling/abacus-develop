#ifndef ABACUS_LR_NACV_RESPONSE_H
#define ABACUS_LR_NACV_RESPONSE_H

#include "nacv_cg.h"
#include "zeq_solver.h"

namespace LR
{
/** @brief INPUT-validated CG controls; tolerance is the absolute response residual in Ry. */
struct NacvResponseControls
{
    double tolerance;
    int max_iterations;
};

/**
 * @brief Dispatch a NAC orbital-response solve using the existing distributed operators.
 * @param[out] Z State-major local solution rows.
 * @param[in] R Corresponding local right-hand sides in Ry.
 * @param[in] ld Local entries per row, including zero on empty ranks.
 * @param[in] nstates Number of right-hand sides.
 * @param[in,out] ops_L Existing A+B operator in the complete response window.
 * @param[in] nspin_x One for restricted, two for unrestricted response.
 * @param[in] zvec_solver Direct solver or NAC-certified CG.
 * @param[in] test_force Existing diagnostic-output selection.
 * @param[in] controls Validated NAC convergence settings.
 */
template <typename T, typename THamL>
inline void solve_nac_zeq_with(
    T *const Z,
    T *const R,
    const int ld,
    const int nstates,
    THamL &ops_L,
    const int nspin_x,
    const std::string &zvec_solver,
    const bool test_force,
    const NacvResponseControls &controls)
{
    for (int i = 0; i < nstates * ld; ++i)
    {
        Z[i] = T(0.0);
    } // clear Z
    if (zvec_solver == "cg")
    {
        const auto apply_response = [&ops_L, ld, nstates](const T *const in, T *const out)
        {
            ops_L.hPsi(
                in,
                out,
                ld,
                nstates);
        };
        solve_nac_cg(
            Z,
            R,
            ld,
            nstates,
            apply_response,
            test_force,
            controls.tolerance,
            controls.max_iterations);
    }
    else if (zvec_solver == "lapack")
    {
        solve_Z_lapack(
            Z,
            R,
            ld,
            nstates,
            ops_L,
            nspin_x,
            test_force);
    }
    else if (zvec_solver == "scalapack")
    {
        solve_Z_scalapack(
            Z,
            R,
            ld,
            nstates,
            ops_L,
            nspin_x);
    }
    else if (zvec_solver == "scalapack_chol")
    {
        solve_Z_scalapack_chol(
            Z,
            R,
            ld,
            nstates,
            ops_L,
            nspin_x);
    }
    else
    {
        throw std::runtime_error("Unsupported Z-vector solver: " + zvec_solver);
    }
}

/**
 * @brief Build the NAC polarized Z-vector response with explicit convergence controls.
 * @param[in] inputs Borrowed geometry, full orbital window, kernels and distributions.
 * @param[in] X State-major local polarized amplitudes.
 * @param[out] Z Response rows in the same layout.
 * @param[in] nstates Number of rows.
 * @param[in] pot LR response potential for the selected spin channel.
 * @param[in] spin_type Existing adaptation label.
 * @param[in] in_dir Existing input directory used by the response operators.
 * @param[in] openshell Whether amplitudes concatenate alpha and beta blocks.
 * @param[in] zvec_solver NAC response solver.
 * @param[in] controls Validated residual tolerance and update limit.
 */
template <typename T>
void solve_nac_zvector(
    const GradientInputs<T> &inputs,
    const T *X,
    T *Z,
    const int &nstates,
    std::weak_ptr<PotHxcLR> pot,
    const std::string &spin_type,
    const std::string &in_dir,
    bool openshell,
    const std::string &zvec_solver,
    const NacvResponseControls &controls);
} // namespace LR

#endif
