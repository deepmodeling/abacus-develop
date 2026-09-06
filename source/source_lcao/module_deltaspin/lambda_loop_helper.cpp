#include "lambda_loop_helper.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>

#include "basic_funcs.h"
#include "source_base/formatter.h"

/**
 * @file lambda_loop_helper.cpp
 * @brief Free-function implementations of lambda loop helpers.
 *
 * @par History
 * Originally these were explicit specializations of SpinConstrain member
 * functions:
 *   template <> void SpinConstrain<std::complex<double>>::print_termination(...)
 * They have been lifted to free function templates so that the SpinConstrain
 * class no longer has to carry the lambda-loop workflow as part of its
 * interface. The complex-double specialization is provided by instantiating
 * the template with TK = std::complex<double>; a double stub is provided
 * elsewhere to keep the link surface stable.
 */

namespace spinconstrain
{

/**
 * @brief Print final spin and lambda values when the lambda loop terminates.
 *
 * @par Output
 * - "after-optimization spin (uB)": Final magnetic moments Mi per atom
 * - "after-optimization lambda (eV/uB)": Final Lagrange multipliers per atom
 * - "Inner optimization for lambda ends.": Termination marker
 */
template <typename TK>
void print_termination(const SpinConstrain<TK>& sc, std::ostream& ofs_running)
{
    print_2d(" after-optimization spin (uB): (print in the inner loop): ", sc.get_Mi(), sc.get_nspin(), 1.0, ofs_running);
    print_2d(" after-optimization lambda (eV/uB): (print in the inner loop): ",
             sc.get_sc_lambda(),
             sc.get_nspin(),
             ModuleBase::Ry_to_eV,
             ofs_running);
    ofs_running << " Inner optimization for lambda ends." << std::endl;
    ofs_running << " ================================================================================" << std::endl;
}

/**
 * @brief Check if RMS error is below convergence threshold or max steps reached.
 *
 * @par Output
 * Prints step info: "Step (Outer -- Inner) = X -- Y   RMS = Z   TIME(s) = T"
 *
 * @par Termination messages
 * - "Meet convergence criterion": RMS < current_sc_thr_ (successfully converged)
 * - "Reach maximum number of steps": i_step == nsc_ - 1 (did not converge)
 *
 * @par Return value
 * - true: loop should terminate (either converged or max steps)
 * - false: continue optimization
 */
template <typename TK>
bool check_rms_stop(const SpinConstrain<TK>& sc,
                    int outer_step,
                    int i_step,
                    double rms_error,
                    double duration,
                    double total_duration,
                    std::ostream& ofs_running)
{
    ofs_running << " Step (Outer -- Inner) =  " << outer_step << " -- " << std::left << std::setw(5) << i_step + 1
                         << "       RMS = " << rms_error << "     TIME(s) = " << std::setw(11) << duration << std::endl;
    const double current_sc_thr = sc.get_current_sc_thr();
    const int nsc = sc.get_nsc();
    if (rms_error < current_sc_thr || i_step == nsc - 1)
    {
        if (rms_error < current_sc_thr)
        {
            ofs_running << " DeltaSpin: lambda loop converged ( RMS < " << current_sc_thr
                                 << " ), inner steps = " << (i_step + 1)
                                 << ", Total TIME(s) = " << total_duration << std::endl;
            ofs_running << std::endl;
        }
        else if (i_step == nsc - 1)
        {
            std::cout << " DeltaSpin: lambda loop reached max steps ( " << nsc
                      << " ), RMS = " << rms_error
                      << ", Total TIME(s) = " << total_duration << std::endl;
            std::cout << std::endl;
        }
        print_termination(sc, ofs_running);
        return true;
    }
    return false;
}

/// @brief Print header at start of lambda optimization loop
template <typename TK>
void print_header(const SpinConstrain<TK>& sc, std::ostream& ofs_running)
{
    ofs_running << " ================================================================================" << std::endl;
    ofs_running << " Inner optimization for lambda begins ..." << std::endl;
    ofs_running << " Covergence criterion for the iteration: " << sc.get_sc_thr() << std::endl;
}

/**
 * @brief Cap step size to prevent overshooting in lambda optimization.
 *
 * @details If |alpha_trial * max(search)| > restrict_current_, reduce alpha_trial
 * so that the maximum lambda change per step is bounded by restrict_current_.
 *
 * This prevents the optimizer from taking steps that are too large, which
 * could lead to oscillation or divergence.
 *
 * @par Output (when restriction is applied)
 * - "alpha after restrict = X eV/uB^2": The capped step size
 * - "boundary after = X eV/uB": The actual maximum lambda change
 */
template <typename TK>
void check_restriction(const SpinConstrain<TK>& sc,
                       const std::vector<ModuleBase::Vector3<double>>& search,
                       double& alpha_trial,
                       std::ostream& ofs_running)
{
    const double restrict_current = sc.get_sccut();
    double boundary = std::abs(alpha_trial * maxval_abs_2d(search));

    if (restrict_current > 0 && boundary > restrict_current)
    {
        alpha_trial = copysign(1.0, alpha_trial) * restrict_current / maxval_abs_2d(search);
        boundary = std::abs(alpha_trial * maxval_abs_2d(search));
        ofs_running << " alpha after restrict = " << alpha_trial * ModuleBase::Ry_to_eV << std::endl;
        ofs_running << " boundary after = " << boundary * ModuleBase::Ry_to_eV << std::endl;
    }
}

/**
 * @brief Compute optimal step size via linear interpolation.
 *
 * @par Algorithm
 * Uses the two-point linear interpolation (secant method) to find the
 * step size that would drive Mi to M_target:
 *
 *   alpha_opt = sum_k / sum_k2 * alpha_trial
 *
 * where:
 *   sum_k  = sum((target - spin) . (spin_plus - spin))   over constrained components
 *   sum_k2 = sum(|spin - spin_plus|^2)                   over constrained components
 *
 * @par Edge case handling
 * - If |sum_k2| < 1e-30: spin and spin_plus are nearly identical, meaning
 *   the lambda change has no effect on Mi. Return alpha_trial as fallback.
 */
template <typename TK>
double cal_alpha_opt(const SpinConstrain<TK>& sc,
                     std::vector<ModuleBase::Vector3<double>> spin,
                     std::vector<ModuleBase::Vector3<double>> spin_plus,
                     const double alpha_trial)
{
    int nat = sc.get_nat();
    const bool print = false;
    const double zero = 0.0;

    // Mask to only constrained components
    std::vector<ModuleBase::Vector3<double>> spin_mask(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> target_spin_mask(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> spin_plus_mask(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> temp_1(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> temp_2(nat, 0.0);
    where_fill_scalar_else_2d(sc.get_constrain(), 0, zero, sc.get_target_mag(), target_spin_mask);
    where_fill_scalar_else_2d(sc.get_constrain(), 0, zero, spin, spin_mask);
    where_fill_scalar_else_2d(sc.get_constrain(), 0, zero, spin_plus, spin_plus_mask);

    // Compute dot products for linear interpolation
    for (int ia = 0; ia < nat; ia++)
    {
        for (int ic = 0; ic < 3; ic++)
        {
            // sum_k: (target - current) . (trial - current)
            temp_1[ia][ic]
                = (target_spin_mask[ia][ic] - spin_mask[ia][ic]) * (spin_plus_mask[ia][ic] - spin_mask[ia][ic]);
            // sum_k2: |current - trial|^2
            temp_2[ia][ic] = std::pow(spin_mask[ia][ic] - spin_plus_mask[ia][ic], 2);
        }
    }
    double sum_k = sum_2d(temp_1);
    double sum_k2 = sum_2d(temp_2);

    // Debug output (controlled by print flag)
    for(int ia=0; ia<std::min(2,(int)nat); ++ia) {
        if (print) {
        printf("[ALPHA-OPT] nat=%d sum_k=%.6e sum_k2=%.6e alpha_trial=%.6e\n", nat, sum_k, sum_k2, alpha_trial);
        printf("[ALPHA-OPT] spin[%d]=(%.4f,%.4f,%.4f) spin_plus[%d]=(%.4f,%.4f,%.4f)\n",
                ia, spin[ia].x, spin[ia].y, spin[ia].z,
                ia, spin_plus[ia].x, spin_plus[ia].y, spin_plus[ia].z);
        }
    }

    // Guard against division by zero
    if (std::abs(sum_k2) < 1e-30) {
        if (print) {
        printf("[ALPHA-OPT] WARNING: sum_k2 too small, returning alpha_trial\n");
        }
        fflush(stdout);
        return alpha_trial;
    }
    fflush(stdout);
    return sum_k * alpha_trial / sum_k2;
}

/**
 * @brief Check if the magnetic susceptibility gradient dM/dlambda has decayed below threshold.
 *
 * @par Algorithm
 * 1. Compute spin_change = new_spin - spin (change in magnetic moments)
 * 2. Compute nu_change = delta_lambda - dnu_last_step (change in lambda)
 * 3. Compute full gradient matrix: dM[ia][ic]/dlambda[ja][jc] = spin_change[ia][ic] / nu_change[ja][jc]
 * 4. Extract diagonal: dM[ia][ic]/dlambda[ia][ic] (self-susceptibility)
 * 5. Find max diagonal gradient per atom type
 * 6. If max_gradient[itype] < decay_grad[itype], return true (early termination)
 *
 * @par Physical meaning
 * The diagonal gradient dM/dlambda represents how sensitive the magnetic moment
 * is to changes in the Lagrange multiplier. When this gradient becomes very small,
 * further increases in lambda produce diminishing returns in Mi, indicating that
 * the optimization has reached its practical limit.
 *
 * @par Output (when triggered)
 * "Reach limitation of current step ( maximum gradient < X uB^2/eV in atom type Y ), exit."
 */
template <typename TK>
bool check_gradient_decay(const SpinConstrain<TK>& sc,
                          std::vector<ModuleBase::Vector3<double>> new_spin,
                          std::vector<ModuleBase::Vector3<double>> spin,
                          std::vector<ModuleBase::Vector3<double>> delta_lambda,
                          std::vector<ModuleBase::Vector3<double>> dnu_last_step,
                          bool print,
                          std::ostream& ofs_running)
{
    const double one = 1.0;
    const double zero = 0.0;
    int nat = sc.get_nat();
    int ntype = sc.get_ntype();

    // Change in magnetic moments and lambda
    std::vector<ModuleBase::Vector3<double>> spin_change(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> nu_change(nat, 1.0);

    // Full gradient matrix: dM[ia][ic]/dlambda[ja][jc]
    std::vector<std::vector<std::vector<std::vector<double>>>> spin_nu_gradient(
        nat,
        std::vector<std::vector<std::vector<double>>>(
            3,
            std::vector<std::vector<double>>(nat, std::vector<double>(3, 0.0))));
    // Diagonal gradient: dM[ia][ic]/dlambda[ia][ic] (self-susceptibility)
    std::vector<ModuleBase::Vector3<double>> spin_nu_gradient_diag(nat, 0.0);
    std::vector<std::pair<int, int>> max_gradient_index(ntype, std::make_pair(0, 0));
    std::vector<double> max_gradient(ntype, 0.0);

    subtract_2d(new_spin, spin, spin_change);
    subtract_2d(delta_lambda, dnu_last_step, nu_change);

    const auto& constrain = sc.get_constrain();
    // Mask unconstrained components
    where_fill_scalar_2d(constrain, 0, zero, spin_change);
    where_fill_scalar_2d(constrain, 0, one, nu_change);

    // Calculate full gradient matrix
    for (int ia = 0; ia < nat; ia++)
    {
        for (int ic = 0; ic < 3; ic++)
        {
            for (int ja = 0; ja < nat; ja++)
            {
                for (int jc = 0; jc < 3; jc++)
                {
                    if (std::abs(nu_change[ja][jc]) < 1e-30) {
                        printf("[GRAD-DECAY] WARNING: nu_change[%d][%d] too small! delta_lambda=(%.6e,%.6e,%.6e) dnu_last=(%.6e,%.6e,%.6e)\n",
                               ja, jc, delta_lambda[ja].x, delta_lambda[ja].y, delta_lambda[ja].z,
                               dnu_last_step[ja].x, dnu_last_step[ja].y, dnu_last_step[ja].z);
                        fflush(stdout);
                        nu_change[ja][jc] = 1e-30;
                    }
                    spin_nu_gradient[ia][ic][ja][jc] = spin_change[ia][ic] / nu_change[ja][jc];
                }
            }
        }
    }

    const auto& atom_counts = sc.get_atomCounts();
    const auto& decay_grad = sc.get_decay_grad();
    const int nspin = sc.get_nspin();

    // Extract diagonal gradient and find max per atom type
    for (const auto& sc_elem: atom_counts)
    {
        int it = sc_elem.first;
        int nat_it = sc_elem.second;
        max_gradient[it] = 0.0;
        for (int ia = 0; ia < nat_it; ia++)
        {
            for (int ic = 0; ic < 3; ic++)
            {
                spin_nu_gradient_diag[ia][ic] = spin_nu_gradient[ia][ic][ia][ic];
                if (std::abs(spin_nu_gradient_diag[ia][ic]) > std::abs(max_gradient[it]))
                {
                    max_gradient[it] = spin_nu_gradient_diag[ia][ic];
                    max_gradient_index[it].first = ia;
                    max_gradient_index[it].second = ic;
                }
            }
        }
    }

    if (print)
    {
        print_2d(" diagonal gradient: ", spin_nu_gradient_diag, nspin, 1.0, ofs_running);
        ofs_running << " maximum gradient appears at: " << std::endl;
        for (int it = 0; it < ntype; it++)
        {
            ofs_running << " ( " << max_gradient_index[it].first << ", " << max_gradient_index[it].second << " )"
                                 << std::endl;
        }
        ofs_running << " maximum gradient: " << std::endl;
        for (int it = 0; it < ntype; it++)
        {
            ofs_running << " " << max_gradient[it]/ModuleBase::Ry_to_eV << std::endl;
        }
    }

    // Check if any atom type's gradient has decayed below threshold
    for (int it = 0; it < ntype; it++)
    {
        if (decay_grad[it] > 0 && std::abs(max_gradient[it]) < decay_grad[it])
        {
            std::cout << " DeltaSpin: lambda loop early-terminated ( maximum gradient < "
                      << decay_grad[it]/ModuleBase::Ry_to_eV // uB^2/Ry to uB^2/eV
                      << " in atom type " << it << " )" << std::endl;
            std::cout << std::endl;
            return true;
        }
    }
    return false;
}

/**
 * @brief Print atomic magnetic moments Mi in a formatted table.
 *
 * @par Output format
 * - nspin=2: "Total Magnetism (uB)" with single z-component column
 * - nspin=4: three columns (Mx, My, Mz)
 *
 * Lifted from SpinConstrain<TK>::print_Mi; accesses state via getters.
 */
template <typename TK>
void print_Mi(const SpinConstrain<TK>& sc, std::ostream& ofs_running)
{
    sc.check_atomCounts();
    const int nat = sc.get_nat();
    const int nspin = sc.get_nspin();
    const auto& Mi = sc.get_Mi();
    const auto& atomLabel = sc.get_atomLabels();
    std::vector<double> mag_x(nat, 0.0);
    std::vector<double> mag_y(nat, 0.0);
    std::vector<double> mag_z(nat, 0.0);
    if (nspin == 2)
    {
        const std::vector<std::string> title = {"Total Magnetism (uB)", ""};
        const std::vector<std::string> fmts = {"%-26s", "%20.10f"};
        FmtTable table(/*titles=*/title,
                       /*nrows=*/nat,
                       /*formats=*/fmts,
                       /*indent=*/0,
                       /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::LEFT});
        for (int iat = 0; iat < nat; ++iat)
        {
            mag_z[iat] = Mi[iat].z;
        }
        table << atomLabel << mag_z;
        ofs_running << table.str() << std::endl;
    }
    else if (nspin == 4)
    {
        const std::vector<std::string> title = {"Total Magnetism (uB)", "", "", ""};
        const std::vector<std::string> fmts = {"%-26s", "%20.10f", "%20.10f", "%20.10f"};
        FmtTable table(/*titles=*/title,
                       /*nrows=*/nat,
                       /*formats=*/fmts,
                       /*indent=*/0,
                       /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::LEFT});
        for (int iat = 0; iat < nat; ++iat)
        {
            mag_x[iat] = Mi[iat].x;
            mag_y[iat] = Mi[iat].y;
            mag_z[iat] = Mi[iat].z;
        }
        table << atomLabel << mag_x << mag_y << mag_z;
        ofs_running << table.str() << std::endl;
    }
}

/**
 * @brief Print the magnetic force (-lambda) per atom in eV/uB.
 *
 * Lifted from SpinConstrain<TK>::print_Mag_Force; accesses state via getters.
 * lambda is read via get_sc_lambda() and converted from Ry to eV on output.
 */
template <typename TK>
void print_Mag_Force(const SpinConstrain<TK>& sc, std::ostream& ofs_running)
{
    sc.check_atomCounts();
    const int nat = sc.get_nat();
    const int nspin = sc.get_nspin();
    const auto& lambda = sc.get_sc_lambda();
    const auto& atomLabel = sc.get_atomLabels();
    std::vector<double> mag_force_x(nat, 0.0);
    std::vector<double> mag_force_y(nat, 0.0);
    std::vector<double> mag_force_z(nat, 0.0);
    if (nspin == 2)
    {
        const std::vector<std::string> title = {"Magnetic force (eV/uB)", ""};
        const std::vector<std::string> fmts = {"%-26s", "%20.10f"};
        FmtTable table(/*titles=*/title,
                       /*nrows=*/nat,
                       /*formats=*/fmts,
                       /*indent=*/0,
                       /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::LEFT});
        for (int iat = 0; iat < nat; ++iat)
        {
            mag_force_z[iat] = lambda[iat].z * ModuleBase::Ry_to_eV;
        }
        table << atomLabel << mag_force_z;
        ofs_running << table.str() << std::endl;
    }
    else if (nspin == 4)
    {
        const std::vector<std::string> title = {"Magnetic force (eV/uB)", "", "", ""};
        const std::vector<std::string> fmts = {"%-26s", "%20.10f", "%20.10f", "%20.10f"};
        FmtTable table(/*titles=*/title,
                       /*nrows=*/nat,
                       /*formats=*/fmts,
                       /*indent=*/0,
                       /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::LEFT});
        for (int iat = 0; iat < nat; ++iat)
        {
            mag_force_x[iat] = lambda[iat].x * ModuleBase::Ry_to_eV;
            mag_force_y[iat] = lambda[iat].y * ModuleBase::Ry_to_eV;
            mag_force_z[iat] = lambda[iat].z * ModuleBase::Ry_to_eV;
        }
        table << atomLabel << mag_force_x << mag_force_y << mag_force_z;
        ofs_running << table.str() << std::endl;
    }
}


// Explicit instantiation for the only supported TK = std::complex<double>.
// The double stub is provided by template_helpers.cpp via the existing
// specialization mechanism (kept as a separate file to avoid duplicate symbols).
template void print_termination<std::complex<double>>(const SpinConstrain<std::complex<double>>&, std::ostream&);
template bool check_rms_stop<std::complex<double>>(const SpinConstrain<std::complex<double>>&,
                                                   int, int, double, double, double, std::ostream&);
template void print_header<std::complex<double>>(const SpinConstrain<std::complex<double>>&, std::ostream&);
template void check_restriction<std::complex<double>>(const SpinConstrain<std::complex<double>>&,
                                                       const std::vector<ModuleBase::Vector3<double>>&,
                                                       double&, std::ostream&);
template double cal_alpha_opt<std::complex<double>>(const SpinConstrain<std::complex<double>>&,
                                                    std::vector<ModuleBase::Vector3<double>>,
                                                    std::vector<ModuleBase::Vector3<double>>,
                                                    const double);
template bool check_gradient_decay<std::complex<double>>(const SpinConstrain<std::complex<double>>&,
                                                          std::vector<ModuleBase::Vector3<double>>,
                                                          std::vector<ModuleBase::Vector3<double>>,
                                                          std::vector<ModuleBase::Vector3<double>>,
                                                          std::vector<ModuleBase::Vector3<double>>,
                                                          bool, std::ostream&);

// print_Mi / print_Mag_Force are generic (no per-TK stub needed): the
// template body works for both TK = std::complex<double> (real usage) and
// TK = double (nspin=2 stub path instantiated by template_helpers_test).
// We instantiate both so callers in PW/LCAO ESolvers that hold either TK
// can link against a single definition.
template void print_Mi<std::complex<double>>(const SpinConstrain<std::complex<double>>&, std::ostream&);
template void print_Mag_Force<std::complex<double>>(const SpinConstrain<std::complex<double>>&, std::ostream&);
template void print_Mi<double>(const SpinConstrain<double>&, std::ostream&);
template void print_Mag_Force<double>(const SpinConstrain<double>&, std::ostream&);

} // namespace spinconstrain
