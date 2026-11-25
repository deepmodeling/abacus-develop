#ifndef STO_DOS_QUALITY_H
#define STO_DOS_QUALITY_H

#include <string>
#include <vector>

/**
 * @file sto_dos_quality.h
 * @brief Quality guards for the SDFT/MDFT stochastic DOS pipeline.
 *
 * Pure, global-state free helpers, all unit tested in
 * `unittests/test_sto_dos_quality.cpp`:
 *
 *  1. `sanitize_che_range()` repairs the Chebyshev energy window [Emin, Emax]
 *     so that it covers the known deterministic eigenvalues and the physical
 *     plane-wave spectrum (`ecutwfc`).  Without it a positive lower bound
 *     excludes every occupied state and T_n(|x| > 1) grows exponentially
 *     (the historical 1e26--1e64 DOS values).  The interval is a property of
 *     the spectrum and is deliberately *not* widened to fit a requested DOS
 *     output grid: at a fixed order that would lower the resolution.
 *
 *  2. `check_dos_parameters()` / `check_dos_result()` report DOS settings that
 *     are silently useless at high temperature (grid spacing vs Gaussian
 *     broadening, Chebyshev order vs range/sigma) and a DOS curve that does
 *     not integrate to the number of plane-wave states, goes negative, or was
 *     evaluated outside the Chebyshev window.
 *
 * Both groups only report through their return values; the caller decides what
 * to print, so a long run is never aborted by a quality suggestion.
 */

namespace ModuleIO
{

/// A Chebyshev energy window plus the flags telling whether it was repaired.
struct CheRange
{
    double emin = 0.0;          ///< lower bound [Ry]
    double emax = 0.0;          ///< upper bound [Ry]
    bool emin_repaired = false; ///< whether the lower bound was modified
    bool emax_repaired = false; ///< whether the upper bound was modified
    bool degenerate = false;    ///< input was non-finite/inverted and was reset
};

/**
 * @brief Check and, if necessary, repair the Chebyshev spectral window.
 *
 * @param[in] emin           automatic/user lower bound [Ry]
 * @param[in] emax           automatic/user upper bound [Ry]
 * @param[in] ecutwfc        wavefunction cutoff [Ry]; when > 0 it is used as
 *                           the lower bound of Emax (`1.05 * ecutwfc`)
 * @param[in] det_eigs       known deterministic (Kohn-Sham) eigenvalues [Ry];
 *                           may be empty
 *
 * Rules, applied in this fixed order (kept simple so the behaviour is testable):
 *  1. non-finite or inverted input is reset to a conservative interval;
 *  2. deterministic eigenvalues are covered (margin 1e-3 Ry);
 *  3. a pathological positive Emin is replaced by `-max(5, 0.01*|emax|)`;
 *  4. `Emax >= 1.05 * ecutwfc` when `ecutwfc > 0`;
 *  5. a minimum width of 1e-6 Ry is enforced.
 *
 * Healthy input is returned bit for bit unchanged.
 */
CheRange sanitize_che_range(double emin,
                            double emax,
                            double ecutwfc,
                            const std::vector<double>& det_eigs);

/**
 * @brief Parameter-level DOS checks, called before the DOS loop.
 *
 * @param range_ev width of the Chebyshev interval [eV]
 * @param sigma_ev `dos_sigma` [eV] (may be zero or negative on bad input)
 * @param de_ev    `dos_edelta_ev` [eV]
 * @param dos_nche number of Chebyshev orders used for the DOS expansion
 * @return one human-readable warning per detected problem; empty when sane
 */
std::vector<std::string> check_dos_parameters(double range_ev,
                                              double sigma_ev,
                                              double de_ev,
                                              int dos_nche);

/// Result-level checks, called after the DOS curve has been reduced.
/// @param min_dos   smallest DOS value on the grid
/// @param max_abs_t max |t| of the output grid in the Chebyshev variable
std::vector<std::string> check_dos_result(double min_dos, double max_abs_t);

} // namespace ModuleIO

#endif // STO_DOS_QUALITY_H
