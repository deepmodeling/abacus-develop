#include "sto_dos_quality.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace
{
// ---------------------------------------------------------------------------
// Constants of the Chebyshev-window sanitiser
// ---------------------------------------------------------------------------
const double det_margin = 1.0e-3; ///< safety margin below/above deterministic eigenvalues [Ry]
const double min_width = 1.0e-6;  ///< minimal window width [Ry], avoids division by zero
const double ecut_margin = 1.05;  ///< Emax floor relative to ecutwfc

/// Conservative fallback for a pathological lower bound.
double fallback_emin(double emax)
{
    return -std::max(5.0, 0.01 * std::fabs(emax));
}

// ---------------------------------------------------------------------------
// Constants of the DOS parameter / result checks
// ---------------------------------------------------------------------------
const double negative_dos_abs_tol = 1.0e-6; ///< |min DOS| below this is ringing noise
const double de_over_sigma_min = 0.02;        ///< grid much finer than the broadening (waste)
const double de_over_sigma_max = 3.0;         ///< grid coarser than the broadening (information loss)
const double order_prefactor = std::sqrt(2.0 * std::log(1.0e6)); ///< target 1e-6 in the Gaussian
const double max_abs_t_tol = 1.0e-9;          ///< tolerance on |t| <= 1

std::string format_double(const double value)
{
    std::ostringstream os;
    os << value;
    return os.str();
}

// ---------------------------------------------------------------------------
// Rules of the Chebyshev-window sanitiser, in the order they are applied
// ---------------------------------------------------------------------------

/// Rule 1: reset non-finite or inverted windows.
void repair_degenerate(double emin, double emax, double ecut_floor, ModuleIO::CheRange& out)
{
    const bool finite_in = std::isfinite(emin) && std::isfinite(emax);
    if (finite_in && (emax - emin) > 0.0)
    {
        return;
    }
    const double abs_emax = std::isfinite(emax) ? std::fabs(emax) : 0.0;
    out.emax = std::max(ecut_floor, abs_emax);
    out.emin = fallback_emin(out.emax);
    out.degenerate = true;
    out.emin_repaired = true;
    out.emax_repaired = true;
}

/// Rule 2: the window has to contain every known deterministic eigenvalue.
void cover_deterministic_eigenvalues(const std::vector<double>& det_eigs, ModuleIO::CheRange& out)
{
    if (det_eigs.empty())
    {
        return;
    }
    const double det_min = *std::min_element(det_eigs.begin(), det_eigs.end());
    const double det_max = *std::max_element(det_eigs.begin(), det_eigs.end());
    if (std::isfinite(det_min) && out.emin > det_min - det_margin)
    {
        out.emin = det_min - det_margin;
        out.emin_repaired = true;
    }
    if (std::isfinite(det_max) && out.emax < det_max + det_margin)
    {
        out.emax = det_max + det_margin;
        out.emax_repaired = true;
    }
}

/// Rule 3: a positive lower bound cannot describe a spectrum with occupied
/// states; replace it by a conservative negative value.
void repair_positive_emin(ModuleIO::CheRange& out)
{
    if (out.emin <= 0.0)
    {
        return;
    }
    out.emin = fallback_emin(out.emax);
    out.emin_repaired = true;
}

/// Rule 4: the highest plane-wave eigenvalue is bounded by roughly ecutwfc.
void raise_emax_to_ecut_floor(double ecut_floor, ModuleIO::CheRange& out)
{
    if (ecut_floor > 0.0 && out.emax < ecut_floor)
    {
        out.emax = ecut_floor;
        out.emax_repaired = true;
    }
}

/// Rule 5: guarantee a finite, non-degenerate width.
void enforce_min_width(ModuleIO::CheRange& out)
{
    if (!(out.emax - out.emin >= min_width))
    {
        out.emax = out.emin + min_width;
        out.emax_repaired = true;
    }
}

} // namespace

namespace ModuleIO
{

CheRange sanitize_che_range(double emin, double emax, double ecutwfc, const std::vector<double>& det_eigs)
{
    CheRange out;
    out.emin = emin;
    out.emax = emax;

    const double ecut_floor = (ecutwfc > 0.0) ? ecut_margin * ecutwfc : 0.0;

    repair_degenerate(emin, emax, ecut_floor, out);
    cover_deterministic_eigenvalues(det_eigs, out);
    repair_positive_emin(out);
    raise_emax_to_ecut_floor(ecut_floor, out);
    enforce_min_width(out);

    return out;
}

std::vector<std::string> check_dos_parameters(double range_ev,
                                              double sigma_ev,
                                              double de_ev,
                                              int dos_nche)
{
    std::vector<std::string> warnings;

    if (de_ev <= 0.0)
    {
        warnings.push_back("dos_edelta_ev must be positive; the DOS grid cannot be built.");
        return warnings;
    }
    if (!(sigma_ev > 0.0))
    {
        warnings.push_back("dos_sigma = " + format_double(sigma_ev)
                           + " eV is not positive; a zero-width Gaussian makes the DOS a sum of delta "
                             "spikes that no finite dos_nche can resolve.");
        return warnings;
    }

    // Grid vs broadening: the two must be comparable, otherwise the curve is
    // either oversampled noise or undersampled peaks.
    if (de_ev < de_over_sigma_min * sigma_ev)
    {
        warnings.push_back("dos_edelta_ev = " + format_double(de_ev) + " eV is much smaller than dos_sigma = "
                           + format_double(sigma_ev) + " eV; the extra grid points carry no independent "
                             "information and will show Chebyshev ringing. Consider dos_edelta_ev ~ dos_sigma.");
    }
    if (de_ev > de_over_sigma_max * sigma_ev)
    {
        warnings.push_back("dos_edelta_ev = " + format_double(de_ev) + " eV is much larger than dos_sigma = "
                           + format_double(sigma_ev) + " eV; the Gaussian is narrower than the grid spacing, so "
                             "peaks can fall between grid points. Use dos_edelta_ev <~ dos_sigma.");
    }

    // Chebyshev order needed to resolve a Gaussian of width sigma on this range.
    if (range_ev > 0.0)
    {
        const double nche_needed = 0.5 * range_ev / sigma_ev * order_prefactor;
        if (nche_needed > static_cast<double>(dos_nche))
        {
            std::ostringstream os;
            os << "dos_sigma = " << sigma_ev << " eV is narrower than the Chebyshev expansion can resolve: "
               << "range = " << range_ev << " eV, dos_nche = " << dos_nche
               << ", estimated order needed ~ " << static_cast<long long>(std::ceil(nche_needed))
               << ". Increase dos_nche or use a larger dos_sigma; otherwise the DOS may show spurious "
                  "oscillations or negative values.";
            warnings.push_back(os.str());
        }
        if (range_ev / sigma_ev > 1.0e4)
        {
            warnings.push_back("the Chebyshev range (" + format_double(range_ev)
                               + " eV) is more than 1e4 times dos_sigma; with such a wide range the required "
                                 "Chebyshev order is prohibitive and the DOS will be under-resolved. Consider "
                                 "restricting dos_emin_ev/dos_emax_ev to the energy window of interest.");
        }
    }

    return warnings;
}

std::vector<std::string> check_dos_result(double min_dos, double max_abs_t)
{
    std::vector<std::string> warnings;

    if (max_abs_t > 1.0 + max_abs_t_tol)
    {
        warnings.push_back("an output energy lies outside the Chebyshev interval (max |t| = "
                           + format_double(max_abs_t)
                           + " > 1). T_n(t) grows exponentially there, so the DOS is not trustworthy; "
                             "check dos_emin_ev / dos_emax_ev against the interval printed above.");
    }

    if (std::isfinite(min_dos) && min_dos < 0.0)
    {
        // Allow the tiny ringing of a well-resolved expansion, flag anything visible.
        if (std::fabs(min_dos) > negative_dos_abs_tol)
        {
            warnings.push_back("the DOS takes negative values (minimum " + format_double(min_dos)
                               + "): the Chebyshev expansion is under-resolved. Increase dos_nche or use a "
                                 "larger dos_sigma.");
        }
    }

    return warnings;
}

} // namespace ModuleIO
