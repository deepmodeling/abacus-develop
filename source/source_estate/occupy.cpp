#include "occupy.h"

#include <limits>
#include <sstream>

#include "source_base/constants.h"
#include "source_base/mymath.h"
#include "source_base/parallel_reduce.h"
#include "source_io/module_parameter/parameter.h"

namespace
{
// A chemical potential with an explicit energy origin. Occupations must use
// both parts; the absolute double is only the rounded public/output value.
struct FermiLevel
{
    double reference;
    double offset;
};

FermiLevel find_fermi_level(const ModuleBase::matrix& ekb,
                           const int nband,
                           const int nks,
                           const double& nelec,
                           const std::vector<double>& wk,
                           const double& smearing_sigma,
                           const int ngauss,
                           const int& is,
                           const std::vector<int>& isk);

double sum_occupations(const ModuleBase::matrix& ekb,
                       const int nband,
                       const int nks,
                       const std::vector<double>& wk,
                       const double& smearing_sigma,
                       const int ngauss,
                       const double reference,
                       const double offset,
                       const int& is,
                       const std::vector<int>& isk,
                       double* roundoff);

// First-order scalar error estimate, used only when the Fermi search stalls.
// The 8u/4u allowances for erfc/exp are tested libm assumptions, not a portable
// correctly-rounded guarantee. Keep MP component magnitudes before cancellation.
double smearing_roundoff(const double x, const int n, const double value, const double x_error)
{
    const double u = 0.5 * std::numeric_limits<double>::epsilon();
    const double tiny = std::numeric_limits<double>::denorm_min();
    if (n == -99)
    {
        const double magnitude = std::abs(x);
        if (magnitude > 200.0)
        {
            const double tail_argument = -magnitude;
            return std::exp(tail_argument) + tiny;
        }
        const double slope = value * (1.0 - value);
        return u * (2.0 * std::abs(value) + 4.0 * slope) + slope * x_error + tiny;
    }
    const double negative_x = -x;
    const double gaussian = 0.5 * std::erfc(negative_x);
    const double x2 = x * x;
    const double negative_x2 = -x2;
    const double density = std::exp(negative_x2) / std::sqrt(ModuleBase::PI);
    if (n == 0)
    {
        return 8.0 * u * gaussian + density * x_error + tiny;
    }
    if (n == -1)
    {
        const double shift = 1.0 / ModuleBase::SQRT2;
        const double xp = x - shift;
        const double xp2 = xp * xp;
        const double negative_xp = -xp;
        const double base = 0.5 * std::erfc(negative_xp);
        const double argument = std::min(200.0, xp2);
        const double negative_argument = -argument;
        const double correction = std::exp(negative_argument) / std::sqrt(ModuleBase::TWO_PI);
        if (xp2 > 200.0)
        {
            // The existing exponential clamp leaves a nonzero tail.
            return 8.0 * u * base + correction + tiny;
        }
        const double derivative = (ModuleBase::SQRT2 - 2.0 * xp) * correction;
        const double slope = std::abs(derivative);
        const double shifted_error = x_error + u * (std::abs(xp) + 2.0 * shift);
        // Square, exp, normalization and addition; include the shift's rounding.
        return u * (8.0 * base + (xp2 + 10.0) * correction + std::abs(value)) + slope * shifted_error + tiny;
    }
    if (n != 1 && n != 2)
    {
        return std::numeric_limits<double>::infinity();
    }
    const double magnitude = std::abs(x);
    double polynomial = 0.5 * magnitude;
    double derivative = 1.5 - x2;
    double slope = std::abs(derivative) * density;
    double operations = 8.0;
    if (n == 2)
    {
        polynomial = magnitude * (0.875 + 0.25 * x2);
        derivative = 1.875 - 2.5 * x2 + 0.5 * x2 * x2;
        slope = std::abs(derivative) * density;
        operations = 12.0;
    }
    const double argument = std::min(200.0, x2);
    const double negative_argument = -argument;
    const double correction_scale = polynomial * std::exp(negative_argument) / std::sqrt(ModuleBase::PI);
    if (x2 > 200.0)
    {
        return 8.0 * u * gaussian + correction_scale + tiny;
    }
    // 8/12 count first-order rounding allowances in the MP1/MP2 recurrence,
    // including the exponential and normalization, not an electron tolerance.
    return u * (8.0 * gaussian + (operations + x2) * correction_scale + n * std::abs(value))
           + slope * x_error + tiny;
}
} // namespace

Occupy::Occupy()
{
}
Occupy::~Occupy()
{
}

//===========================================================
// Four smearing methods:
// (1) do nothing
// (2) gaussian_broadening method (need 'smearing_sigma' value).
// (3) tetrahedron method
// (4) fixed_occupations
//===========================================================

bool Occupy::use_gaussian_broadening = false;
int Occupy::gaussian_type = 0;
double Occupy::gaussian_parameter = 0.01;
bool Occupy::fixed_occupations = false;

void Occupy::decision(const std::string& name, const std::string& smearing_method, const double& smearing_sigma)
{
    ModuleBase::TITLE("Occupy", "decision");
    use_gaussian_broadening = false;
    fixed_occupations = false;

    gaussian_type = 0;
    gaussian_parameter = smearing_sigma;

    if (name == "fixed")
    {
        if (gaussian_parameter != 0.0)
        {
            ModuleBase::WARNING("smearing_decision", "Fixed occupations,gauss broadening ignored");
            ModuleBase::GlobalFunc::AUTO_SET("gaussian_parameter", 0.0);
            gaussian_parameter = 0.0;
        }
    }
    else if (name == "smearing" && smearing_method == "fixed")
    {
        if (gaussian_parameter != 0.0)
        {
            ModuleBase::WARNING("smearing_decision", "Fixed occupations,gauss broadening ignored");
            ModuleBase::GlobalFunc::AUTO_SET("gaussian_parameter", 0.0);
            gaussian_parameter = 0.0;
        }
    }

    // there are four types of smearing methods:
    // (1) gaussian
    // (2) methfessel-paxton
    // (3) Marzari-Vanderbilt
    // (4) Fermi-Dirac
    else if (name == "smearing")
    {
        use_gaussian_broadening = true;
        if (gaussian_parameter == 0.0)
        {
            ModuleBase::WARNING_QUIT(
                "smearing_decision",
                "Smearing requires gaussian broadening,but gaussian_parameter = 0(default value = 0.01)");
        }
        if (smearing_method == "gaussian" || smearing_method == "gauss")
        {
            gaussian_type = 0; //  0: gaussian
        }
        else if (smearing_method == "methfessel-paxton" || smearing_method == "mp")
        {
            gaussian_type = 1; // >0 Methfessel-Paxton method.
        }
        else if (smearing_method == "mp2")
        {
            gaussian_type = 2; // 2nd Methfessel-Paxton method.
        }
        else if (smearing_method == "mp3")
        {
            // acually any order Methfessel-Paxton method can be supported in occupy_smearing::w1gauss()
            // however the parameter is string instead of int
            ModuleBase::WARNING_QUIT(
                "occupy",
                "Some refactor of smearing shoule be done before supporting any order of Methfessel-Paxton method!");
        }

        else if (smearing_method == "marzari-vanderbilt" || smearing_method == "cold" || smearing_method == "mv")
        {
            gaussian_type = -1;
        }
        else if (smearing_method == "fermi-dirac" || smearing_method == "fd")
        {
            gaussian_type = -99;
        }
    }
    else if (name == "tetrahedra")
    {
        ModuleBase::WARNING_QUIT("occupy", "not implemented yet!");
    }
    else if (name == "from_input")
    {
        fixed_occupations = true;
    }
    else
    {
        ModuleBase::WARNING_QUIT("occupy_decision", "occupations, not implemented");
    }
    return;
}

/**
 * @brief calculates weights and fermi energy for semiconductors and insulators
 *
 * @note It is only applicable to semiconductors and insulators that have an energy gap.
 *       Here, bands are either empty or filled, and there is no electronic transfer among different k-points.
 *
 * @param nks number of k points.
 * @param wk weight of each k point (consider symmetry).
 * @param nbands number of locally owned bands.
 * @param band_offset global index of the first locally owned band.
 * @param nelec number of electrons for this spin direction.
 * @param ekb the array save the band energy.
 * @param ef output: the highest occupied Kohn-Sham level.
 * @param wg output: weight for each k, each band.
 * @param is the spin index now.
 * @param isk distinguish k point belong to which spin.
 * @param nspin number of spin channels; 2 selects a single spin through isk when
 *              is != -1, and 4 halves the spin degeneracy.
 */
void Occupy::iweights(
    const int nks,
    const std::vector<double>& wk,
    const int nbands,
    const int band_offset,
    const double& nelec,
    const ModuleBase::matrix& ekb,
    double& ef,
    ModuleBase::matrix& wg,
    const int& is, //<- is should be -1, 0, or 1. -1 means set all spins, and 0 means spin up, 1 means spin down.
    const std::vector<int>& isk,
    const int nspin)
{
    assert(is < 2);
    double degspin = 2.0;
    if (nspin == 4) {
        degspin = 1.0;
}
    if (is != -1) {
        degspin = 1.0;
}

    double ib_mind = nelec / degspin;
    int ib_min = std::ceil(ib_mind);
    if (ib_min != int(ib_mind))
    {
        ModuleBase::WARNING_QUIT("iweights",
                                 "It is not a semiconductor or insulator. Please do not set 'smearing_method=fixed', "
                                 "and try other options.");
    }
    ef = -1e+10;

    for (int ik = 0; ik < nks; ++ik)
    {
        // when NSPIN=2, only calculate spin up or spin down with TWO_FERMI mode(nupdown != 0)
        if (nspin == 2 && isk[ik] != is && is != -1)
        {
            continue;
        }

        for (int ib = 0; ib < nbands; ++ib)
        {
            const int global_band = band_offset + ib;
            if (global_band < ib_min)
            {
                wg(ik, ib) = wk[ik];
                ef = std::max(ef, ekb(ik, ib));
            }
            else
            {
                wg(ik, ib) = 0.0;
            }
        }
    }
    #ifdef __MPI
    Parallel_Reduce::reduce_max(ef);
    #endif
    return;
}

/**
 * @brief calculates weights with the gaussian spreading technique
 *
 * @param nks number of k points.
 * @param wk weight of each k point(symmetry considered).
 * @param nband number of bands.
 * @param nelec number of electrons.
 * @param smearing_sigma parameter input by user.
 * @param ngauss which type of smearing.
 * @param ekb band energy.
 * @param ef output: rounded absolute Fermi level
 * @param demet output: energy correction for metal
 * @param wg output: weight of each band at each k point
 * @param is spin
 * @param isk array to point out each k belong to which spin
 */
void Occupy::gweights(const int nks,
                      const std::vector<double>& wk,
                      const int nband,
                      const double& nelec,
                      const double& smearing_sigma,
                      const int ngauss,
                      const ModuleBase::matrix& ekb,
                      double& ef,
                      double& demet,
                      ModuleBase::matrix& wg,
                      const int& is,
                      const std::vector<int>& isk)
{
    // ModuleBase::TITLE("Occupy","gweights");
    //===============================
    //  Calculate the Fermi energy ef
    //===============================
    // Keep the centered solution when constructing the occupations.
    const FermiLevel level = find_fermi_level(ekb, nband, nks, nelec, wk, smearing_sigma, ngauss, is, isk);
    ef = level.reference + level.offset;
    demet = 0.0;

    for (int ik = 0; ik < nks; ik++)
    {
        // mohan add 2011-04-03
        if (is != -1 && is != isk[ik]) {
            continue;
}

        for (int ib = 0; ib < nband; ib++)
        {
            //================================
            // Calculate the gaussian weights
            //================================
            // call wgauss
            const double relative_energy = ekb(ik, ib) - level.reference;
            const double argument = (level.offset - relative_energy) / smearing_sigma;
            wg(ik, ib) = wk[ik] * occupy_smearing::wgauss(argument, ngauss);

            //====================================================================
            // The correct form of the band energy is  \int e n(e) de   for e<ef
            // demet is the correction to add to the sum of eigenvalues
            //====================================================================
            // Mohan fix bug 2010-1-9
            // call w1gauss
            demet += wk[ik] * smearing_sigma * occupy_smearing::w1gauss(argument, ngauss);
        }
    }

    return;
} // end subroutine gweights

/**
 * @brief Finds the Fermi energy in a centered coordinate
 *
 * @param ekb band energy.
 * @param nband number of bands.
 * @param nks nks number of k points.
 * @param nelec number of electrons.
 * @param wk wk weight of each k point(symmetry considered).
 * @param smearing_sigma parameter input by user.
 * @param ngauss which type of smearing.
 * @return energy reference and converged relative Fermi level
 * @param is spin
 * @param isk array to point out each k belong to which spin
 */
namespace
{
FermiLevel find_fermi_level(const ModuleBase::matrix& ekb,
                           const int nband,
                           const int nks,
                           const double& nelec,
                           const std::vector<double>& wk,
                           const double& smearing_sigma,
                           const int ngauss,
                           const int& is,
                           const std::vector<int>& isk)
{
    // Resolve particle number to floating-point accuracy, independently of the
    // outer SCF threshold. A fixed absolute tolerance leaves a charge floor.
    const double electron_magnitude = std::abs(nelec);
    const double scale = std::max(1.0, electron_magnitude);
    const double tolerance = 8.0 * std::numeric_limits<double>::epsilon() * scale;
    if (!(smearing_sigma > 0.0) || !std::isfinite(smearing_sigma) || !std::isfinite(nelec))
    {
        ModuleBase::WARNING_QUIT("occupy_smearing::efermig", "Invalid smearing width or electron number");
    }

    // Empty local partitions contribute neutral extrema to the global bounds.
    double lower = std::numeric_limits<double>::infinity();
    double upper = -std::numeric_limits<double>::infinity();
    for (int ik = 0; ik < nks; ++ik)
    {
        if (is != -1 && is != isk[ik])
        {
            continue;
        }
        for (int ib = 0; ib < nband; ++ib)
        {
            const double eigenvalue = ekb(ik, ib);
            lower = std::min(lower, eigenvalue);
            upper = std::max(upper, eigenvalue);
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_min(lower);
    Parallel_Reduce::reduce_max(upper);
#endif
    if (!std::isfinite(lower) || !std::isfinite(upper))
    {
        ModuleBase::WARNING_QUIT("occupy_smearing::efermig", "No finite eigenvalue bounds");
    }

    // Keep the energy origin separate from the small chemical-potential offset.
    // Rounding reference + offset before evaluating occupations can destroy
    // particle-number accuracy even for a harmless constant shift of all levels.
    double reference = 0.5 * lower + 0.5 * upper;
    lower -= reference;
    upper -= reference;
    double step = 2.0 * smearing_sigma;
    lower -= step;
    upper += step;
    double low_count = sum_occupations(ekb, nband, nks, wk, smearing_sigma, ngauss, reference, lower, is, isk, nullptr);
    double high_count = sum_occupations(ekb, nband, nks, wk, smearing_sigma, ngauss, reference, upper, is, isk, nullptr);
    int expansion = 0;
    while (low_count > nelec + tolerance || high_count < nelec - tolerance)
    {
        if (++expansion > 100 || !std::isfinite(step))
        {
            ModuleBase::WARNING_QUIT("occupy_smearing::efermig", "Cannot bracket the requested electron number");
        }
        step *= 2.0;
        if (low_count > nelec + tolerance)
        {
            lower -= step;
            low_count = sum_occupations(ekb, nband, nks, wk, smearing_sigma, ngauss, reference, lower, is, isk, nullptr);
        }
        if (high_count < nelec - tolerance)
        {
            upper += step;
            high_count = sum_occupations(ekb, nband, nks, wk, smearing_sigma, ngauss, reference, upper, is, isk, nullptr);
        }
    }

    if (!std::isfinite(low_count) || !std::isfinite(high_count)
        || !std::isfinite(lower) || !std::isfinite(upper))
    {
        ModuleBase::WARNING_QUIT("occupy_smearing::efermig", "Non-finite Fermi bracket");
    }

    // Keep the best evaluated point, including the saturated endpoints. For
    // empty/full bands or a resolved gap, the chemical potential is not unique.
    FermiLevel result{reference, lower};
    double best_error = std::abs(low_count - nelec);
    const double high_error = std::abs(high_count - nelec);
    if (high_error < best_error)
    {
        result.offset = upper;
        best_error = high_error;
    }
    // Inherited safety cap, not an accuracy criterion. Stagnation is checked
    // separately; exhausting this cap must never report convergence.
    const int maxiter = 300;
    for (int iteration = 0; iteration < maxiter; ++iteration)
    {
        if (best_error <= tolerance)
        {
            return result;
        }
        // Follow the shrinking bracket: a fixed spectral midpoint can be far
        // from the Fermi surface when there are high unoccupied bands. Keep
        // the small interval in offsets even when the absolute origin rounds.
        const double trial_offset = 0.5 * lower + 0.5 * upper;
        const double next_reference = reference + trial_offset;
        const double shift = next_reference - reference;
        if (shift != 0.0 && std::isfinite(next_reference))
        {
            lower -= shift;
            upper -= shift;
            reference = next_reference;
        }
        const double middle = 0.5 * lower + 0.5 * upper;
        if (middle == lower || middle == upper)
        {
            // Audit only the best candidate, in its original centered coordinates.
            // A large error estimate must not buy a loose electron tolerance.
            double roundoff = 0.0;
            const double count = sum_occupations(ekb, nband, nks, wk, smearing_sigma, ngauss,
                                                 result.reference, result.offset, is, isk, &roundoff);
            const double count_error = count - nelec;
            const double residual = std::abs(count_error);
            std::ostringstream message;
            message.precision(17);
            if (std::isfinite(roundoff) && roundoff <= tolerance && residual <= tolerance + roundoff)
            {
                message << "Fermi search reached rounding accuracy; electron residual = " << residual
                        << ", target = " << tolerance << ", estimated roundoff = " << roundoff;
                const std::string diagnostic = message.str();
                ModuleBase::WARNING("occupy_smearing::efermig", diagnostic);
                return result;
            }
            message << "Fermi bracket reached floating-point resolution; electron residual = "
                    << residual << ", requested tolerance = " << tolerance << ", estimated roundoff = " << roundoff;
            const std::string diagnostic = message.str();
            ModuleBase::WARNING_QUIT("occupy_smearing::efermig", diagnostic);
            return result;
        }
        const double count = sum_occupations(ekb, nband, nks, wk, smearing_sigma, ngauss, reference, middle, is, isk, nullptr);
        if (!std::isfinite(count))
        {
            ModuleBase::WARNING_QUIT("occupy_smearing::efermig", "Non-finite electron number");
        }
        const double error = std::abs(count - nelec);
        if (error < best_error)
        {
            result.reference = reference;
            result.offset = middle;
            best_error = error;
        }
        if (best_error <= tolerance)
        {
            return result;
        }
        if (count < nelec)
        {
            lower = middle;
        }
        else
        {
            upper = middle;
        }
    }
    ModuleBase::WARNING_QUIT("occupy_smearing::efermig", "Fermi search exceeded its iteration limit");
    return result;
}

} // namespace

void occupy_smearing::efermig(const ModuleBase::matrix& ekb,
                             const int nband,
                             const int nks,
                             const double& nelec,
                             const std::vector<double>& wk,
                             const double& smearing_sigma,
                             const int ngauss,
                             double& ef,
                             const int& is,
                             const std::vector<int>& isk)
{
    const FermiLevel level = find_fermi_level(ekb, nband, nks, nelec, wk, smearing_sigma, ngauss, is, isk);
    // Public output is rounded to an absolute double; gweights retains both parts.
    ef = level.reference + level.offset;
}

/**
 * @brief This function computes the number of states under a given energy e
 *
 * @param ekb band energy.
 * @param nband number of bands.
 * @param nks number of k points.
 * @param wk weight of each k point(symmetry considered).
 * @param smearing_sigma parameter input by user
 * @param ngauss which type of smearing.
 * @param reference common energy origin
 * @param offset Fermi level relative to the common origin
 * @param is spin
 * @param isk array to point out each k belong to which spin
 * @return (double) the number of states
 */
namespace
{
double sum_occupations(const ModuleBase::matrix& ekb,
                       const int nband,
                       const int nks,
                       const std::vector<double>& wk,
                       const double& smearing_sigma,
                       const int ngauss,
                       const double reference,
                       const double offset,
                       const int& is,
                       const std::vector<int>& isk,
                       double* roundoff)
{
    // ModuleBase::TITLE("Occupy","sumkg");
    double sum2 = 0.0;
    double compensation = 0.0;
    const double u = 0.5 * std::numeric_limits<double>::epsilon();
    double absolute_sum = 0.0;
    double evaluation_error = 0.0;
    for (int ik = 0; ik < nks; ++ik)
    {
        if (is != -1 && is != isk[ik])
        {
            continue;
        }
        for (int ib = 0; ib < nband; ++ib)
        {
            // Sum the same weighted occupations as gweights, compensating the
            // local accumulation before the existing pool/band reduction.
            const double relative_energy = ekb(ik, ib) - reference;
            const double argument = (offset - relative_energy) / smearing_sigma;
            const double value = occupy_smearing::wgauss(argument, ngauss);
            const double occupation = wk[ik] * value;
            const double corrected = occupation - compensation;
            const double next = sum2 + corrected;
            compensation = (next - sum2) - corrected;
            sum2 = next;
            if (roundoff != nullptr)
            {
                const double magnitude = std::abs(occupation);
                absolute_sum += magnitude;
                const double difference = offset - relative_energy;
                const double coordinate_scale = std::abs(relative_energy) / smearing_sigma
                                                + std::abs(difference) / smearing_sigma + std::abs(argument);
                const double x_error = u * coordinate_scale;
                const double scalar_error = smearing_roundoff(argument, ngauss, value, x_error);
                evaluation_error += std::abs(wk[ik]) * scalar_error + u * magnitude;
            }
        }
    }

    // GlobalV::ofs_running << "\n sum2 before reduce = " << sum2 << std::endl;

    double totals[3] = {sum2, absolute_sum, evaluation_error};
    double reduction_error = 0.0;
#ifdef __MPI
    const int npool = GlobalV::KPAR * PARAM.inp.bndpar;
    const int count = roundoff == nullptr ? 1 : 3;
    Parallel_Reduce::reduce_double_allpool(npool, GlobalV::NPROC_IN_POOL, totals, count);
    if (roundoff != nullptr && npool > 1)
    {
        int participants = 1;
        Parallel_Reduce::reduce_all(participants);
        const double additions = (participants - 1) * u;
        // Pool scaling plus an order-independent reduction bound; used only
        // for the capped stagnation check, never for early termination.
        reduction_error = (u + additions / (1.0 - additions)) * totals[1];
    }
#endif
    if (roundoff != nullptr)
    {
        // First-order Kahan bound, function/argument evaluation, MPI and the
        // final count subtraction. Subnormal widths invalidate this model.
        *roundoff = totals[2] + 2.0 * u * totals[1] + reduction_error + u * std::abs(totals[0]);
        if (smearing_sigma < std::numeric_limits<double>::min())
        {
            *roundoff = std::numeric_limits<double>::infinity();
        }
    }

    // GlobalV::ofs_running << "\n sum2 after reduce = " << sum2 << std::endl;

    return totals[0];
}
} // namespace

double occupy_smearing::sumkg(const ModuleBase::matrix& ekb,
                             const int nband,
                             const int nks,
                             const std::vector<double>& wk,
                             const double& smearing_sigma,
                             const int ngauss,
                             const double& e,
                             const int& is,
                             const std::vector<int>& isk)
{
    const double reference = 0.0;
    return sum_occupations(ekb, nband, nks, wk, smearing_sigma, ngauss, reference, e, is, isk, nullptr);
}

double occupy_smearing::wgauss(const double& x, const int n)
{
    // ModuleBase::TITLE("Occupy","wgauss");
    //=====================================================================
    //  This function computes the approximate theta function for the
    //  iven order n, at the point x.
    //  --> (n>=0) : Methfessel-Paxton case. See PRB 40, 3616 (1989).
    //  --> (n=-1 ): Cold smearing
    //               (Marzari-Vanderbilt). See PRL 82, 3296 (1999)
    //        1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
    //  --> (n=-99): Fermi-Dirac case: 1.0/(1.0+exp(-x)).
    //=====================================================================
    //============================
    //  The value of this function
    //============================
    double wga = 0.0;
    const double maxarg = 200.0;

    //===========================
    // Fermi-Dirac(fd) smearing
    //===========================
    if (n == -99)
    {
        if (x < -maxarg)
        {
            wga = 0.0;
        }
        else if (x > maxarg)
        {
            wga = 1.0;
        }
        else
        {
            wga = 1.00 / (1.0 + std::exp(-x));
        }
        return wga;
    }

    //===================
    // Cold smearing(mv)
    //===================
    if (n == -1)
    {
        const double xp = x - 1.00 / ModuleBase::SQRT2;
        const double arg = std::min(maxarg, xp * xp);
        const double negative_xp = -xp;
        wga = 0.50 * std::erfc(negative_xp) + 1.00 / sqrt(ModuleBase::TWO_PI) * std::exp(-arg);
        return wga;
    }

    //====================
    // Methfessel-Paxton       //pengfei 2014-10-13
    //====================
    const double negative_x = -x;
    wga = 0.5 * std::erfc(negative_x);
    // wga = gauss_freq(x * ModuleBase::SQRT2);
    //	std::cout<<"\n x="<<x<<" wga="<<wga;
    if (n == 0)
    {
        return wga;
    }

    // double hd = 0.0;

    int ni = 0;
    const double arg = std::min(maxarg, x * x);
    double hp = std::exp(-arg);
    double h0 = 1.00;
    double h1 = -2.00 * x;
    double a = 1.0 / sqrt(ModuleBase::PI);
    for (int i = 0; i < n; i++)
    {
        a = -a / (static_cast<double>(i + 1) * 4.00);
        wga = wga + a * h1 * hp;
        ++ni;
        // std::cout << " wga = " <<wga<<std::endl;
        h0 = 2.00 * (-x) * h1 - 2.00 * static_cast<double>(ni) * h0;
        ++ni;

        h1 = 2.00 * (-x) * h0 - 2.00 * static_cast<double>(ni) * h1;
    }
    return wga;
} // end function wgauss

double occupy_smearing::w1gauss(const double& x, const int n)
{
    //========================================================================
    //    w1gauss(x,n) = \int_{-\infty}^x   y delta(y) dy
    //    where delta(x) is the current approximation for the delta function,
    //    as obtained from w0gauss(x,n)
    //
    // --> (n>=0) : Methfessel-Paxton case
    //
    // --> (n=-1): Cold smearing (Marzari-Vanderbilt)
    //     w1gauss = 1/sqrt(2*pi)*(x-1/sqrt(2))*exp(-(x-1/sqrt(2))**2)
    //
    // --> (n=-99): Fermi-Dirac case. In this case w1gauss corresponds
    //     to the negative of the electronic entropy.
    //========================================================================
    double w1 = 0.0;

    //=======================
    // Fermi-Dirac smearing
    //=======================
    if (n == -99)
    {
        if (std::abs(x) <= 36.0)
        {
            const double f = 1.00 / (1.00 + exp(-x));
            const double onemf = 1.00 - f;
            w1 = f * log(f) + onemf * log(onemf);
            //==================================================
            // in order to avoid problems for large values of x
            //==================================================
        }
        else
        {
            //=============================================
            // neglect w1gauss when abs(w1gauss) < 1.0d-14
            //=============================================
            w1 = 0.00;
        }
        return w1;
    }
    //===============
    // Cold smearing
    //===============

    if (n == -1)
    {
        const double xp = x - 1.00 / ModuleBase::SQRT2;
        const double arg = std::min(200.0, xp * xp);
        w1 = 1.00 / sqrt(ModuleBase::TWO_PI) * xp * std::exp(-arg);
        return w1;
    }

    //====================
    // Methfessel-Paxton
    //====================
    const double arg = std::min(200.0, x * x);
    w1 = -0.50 * std::exp(-arg) / sqrt(ModuleBase::PI);

    // std::cout << "\n x=" << x << " w1=" << w1;
    if (n == 0) // specific case : gaussian smearing.
    {
        return w1;
    }

    /*    double hd = 0.0;
        double hp = exp(- arg);
        int ni = 0;
        double a = 1.0 / sqrt(ModuleBase::PI);

        for (int i = 0;i < n;i++)
        {
            hd = 2.00 * x * hp - 2.00 * static_cast<double>(ni) * hd;//dble(ni)
            ni ++;
            const double hpm1 = hp;
            hp = 2.00 * x * hd - 2.00 * static_cast<double>(ni) * hp;//dble(ni)
            ni ++;
            // mohan fixed bug 2010-1-10
            // description of bug: i should not be 0.
            a = - a / ( static_cast<double>(i+1) * 4.00);//dble(i)
            std::cout << " w1 == "<<w1<<std::endl;
            w1 -= a * (0.50 * hp + static_cast<double>(ni) * hpm1);//dble(ni)
            std::cout << " w1 == "<<w1<<std::endl;
        }*/

    int ni = 0;
    double hp = std::exp(-arg);
    double h0 = 1.00;
    double h1 = 2.00 * x;
    double a = 1.0 / sqrt(ModuleBase::PI);
    for (int i = 0; i < n; i++)
    {
        a = -a / (static_cast<double>(i + 1) * 4.00);
        ni++;
        const double h0m1 = h0;
        h0 = 2.00 * x * h1 - 2.00 * static_cast<double>(ni) * h0;
        ni++;
        h1 = 2.00 * x * h0 - 2.00 * static_cast<double>(ni) * h1;
        w1 = w1 - a * (0.50 * h0 + static_cast<double>(ni) * h0m1) * hp;
    }

    return w1;
} // end function w1gauss
