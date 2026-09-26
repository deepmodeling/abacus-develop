#ifndef SOURCE_DFTB_SKF_DATA_H
#define SOURCE_DFTB_SKF_DATA_H

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace ModuleDFTB
{

/** Values and first radial derivative of a pair repulsive potential. */
struct SkfRepulsiveValue
{
    double energy_hartree = 0.0;
    double derivative_hartree_per_bohr = 0.0;
};

/** Values and the first two radial derivatives of all ten legacy SK channels. */
struct SkfIntegralEvaluation
{
    std::array<double, 10> hamiltonian{};
    std::array<double, 10> overlap{};
    std::array<double, 10> d_hamiltonian_dr{};
    std::array<double, 10> d_overlap_dr{};
    std::array<double, 10> d2_hamiltonian_dr2{};
    std::array<double, 10> d2_overlap_dr2{};
};

/**
 * Repulsive spline in the legacy DFTB Slater-Koster format.
 * Distances are in Bohr and energies/coefficients in Hartree.
 */
struct SkfRepulsiveSpline
{
    double cutoff_bohr = 0.0;
    std::array<double, 3> exponential_coefficients{};
    std::vector<double> interval_starts_bohr;
    std::vector<double> interval_ends_bohr;
    std::vector<std::array<double, 4>> cubic_coefficients;
    std::array<double, 6> tail_coefficients{};

    SkfRepulsiveValue evaluate(double distance_bohr) const;
};

/**
 * One legacy (s/p/d) SKF file. The ten on-disk integral columns follow the
 * v1.0 legacy order: dd-sigma, dd-pi, dd-delta, pd-sigma, pd-pi, pp-sigma,
 * pp-pi, sd-sigma, sp-sigma, ss-sigma. SKF distance/energy units are
 * Bohr/Hartree.
 *
 * This reader intentionally rejects the newer '@' extended (s/p/d/f) format
 * instead of silently interpreting its different column layout. Integral
 * interpolation follows DFTB+'s default (new) 8-point polynomial scheme and
 * quintic cutoff smoothing with a one-Bohr distance extension.
 */
struct SkfData
{
    std::string filename;
    bool homonuclear = false;
    double grid_spacing_bohr = 0.0;
    std::size_t declared_grid_points = 0;
    std::vector<std::array<double, 10>> hamiltonian;
    std::vector<std::array<double, 10>> overlap;

    // Atomic data appear only in homonuclear files. Shell order is s, p, d.
    std::array<double, 3> onsite_hartree{};
    std::array<double, 3> hubbard_u_hartree{};
    std::array<double, 3> reference_occupations{};
    double mass_amu = 0.0;
    bool has_atomic_data = false;

    bool has_repulsive_spline = false;
    SkfRepulsiveSpline repulsive;

    static SkfData read_legacy(const std::string& filename, bool homonuclear);

    /** Neutral-atom valence count from the s, p, d reference shell populations. */
    double valence_electron_count() const;

    /** Evaluate all ten H/S channels, including DFTB+ default cutoff handling. */
    SkfIntegralEvaluation evaluate_integrals(double distance_bohr) const;
};

/**
 * Validate metadata shared by A-B and B-A files. This mirrors the constraints
 * used by DFTB+ for the legacy format: matching electronic grid and the same
 * pair repulsive function. Directional electronic integrals are not compared.
 */
void validate_pair_directions(const SkfData& ab, const SkfData& ba, double tolerance = 1.0e-8);

} // namespace ModuleDFTB

#endif // SOURCE_DFTB_SKF_DATA_H
