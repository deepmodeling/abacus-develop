#ifndef SOURCE_DFTB_EIGENSOLVER_H
#define SOURCE_DFTB_EIGENSOLVER_H

#include "source_dftb/sk_matrix.h"

#include <complex>
#include <cstddef>
#include <vector>

namespace ModuleDFTB
{

struct DftbEigenSolution
{
    std::size_t dimension = 0;
    std::vector<double> eigenvalues_hartree;
    // Eigenvectors are column-major: coefficient(orbital, band).
    std::vector<std::complex<double>> eigenvectors;

    const std::complex<double>& coefficient(std::size_t orbital, std::size_t band) const;
};

struct DftbKPointSpectrum
{
    double weight = 0.0;
    DftbEigenSolution solution;
};

struct DftbFermiFilling
{
    double fermi_energy_hartree = 0.0;
    double electron_count = 0.0;
    double band_energy_hartree = 0.0;
    double band_free_energy_hartree = 0.0;
    std::vector<std::vector<double>> occupations;
};

/** Solve H C = S C epsilon for a Hermitian pair of dense complex matrices. */
DftbEigenSolution solve_generalized_hermitian(const DftbBlochMatrices& matrices);

/**
 * Fill spin-degenerate states (0..2 electrons/band) over a normalized k-point
 * mesh. thermal_energy_hartree is k_B*T; zero selects equal filling of
 * degenerate states at the Fermi level.
 */
DftbFermiFilling fill_fermi_occupations(const std::vector<DftbKPointSpectrum>& spectra,
                                         double electron_count,
                                         double thermal_energy_hartree);

} // namespace ModuleDFTB

#endif // SOURCE_DFTB_EIGENSOLVER_H
