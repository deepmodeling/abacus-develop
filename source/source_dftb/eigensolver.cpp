#include "source_dftb/eigensolver.h"

#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace ModuleDFTB
{
namespace
{
double fermi_occupation(double energy, double chemical_potential, double thermal_energy)
{
    const double scaled = (energy - chemical_potential) / thermal_energy;
    if (scaled > 40.0) return 2.0 * std::exp(-scaled);
    if (scaled < -40.0) return 2.0;
    return 2.0 / (1.0 + std::exp(scaled));
}
} // namespace

const std::complex<double>& DftbEigenSolution::coefficient(std::size_t orbital, std::size_t band) const
{
    if (orbital >= dimension || band >= dimension)
    {
        throw std::out_of_range("DFTB eigenvector index is out of range");
    }
    return eigenvectors[band * dimension + orbital];
}

DftbEigenSolution solve_generalized_hermitian(const DftbBlochMatrices& matrices)
{
    if (matrices.dimension == 0
        || matrices.hamiltonian.size() != matrices.dimension * matrices.dimension
        || matrices.overlap.size() != matrices.dimension * matrices.dimension
        || matrices.dimension > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Invalid dimensions for the DFTB generalized eigensystem");
    }
    const int n = static_cast<int>(matrices.dimension);
    const int lda = n;
    const int ldb = n;
    const int itype = 1;
    const char jobz = 'V';
    const char uplo = 'U';
    int info = 0;
    int lwork = -1;
    const std::size_t matrix_size = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
    std::vector<std::complex<double>> h(matrix_size);
    std::vector<std::complex<double>> s(matrix_size);
    for (int row = 0; row < n; ++row)
    {
        for (int col = 0; col < n; ++col)
        {
            h[row + col * n] = matrices.h(static_cast<std::size_t>(row), static_cast<std::size_t>(col));
            s[row + col * n] = matrices.s(static_cast<std::size_t>(row), static_cast<std::size_t>(col));
        }
    }

    std::vector<double> eigenvalues(static_cast<std::size_t>(n));
    std::complex<double> work_query;
    double rwork_query = 0.0;
    zhegv_(&itype, &jobz, &uplo, &n, h.data(), &lda, s.data(), &ldb,
           eigenvalues.data(), &work_query, &lwork, &rwork_query, &info);
    if (info != 0 || !std::isfinite(work_query.real()) || work_query.real() < 1.0)
    {
        throw std::runtime_error("LAPACK zhegv workspace query failed (info=" + std::to_string(info) + ")");
    }
    lwork = std::max(1, static_cast<int>(std::ceil(work_query.real())));
    std::vector<std::complex<double>> work(static_cast<std::size_t>(lwork));
    // ABACUS's LAPACK connector tests also reserve 7*n real work values for
    // vendor implementations whose workspace query is not populated.
    std::vector<double> rwork(std::max<std::size_t>(1, 7 * static_cast<std::size_t>(n)));
    zhegv_(&itype, &jobz, &uplo, &n, h.data(), &lda, s.data(), &ldb,
           eigenvalues.data(), work.data(), &lwork, rwork.data(), &info);
    if (info < 0)
    {
        throw std::runtime_error("LAPACK zhegv rejected argument " + std::to_string(-info));
    }
    if (info > 0)
    {
        throw std::runtime_error("LAPACK zhegv failed: overlap is not positive definite or eigensystem did not converge (info="
                                 + std::to_string(info) + ")");
    }

    DftbEigenSolution result;
    result.dimension = matrices.dimension;
    result.eigenvalues_hartree = std::move(eigenvalues);
    result.eigenvectors = std::move(h);
    return result;
}

DftbFermiFilling fill_fermi_occupations(const std::vector<DftbKPointSpectrum>& spectra,
                                         double electron_count,
                                         double thermal_energy_hartree)
{
    if (spectra.empty() || !std::isfinite(electron_count) || electron_count < 0.0
        || !std::isfinite(thermal_energy_hartree) || thermal_energy_hartree < 0.0)
    {
        throw std::invalid_argument("Invalid k-point spectrum, electron count, or electronic temperature");
    }
    double weight_sum = 0.0;
    double capacity = 0.0;
    double min_energy = std::numeric_limits<double>::infinity();
    double max_energy = -std::numeric_limits<double>::infinity();
    for (const auto& spectrum : spectra)
    {
        if (!(spectrum.weight >= 0.0) || !std::isfinite(spectrum.weight) || spectrum.solution.dimension == 0
            || spectrum.solution.eigenvalues_hartree.size() != spectrum.solution.dimension)
        {
            throw std::invalid_argument("Invalid DFTB k-point spectrum or weight");
        }
        weight_sum += spectrum.weight;
        capacity += 2.0 * spectrum.weight * static_cast<double>(spectrum.solution.dimension);
        for (const double energy : spectrum.solution.eigenvalues_hartree)
        {
            if (!std::isfinite(energy)) throw std::invalid_argument("DFTB eigenvalues must be finite");
            if (spectrum.weight > 0.0)
            {
                min_energy = std::min(min_energy, energy);
                max_energy = std::max(max_energy, energy);
            }
        }
    }
    if (std::abs(weight_sum - 1.0) > 1.0e-10)
    {
        throw std::invalid_argument("DFTB k-point weights must sum to one");
    }
    if (electron_count > capacity + 1.0e-12)
    {
        throw std::invalid_argument("Requested electron count exceeds the number of spin-degenerate bands");
    }

    DftbFermiFilling result;
    result.electron_count = electron_count;
    result.occupations.resize(spectra.size());
    for (std::size_t k = 0; k < spectra.size(); ++k)
    {
        result.occupations[k].resize(spectra[k].solution.dimension);
    }

    bool gap_fermi_found = false;
    double gap_fermi = 0.0;
    if (thermal_energy_hartree > 0.0 && electron_count > 0.0 && electron_count < capacity)
    {
        std::vector<std::pair<double, double>> levels;
        for (const auto& spectrum : spectra)
        {
            if (spectrum.weight == 0.0) continue;
            for (const double energy : spectrum.solution.eigenvalues_hartree)
            {
                levels.emplace_back(energy, 2.0 * spectrum.weight);
            }
        }
        std::sort(levels.begin(), levels.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
            return a.first < b.first;
        });
        double occupied_capacity = 0.0;
        for (std::size_t first = 0; first < levels.size();)
        {
            std::size_t last = first + 1;
            const double tolerance = 1.0e-12 * std::max(1.0, std::abs(levels[first].first));
            double group_capacity = levels[first].second;
            while (last < levels.size() && std::abs(levels[last].first - levels[first].first) <= tolerance)
            {
                group_capacity += levels[last].second;
                ++last;
            }
            occupied_capacity += group_capacity;
            const double electron_tolerance = 1.0e-10 * std::max(1.0, electron_count);
            if (last < levels.size() && std::abs(occupied_capacity - electron_count) <= electron_tolerance
                && levels[last].first - levels[first].first > 80.0 * thermal_energy_hartree)
            {
                gap_fermi = 0.5 * (levels[first].first + levels[last].first);
                gap_fermi_found = true;
                break;
            }
            first = last;
        }
    }

    if (thermal_energy_hartree > 0.0)
    {
        if (electron_count == 0.0)
        {
            result.fermi_energy_hartree = min_energy - 1.0;
        }
        else if (electron_count >= capacity - 1.0e-12)
        {
            result.fermi_energy_hartree = max_energy + 1.0;
        }
        else if (gap_fermi_found)
        {
            // In a wide-gap system, finite-T occupations are exactly clipped to 0/2
            // in double precision over a range of chemical potentials. Use the
            // reproducible mid-gap value instead of a bisection-dependent plateau point.
            result.fermi_energy_hartree = gap_fermi;
        }
        else
        {
            double lower = min_energy - std::max(1.0, 100.0 * thermal_energy_hartree);
            double upper = max_energy + std::max(1.0, 100.0 * thermal_energy_hartree);
            for (int iteration = 0; iteration < 256; ++iteration)
            {
                const double middle = 0.5 * (lower + upper);
                double count = 0.0;
                for (const auto& spectrum : spectra)
                {
                    for (const double energy : spectrum.solution.eigenvalues_hartree)
                    {
                        count += spectrum.weight * fermi_occupation(energy, middle, thermal_energy_hartree);
                    }
                }
                if (count < electron_count) lower = middle;
                else upper = middle;
                if (upper - lower < 1.0e-13) break;
            }
            result.fermi_energy_hartree = 0.5 * (lower + upper);
        }
        for (std::size_t k = 0; k < spectra.size(); ++k)
        {
            if (spectra[k].weight == 0.0) continue;
            for (std::size_t band = 0; band < spectra[k].solution.dimension; ++band)
            {
                const double energy = spectra[k].solution.eigenvalues_hartree[band];
                const double occupation = fermi_occupation(energy, result.fermi_energy_hartree,
                                                           thermal_energy_hartree);
                result.occupations[k][band] = occupation;
                result.band_energy_hartree += spectra[k].weight * occupation * energy;
                const double spin_occupation = 0.5 * occupation;
                if (spin_occupation > 0.0 && spin_occupation < 1.0)
                {
                    result.band_free_energy_hartree += 2.0 * spectra[k].weight * thermal_energy_hartree
                        * (spin_occupation * std::log(spin_occupation)
                           + (1.0 - spin_occupation) * std::log(1.0 - spin_occupation));
                }
            }
        }
        result.band_free_energy_hartree += result.band_energy_hartree;
        return result;
    }

    struct State
    {
        double energy;
        std::size_t k;
        std::size_t band;
        double capacity;
    };
    std::vector<State> states;
    for (std::size_t k = 0; k < spectra.size(); ++k)
    {
        if (spectra[k].weight == 0.0) continue;
        for (std::size_t band = 0; band < spectra[k].solution.dimension; ++band)
        {
            states.push_back({spectra[k].solution.eigenvalues_hartree[band], k, band, 2.0 * spectra[k].weight});
        }
    }
    std::sort(states.begin(), states.end(), [](const State& a, const State& b) { return a.energy < b.energy; });
    double remaining = electron_count;
    std::size_t first = 0;
    while (first < states.size() && remaining > 1.0e-12)
    {
        std::size_t last = first + 1;
        const double degeneracy_tolerance = 1.0e-12 * std::max(1.0, std::abs(states[first].energy));
        while (last < states.size() && std::abs(states[last].energy - states[first].energy) <= degeneracy_tolerance)
        {
            ++last;
        }
        double group_capacity = 0.0;
        for (std::size_t i = first; i < last; ++i) group_capacity += states[i].capacity;
        const double fraction = std::min(1.0, remaining / group_capacity);
        for (std::size_t i = first; i < last; ++i)
        {
            const auto& state = states[i];
            const double occupation = 2.0 * fraction;
            result.occupations[state.k][state.band] = occupation;
            result.band_energy_hartree += spectra[state.k].weight * occupation * state.energy;
        }
        remaining -= fraction * group_capacity;
        result.fermi_energy_hartree = states[first].energy;
        first = last;
    }
    result.band_free_energy_hartree = result.band_energy_hartree;
    return result;
}

} // namespace ModuleDFTB
