#include "dftu_full_dc.h"

#include <cmath>
#include <stdexcept>

namespace DFTU_LCAO
{
template <typename T>
FullDcEnergy evaluate_full_dc(const std::vector<T>& occupation,
                             const int dimension,
                             const double u_eff,
                             const double alpha,
                             std::vector<T>& potential)
{
    if (dimension <= 0 || occupation.size() / dimension != static_cast<std::size_t>(dimension)
        || occupation.size() % dimension != 0)
    {
        throw std::invalid_argument("Full DFT+U requires one square occupation matrix");
    }
    if (!std::isfinite(u_eff) || !std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0)
    {
        throw std::invalid_argument("Full DFT+U requires finite U and alpha in [0,1]");
    }

    double trace = 0.0;
    double square_trace = 0.0;
    for (int i = 0; i < dimension; ++i)
    {
        trace += std::real(occupation[i * dimension + i]);
        for (int j = 0; j < dimension; ++j)
        {
            const T nij = occupation[i * dimension + j];
            const T nji = occupation[j * dimension + i];
            const auto hermiticity_residual = nij - std::conj(nji);
            if (!std::isfinite(std::real(nij)) || !std::isfinite(std::imag(nij))
                || std::abs(hermiticity_residual) > 1.0e-10)
            {
                throw std::invalid_argument("Full DFT+U requires a finite Hermitian occupation matrix");
            }
            square_trace += std::real(nij * nji);
        }
    }

    const double mean = trace / dimension;
    const double shift = 0.5 * alpha + (1.0 - alpha) * mean;
    std::vector<T> result(occupation.size());
    for (int i = 0; i < dimension; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            const double diagonal = (i == j) ? shift : 0.0;
            result[i * dimension + j] = u_eff * (diagonal - occupation[i * dimension + j]);
        }
    }
    const double energy = 0.5 * u_eff * (alpha * trace + (1.0 - alpha) * trace * mean - square_trace);
    const double correction = 0.5 * u_eff * (square_trace - (1.0 - alpha) * trace * mean);
    potential.swap(result);
    return FullDcEnergy{energy, correction};
}

template FullDcEnergy evaluate_full_dc<double>(const std::vector<double>&,
                                              int, double, double, std::vector<double>&);
template FullDcEnergy evaluate_full_dc<std::complex<double>>(
    const std::vector<std::complex<double>>&, int, double, double,
    std::vector<std::complex<double>>&);
} // namespace DFTU_LCAO
