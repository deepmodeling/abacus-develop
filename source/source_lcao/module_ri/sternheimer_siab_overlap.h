#ifndef STERNHEIMER_SIAB_OVERLAP_H
#define STERNHEIMER_SIAB_OVERLAP_H

#include <complex>
#include <vector>

namespace module_ri
{
namespace sternheimer_siab
{

double norm(const std::vector<std::complex<double>>& y, double delta_omega);

std::vector<std::complex<double>> overlap_q(
    const std::vector<std::complex<double>>& y,
    const std::vector<std::vector<std::complex<double>>>& primitives,
    double delta_omega);

std::vector<std::complex<double>> overlap_s(
    const std::vector<std::vector<std::complex<double>>>& primitives,
    double delta_omega);

} // namespace sternheimer_siab
} // namespace module_ri

#endif
