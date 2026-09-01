#ifndef STERNHEIMER_RESPONSE_GRID_H
#define STERNHEIMER_RESPONSE_GRID_H

#include <memory>
#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleRI
{

struct SternheimerResponseGrid
{
    SternheimerResponseGrid();
    ~SternheimerResponseGrid();
    SternheimerResponseGrid(SternheimerResponseGrid&&) noexcept;
    SternheimerResponseGrid& operator=(SternheimerResponseGrid&&) noexcept;
    SternheimerResponseGrid(const SternheimerResponseGrid&) = delete;
    SternheimerResponseGrid& operator=(const SternheimerResponseGrid&) = delete;

    std::unique_ptr<ModulePW::PW_Basis> serial_fine_basis;
    std::unique_ptr<ModulePW::PW_Basis> serial_response_basis;
    const ModulePW::PW_Basis* basis = nullptr;
    bool independent = false;
    double requested_ecutwfc = 0.0;
};

bool sternheimer_uses_independent_response_grid(double response_ecutwfc, double pbe_ecutwfc);

SternheimerResponseGrid make_sternheimer_response_grid(const ModulePW::PW_Basis& pbe_basis,
                                                       double pbe_ecutwfc,
                                                       double response_ecutwfc,
                                                       int fft_mode);

std::vector<double> restrict_sternheimer_real_field(const ModulePW::PW_Basis& fine_basis,
                                                    const ModulePW::PW_Basis& coarse_basis,
                                                    const std::vector<double>& fine_values);

} // namespace ModuleRI

#endif
