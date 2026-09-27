#ifndef STERNHEIMER_RESPONSE_GRID_H
#define STERNHEIMER_RESPONSE_GRID_H

#include <array>
#include <memory>
#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleRI
{

enum class SternheimerResponseGridSource
{
    Pbe,
    Cutoff,
    Explicit
};

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
    std::array<int, 3> requested_dimensions = {0, 0, 0};
    SternheimerResponseGridSource source = SternheimerResponseGridSource::Pbe;
};

bool sternheimer_uses_independent_response_grid(double response_ecutwfc, double pbe_ecutwfc);

const char* sternheimer_response_grid_source_name(SternheimerResponseGridSource source);

SternheimerResponseGrid make_sternheimer_response_grid(const ModulePW::PW_Basis& pbe_basis,
                                                       double pbe_ecutwfc,
                                                       double response_ecutwfc,
                                                       int fft_mode,
                                                       std::array<int, 3> response_dimensions = {0, 0, 0},
                                                       int fd_order = 2);

std::vector<double> restrict_sternheimer_real_field(const ModulePW::PW_Basis& fine_basis,
                                                    const ModulePW::PW_Basis& coarse_basis,
                                                    const std::vector<double>& fine_values);

std::vector<double> restrict_sternheimer_real_field_rectangular(const ModulePW::PW_Basis& fine_basis,
                                                                const ModulePW::PW_Basis& coarse_basis,
                                                                const std::vector<double>& fine_values);

} // namespace ModuleRI

#endif
