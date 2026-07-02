#include "source_lcao/module_ri/sternheimer_abacus_fd_smoke.h"

#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_io/module_parameter/parameter.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

namespace ModuleRI
{
namespace
{

constexpr const char* kSmokeEnv = "ABACUS_STERNHEIMER_FD_ZERO_ORDER_SMOKE";
constexpr const char* kOutputEnv = "ABACUS_STERNHEIMER_FD_ZERO_ORDER_OUT";
constexpr const char* kBandsEnv = "ABACUS_STERNHEIMER_FD_ZERO_ORDER_BANDS";
constexpr const char* kMaxDenseEnv = "ABACUS_STERNHEIMER_FD_ZERO_ORDER_MAX_DENSE";
constexpr const char* kToleranceEnv = "ABACUS_STERNHEIMER_FD_ZERO_ORDER_TOLERANCE";

std::string lower_string(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool env_is_true(const char* name)
{
    const char* raw = std::getenv(name);
    if (raw == nullptr)
    {
        return false;
    }
    const std::string value = lower_string(raw);
    return !(value.empty() || value == "0" || value == "false" || value == "off" || value == "no");
}

int positive_int_from_env(const char* name, const int default_value)
{
    const char* raw = std::getenv(name);
    if (raw == nullptr)
    {
        return default_value;
    }
    std::size_t parsed = 0;
    const int value = std::stoi(raw, &parsed);
    if (raw[parsed] != '\0' || value <= 0)
    {
        throw std::invalid_argument(std::string("Invalid positive integer in ") + name + ".");
    }
    return value;
}

double nonnegative_double_from_env(const char* name, const double default_value)
{
    const char* raw = std::getenv(name);
    if (raw == nullptr)
    {
        return default_value;
    }
    std::size_t parsed = 0;
    const double value = std::stod(raw, &parsed);
    if (raw[parsed] != '\0' || value < 0.0)
    {
        throw std::invalid_argument(std::string("Invalid non-negative floating-point value in ") + name + ".");
    }
    return value;
}

std::string default_report_path(const std::string& output_dir)
{
    const char* explicit_path = std::getenv(kOutputEnv);
    if (explicit_path != nullptr && explicit_path[0] != '\0')
    {
        return explicit_path;
    }
    if (output_dir.empty())
    {
        return "STERNHEIMER_FD_ZERO_ORDER.dat";
    }
    if (output_dir.back() == '/')
    {
        return output_dir + "STERNHEIMER_FD_ZERO_ORDER.dat";
    }
    return output_dir + "/STERNHEIMER_FD_ZERO_ORDER.dat";
}

void write_failure_report(std::ofstream& out, const std::string& reason)
{
    out << "# ABACUS Sternheimer FD zero-order smoke test\n";
    out << "status failed\n";
    out << "reason " << reason << '\n';
}

} // namespace

bool sternheimer_fd_zero_order_smoke_enabled()
{
    return env_is_true(kSmokeEnv);
}

void run_sternheimer_fd_zero_order_smoke(const elecstate::Potential& potential,
                                         const ModulePW::PW_Basis& pw_basis,
                                         const UnitCell& ucell,
                                         const elecstate::ElecState& elec_state,
                                         const std::string& output_dir)
{
    if (!sternheimer_fd_zero_order_smoke_enabled() || GlobalV::MY_RANK != 0)
    {
        return;
    }

    const std::string report_path = default_report_path(output_dir);
    std::ofstream out(report_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out)
    {
        GlobalV::ofs_running << " Sternheimer FD zero-order smoke: failed to open " << report_path << std::endl;
        return;
    }

    try
    {
        if (GlobalV::NPROC != 1)
        {
            throw std::runtime_error(
                "The current dense FD smoke test requires a single MPI rank so pw_basis.nrxx is the full grid.");
        }
        if (elec_state.ekb.nc <= 0)
        {
            throw std::runtime_error("ABACUS DFT eigenvalues are not available.");
        }

        const int requested_bands = positive_int_from_env(kBandsEnv, 4);
        const int num_bands = std::min(requested_bands, elec_state.ekb.nc);
        const int max_dense_size = positive_int_from_env(kMaxDenseEnv, 4096);
        const double eigenvalue_tolerance = nonnegative_double_from_env(kToleranceEnv, 1.0e-2);

        const SternheimerABACUSDFTZeroOrderResult result = compare_sternheimer_abacus_fd_zero_order_to_dft(
            potential,
            pw_basis,
            ucell,
            elec_state,
            0,
            0,
            num_bands,
            eigenvalue_tolerance,
            max_dense_size,
            1.0);

        out << format_sternheimer_fd_zero_order_report(result);
        GlobalV::ofs_running << " Sternheimer FD zero-order smoke report: " << report_path << std::endl;
    }
    catch (const std::exception& error)
    {
        write_failure_report(out, error.what());
        GlobalV::ofs_running << " Sternheimer FD zero-order smoke failed: " << error.what() << std::endl;
        GlobalV::ofs_running << " Sternheimer FD zero-order smoke report: " << report_path << std::endl;
    }
}

} // namespace ModuleRI
