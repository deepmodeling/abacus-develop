#include "vdwd4.h"

#include "source_base/constants.h"
#include "source_base/element_name.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstring>
#include <string>
#include <vector>

#ifdef __DFTD4
#include <dftd4.h>
#endif

namespace vdw
{

namespace
{

std::string to_lower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

void check_dftd4_error(dftd4_error error, const std::string& where)
{
    if (dftd4_check_error(error))
    {
        char buffer[1024];
        std::memset(buffer, 0, sizeof(buffer));
        int buffersize = static_cast<int>(sizeof(buffer));
        dftd4_get_error(error, buffer, &buffersize);

        ModuleBase::WARNING_QUIT("Vdwd4::" + where, std::string(buffer));
    }
}

int atomic_number_from_symbol(const std::string& symbol)
{
    for (int i = 0; i < static_cast<int>(ModuleBase::element_name.size()); ++i)
    {
        if (symbol == ModuleBase::element_name[i])
        {
            return i + 1; // DFT-D4 expects the true atomic number.
        }
    }

    ModuleBase::WARNING_QUIT("Vdwd4::atomic_number_from_symbol", "Unknown element symbol: " + symbol);
    return 0;
}

double length_to_bohr(double value, const std::string& unit)
{
    if (unit == "A")
    {
        return value / ModuleBase::BOHR_TO_A;
    }
    if (unit == "Bohr")
    {
        return value;
    }

    ModuleBase::WARNING_QUIT("Vdwd4::length_to_bohr", "Unsupported length unit: " + unit);
    return value;
}

double cutoff_to_bohr(const std::string& value, const std::string& unit)
{
    return length_to_bohr(std::stod(value), unit);
}

} // namespace

Vdwd4::Vdwd4(const UnitCell& unit_in, const std::string& xc_name, const Input_para& input)
    : Vdw(unit_in), xc_name_(to_lower(xc_name)), model_name_(to_lower(input.vdw_d4_model))
{
    cutoff_disp2_ = cutoff_to_bohr(input.vdw_cutoff_radius, input.vdw_radius_unit);
    cutoff_disp3_ = std::min(40.0, cutoff_disp2_);
    cutoff_cn_ = length_to_bohr(input.vdw_cn_thr, input.vdw_cn_thr_unit);
}

void Vdwd4::build_structure(std::vector<int>& numbers,
                            std::vector<double>& positions,
                            std::vector<double>& lattice,
                            std::array<bool, 3>& periodic) const
{
    numbers.clear();
    positions.clear();
    lattice.clear();

    numbers.reserve(ucell_.nat);
    positions.reserve(3 * ucell_.nat);
    lattice.reserve(9);

    for (int it = 0; it < ucell_.ntype; ++it)
    {
        const int atomic_number = atomic_number_from_symbol(ucell_.atoms[it].ncpp.psd);

        for (int ia = 0; ia < ucell_.atoms[it].na; ++ia)
        {
            const ModuleBase::Vector3<double> position = ucell_.atoms[it].tau[ia] * ucell_.lat0;

            numbers.push_back(atomic_number);
            positions.push_back(position.x);
            positions.push_back(position.y);
            positions.push_back(position.z);
        }
    }

    const ModuleBase::Vector3<double> a1 = ucell_.a1 * ucell_.lat0;
    const ModuleBase::Vector3<double> a2 = ucell_.a2 * ucell_.lat0;
    const ModuleBase::Vector3<double> a3 = ucell_.a3 * ucell_.lat0;

    // DFT-D4 C API documents lattice as lattice[3][3] in Bohr.
    // Keep this order aligned with positions[natoms][3], and verify it with
    // a non-orthogonal-cell regression test before relying on stress values.
    lattice.push_back(a1.x);
    lattice.push_back(a1.y);
    lattice.push_back(a1.z);
    lattice.push_back(a2.x);
    lattice.push_back(a2.y);
    lattice.push_back(a2.z);
    lattice.push_back(a3.x);
    lattice.push_back(a3.y);
    lattice.push_back(a3.z);

    periodic = {{true, true, true}};
}

void Vdwd4::compute(double& energy_ha,
                    std::vector<double>* gradient_ha_bohr,
                    std::array<double, 9>* sigma_ha)
{
#ifdef __DFTD4
    std::vector<int> numbers;
    std::vector<double> positions;
    std::vector<double> lattice;
    std::array<bool, 3> periodic;

    build_structure(numbers, positions, lattice, periodic);

    if (gradient_ha_bohr != nullptr
        && gradient_ha_bohr->size() != static_cast<std::size_t>(3 * ucell_.nat))
    {
        ModuleBase::WARNING_QUIT("Vdwd4::compute",
                                 "gradient_ha_bohr must have size 3 * nat when requested.");
    }

    // These vectors own all arrays passed to DFT-D4. Their data() pointers
    // remain valid until all DFT-D4 handles created below are deleted.
    dftd4_error error = dftd4_new_error();

    dftd4_structure mol = dftd4_new_structure(error,
                                              ucell_.nat,
                                              numbers.data(),
                                              positions.data(),
                                              nullptr,
                                              lattice.data(),
                                              periodic.data());
    check_dftd4_error(error, "dftd4_new_structure");

    dftd4_model model = nullptr;
    if (model_name_ == "d4")
    {
        model = dftd4_new_d4_model(error, mol);
        check_dftd4_error(error, "dftd4_new_d4_model");
    }
    else if (model_name_ == "d4s")
    {
        model = dftd4_new_d4s_model(error, mol);
        check_dftd4_error(error, "dftd4_new_d4s_model");
    }
    else
    {
        ModuleBase::WARNING_QUIT("Vdwd4::compute", "Unsupported DFT-D4 model: " + model_name_);
    }

    dftd4_set_model_realspace_cutoff(error, model, cutoff_disp2_, cutoff_disp3_, cutoff_cn_);
    check_dftd4_error(error, "dftd4_set_model_realspace_cutoff");

    std::vector<char> method(xc_name_.begin(), xc_name_.end());
    method.push_back('\0');

    // Use the DFT-D4 library's internal rational damping parameters.
    // The final boolean selects the three-body-specific parameterization in the C API.
    const bool atm = true;
    dftd4_param param = dftd4_load_rational_damping(error, method.data(), atm);
    check_dftd4_error(error, "dftd4_load_rational_damping");

    // These are borrowed output buffers for this synchronous C API call only.
    double* gradient = gradient_ha_bohr ? gradient_ha_bohr->data() : nullptr;
    double* sigma = sigma_ha ? sigma_ha->data() : nullptr;

    dftd4_get_dispersion(error, mol, model, param, &energy_ha, gradient, sigma);
    check_dftd4_error(error, "dftd4_get_dispersion");

    dftd4_delete_param(&param);
    dftd4_delete_model(&model);
    dftd4_delete_structure(&mol);
    dftd4_delete_error(&error);
#else
    ModuleBase::WARNING_QUIT("Vdwd4::compute", "DFT-D4 correction requires ABACUS to be built with DFTD4.");
#endif
}

void Vdwd4::cal_energy()
{
    ModuleBase::TITLE("Vdwd4", "cal_energy");
    ModuleBase::timer::start("Vdwd4", "cal_energy");

    double energy_ha = 0.0;
    compute(energy_ha, nullptr, nullptr);

    // DFT-D4 returns Hartree; ABACUS vdW energies are stored in Ry.
    energy_ = 2.0 * energy_ha;

    ModuleBase::timer::end("Vdwd4", "cal_energy");
}

void Vdwd4::set_force_from_gradient(const std::vector<double>& gradient_ha_bohr)
{
    force_.clear();
    force_.resize(ucell_.nat);

    for (int iat = 0; iat < ucell_.nat; ++iat)
    {
        // DFT-D4 returns dE/dR in Ha/Bohr; ABACUS forces are -dE/dR in Ry/Bohr.
        force_[iat].x = -2.0 * gradient_ha_bohr[3 * iat + 0];
        force_[iat].y = -2.0 * gradient_ha_bohr[3 * iat + 1];
        force_[iat].z = -2.0 * gradient_ha_bohr[3 * iat + 2];
    }

    has_force_cache_ = true;
}

void Vdwd4::set_stress_from_sigma(const std::array<double, 9>& sigma_ha)
{
    // Tentative mapping consistent with the current D3 convention.
    // Confirm sign, transposition and volume normalization by finite-strain tests.
    stress_ = ModuleBase::Matrix3(2.0 * sigma_ha[0], 2.0 * sigma_ha[1], 2.0 * sigma_ha[2],
                                  2.0 * sigma_ha[3], 2.0 * sigma_ha[4], 2.0 * sigma_ha[5],
                                  2.0 * sigma_ha[6], 2.0 * sigma_ha[7], 2.0 * sigma_ha[8])
              / ucell_.omega;

    has_stress_cache_ = true;
}

void Vdwd4::cal_force()
{
    ModuleBase::TITLE("Vdwd4", "cal_force");
    ModuleBase::timer::start("Vdwd4", "cal_force");

    if (!has_force_cache_ || !has_stress_cache_)
    {
        double energy_ha = 0.0;
        std::vector<double> gradient(3 * ucell_.nat, 0.0);
        std::array<double, 9> sigma;
        sigma.fill(0.0);

        // Request sigma together with the gradient.  The DFT-D4 C API computes
        // sigma internally for gradient calculations anyway, so keep it and
        // avoid a second expensive D4 call when ABACUS subsequently requests stress.
        compute(energy_ha, &gradient, &sigma);
        set_force_from_gradient(gradient);
        set_stress_from_sigma(sigma);
    }

    ModuleBase::timer::end("Vdwd4", "cal_force");
}

void Vdwd4::cal_stress()
{
    ModuleBase::TITLE("Vdwd4", "cal_stress");
    ModuleBase::timer::start("Vdwd4", "cal_stress");

    if (!has_stress_cache_)
    {
        double energy_ha = 0.0;
        std::vector<double> gradient(3 * ucell_.nat, 0.0);
        std::array<double, 9> sigma;
        sigma.fill(0.0);

        // DFT-D4 may require a valid gradient buffer when sigma is requested.
        compute(energy_ha, &gradient, &sigma);
        set_force_from_gradient(gradient);
        set_stress_from_sigma(sigma);
    }

    ModuleBase::timer::end("Vdwd4", "cal_stress");
}

} // namespace vdw
