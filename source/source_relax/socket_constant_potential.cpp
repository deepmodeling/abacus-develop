#include "socket_constant_potential.h"
#include "source_base/constants.h"
#include "source_esolver/esolver.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_parameter/system_parameter.h"

#include <cmath>
#include <iomanip>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>

namespace
{
std::string trim(const std::string& text)
{
    const auto begin = text.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) { return ""; }
    return text.substr(begin, text.find_last_not_of(" \t\r\n") - begin + 1);
}

double number(const std::map<std::string, std::string>& tokens, const std::string& key)
{
    const auto found = tokens.find(key);
    if (found == tokens.end()) { throw std::runtime_error("constant-potential INIT requires " + key); }
    std::istringstream in(found->second);
    in.imbue(std::locale::classic());
    double value = 0.0;
    if (!(in >> value) || !std::isfinite(value))
    { throw std::runtime_error("INIT " + key + " must be finite and numerical"); }
    in >> std::ws;
    if (!in.eof()) { throw std::runtime_error("invalid trailing data in INIT " + key); }
    return value;
}
}

namespace SocketCP
{
ScopedElectronicStep::ScopedElectronicStep(UnitCell& cell, elecstate::ElectronicStep* step) : cell_(cell)
{
    if (cell_.electronic_step != nullptr)
    { throw std::runtime_error("an electronic step is already bound to this cell"); }
    cell_.electronic_step = step;
}

ScopedElectronicStep::~ScopedElectronicStep()
{
    cell_.electronic_step = nullptr;
}

void run_electronic_step(ModuleESolver::ESolver& solver, UnitCell& cell,
                         const Input_para& inp, const System_para& system,
                         int functional_type, int istep)
{
    if (cell.electronic_step == nullptr || cell.electronic_step->parameters == nullptr
        || cell.electronic_step->request_applied || cell.electronic_step->result_ready)
    { throw std::runtime_error("electron control requires a fresh bound electronic step"); }
    const elecstate::ElectronicControl& control = cell.electronic_step->control;
    if (inp.assume_isolated == "makov-payne")
    {
        throw std::runtime_error(
            "constant-potential electron control is incompatible with the energy-only Makov-Payne correction");
    }
    if (inp.esolver_type != "ksdft" || (inp.basis_type != "pw" && inp.basis_type != "lcao")
        || inp.device != "cpu" || inp.precision != "double" || system.use_uspp
        || inp.sc_mag_switch || inp.deepks_scf || inp.rpa || inp.rdmft || inp.dft_plus_dmft)
    {
        throw std::runtime_error(
            "constant-potential control requires CPU double-precision PW/LCAO KSDFT with "
            "norm-conserving pseudopotentials; specialized solver corrections are not supported");
    }
    if (functional_type < 1 || functional_type > 5)
    {
        throw std::runtime_error(
            "constant-potential control requires a supported KS exchange-correlation functional");
    }
    if (inp.dft_plus_u != 0 && inp.dft_plus_u != 1)
    {
        throw std::runtime_error(
            "constant-potential forces require dft_plus_u=1; "
            "the legacy dft_plus_u=2 force implementation is unavailable");
    }
    if (inp.dft_plus_u && (inp.init_occ_mat < 0 || inp.init_occ_mat > 1
                          || inp.yukawa_potential || system.uramping > 0.0))
    {
        throw std::runtime_error(
            "constant-potential DFT+U requires fixed Hubbard U and self-consistent occupations "
            "(init_occ_mat=0/1, yukawa_potential=false, uramping disabled)");
    }
    if (inp.init_chg == "dm_no_renormalize")
    { throw std::runtime_error("electron-number control requires charge-density renormalization"); }
    if (inp.nspin != 1 && inp.nspin != 2)
    { throw std::runtime_error("constant-potential control supports only nspin=1 or 2"); }
    if (inp.ocp || inp.nbands <= 0 || !std::isfinite(control.nelect) || control.nelect <= 0.0)
    { throw std::runtime_error("invalid total electron number, band capacity, or fixed occupation array"); }
    const bool two_fermi = system.two_fermi;
    if (two_fermi && inp.nspin != 2)
    { throw std::runtime_error("two_fermi electron control requires nspin=2"); }
    const double up = 0.5 * (control.nelect + inp.nupdown);
    const double down = 0.5 * (control.nelect - inp.nupdown);
    if ((two_fermi && (!std::isfinite(inp.nupdown) || up <= 0.0 || down <= 0.0
                      || up > inp.nbands || down > inp.nbands))
        || control.nelect > 2.0 * inp.nbands)
    {
        throw std::runtime_error(
            "electron number exceeds band capacity or gives invalid spin populations; increase nbands");
    }
    if (inp.smearing_method == "fixed"
        && (std::abs(control.nelect - std::round(control.nelect)) > 1.0e-8
            || (two_fermi && (std::abs(up-std::round(up)) > 1.0e-8 || std::abs(down-std::round(down)) > 1.0e-8))))
    { throw std::runtime_error("fractional electron populations require a smearing occupation method"); }
    solver.runner(cell, istep);
    if (!cell.electronic_step->request_applied)
    { throw std::runtime_error("selected ESolver did not apply the constant-potential electron request"); }
}

Configuration parse_init(const std::string& parameters)
{
    const std::string keys[] = {"NELECT", "neutral_electrons",
                               "potential_average_axis", "potential_average_min_A", "potential_average_max_A"};
    std::map<std::string, std::string> tokens;
    std::istringstream input(parameters);
    std::string token;
    while (std::getline(input, token, ','))
    {
        token = trim(token);
        const auto colon = token.find(':');
        const std::string key = trim(token.substr(0, colon));
        if (key == "solvation_update_stride")
        { throw std::runtime_error("INIT solvation_update_stride has been removed"); }
        for (const auto& known : keys)
        {
            if (key != known) { continue; }
            if (colon == std::string::npos || tokens.count(key))
            { throw std::runtime_error("missing separator or duplicate constant-potential INIT token: " + key); }
            tokens[key] = trim(token.substr(colon + 1));
        }
    }
    Configuration config;
    if (tokens.empty()) { return config; }
    config.enabled = true;
    if (number(tokens, "NELECT") <= 0.0)
    { throw std::runtime_error("INIT NELECT must be positive"); }
    if (number(tokens, "neutral_electrons") <= 0.0)
    { throw std::runtime_error("INIT neutral_electrons must be positive"); }
    if (tokens.count("potential_average_axis") || tokens.count("potential_average_min_A")
        || tokens.count("potential_average_max_A"))
    {
        const std::string axis = tokens["potential_average_axis"];
        if (axis.size() != 1 || axis[0] < 'x' || axis[0] > 'z')
        { throw std::runtime_error("INIT potential_average_axis must be x, y or z"); }
        config.window.axis = axis[0] - 'x';
        config.window.minimum_A = number(tokens, "potential_average_min_A");
        config.window.maximum_A = number(tokens, "potential_average_max_A");
        if (config.window.minimum_A < 0.0 || config.window.maximum_A <= config.window.minimum_A)
        { throw std::runtime_error("INIT potential window must satisfy 0 <= min < max"); }
    }
    return config;
}

elecstate::ElectronicControl make_control(const Configuration& config, double nelect, int flags)
{
    if (!config.enabled) { throw std::runtime_error("CHGDATA requires constant-potential INIT metadata"); }
    if (!std::isfinite(nelect) || nelect <= 0.0)
    { throw std::runtime_error("CHGDATA total electron number must be finite and positive"); }
    if (flags != 0 && flags != 2)
    { throw std::runtime_error("CHGDATA flags must be 0 or 2; bit 0 is reserved"); }
    if (flags & 2) { throw std::runtime_error("ABACUS Ne-electrode CHGDATA is not supported"); }
    elecstate::ElectronicControl control;
    control.nelect = nelect;
    return control;
}

Observables collect_observables(const Configuration& config, const elecstate::ElectronicStep& step)
{
    if (!config.enabled) { throw std::runtime_error("electronic observables require constant-potential INIT"); }
    if (!step.request_applied || !step.result_ready)
    { throw std::runtime_error("ESolver did not complete this frame's electronic observables"); }
    const elecstate::ElectronicObservables& physical = step.result;
    Observables result;
    result.fermi_eV = physical.fermi_ry * ModuleBase::Ry_to_eV;
    result.nelect = physical.nelect;
    if (!std::isfinite(result.fermi_eV) || !std::isfinite(result.nelect) || result.nelect <= 0.0)
    { throw std::runtime_error("invalid constant-potential electronic observables"); }
    result.has_workfunction = config.window.axis >= 0;
    if (result.has_workfunction)
    {
        result.reference_eV = physical.reference_ry * ModuleBase::Ry_to_eV;
        result.workfunction_eV = result.reference_eV - result.fermi_eV;
        if (!std::isfinite(result.reference_eV) || !std::isfinite(result.workfunction_eV))
        { throw std::runtime_error("non-finite workfunction"); }
    }
    return result;
}

std::string serialize_observables(const Observables& result)
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(17) << "{\"fermi_level_eV\":" << result.fermi_eV
        << ",\"nelect\":" << result.nelect;
    if (result.has_workfunction)
    {
        out << ",\"workfunction_eV\":" << result.workfunction_eV
            << ",\"electrostatic_reference_eV\":" << result.reference_eV;
    }
    out << "}";
    return out.str();
}
}
