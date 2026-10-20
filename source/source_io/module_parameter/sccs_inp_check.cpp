// Parsing and validation of the SCCS INPUT keys; the keys themselves are
// registered in read_inp_sccs.cpp.
#include "read_inp_sccs.h"
#include "input_parameter.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <utility>

namespace ModuleIO
{
bool parse_imp_sol(const std::string& value, int& model, std::string& error)
{
    error.clear();
    std::string normalized = value;
    for (char& c : normalized)
    {
        const unsigned char byte = static_cast<unsigned char>(c);
        const int lower = std::tolower(byte);
        c = static_cast<char>(lower);
    }
    const std::vector<std::string> yes = {"true", "1", "t", "yes", "y", "on", ".true."};
    const std::vector<std::string> no = {"false", "0", "f", "no", "n", "off", ".false."};
    if (normalized == "2")
    {
        model = 2;
    }
    else if (std::find(yes.begin(), yes.end(), normalized) != yes.end())
    {
        model = 1;
    }
    else if (std::find(no.begin(), no.end(), normalized) != no.end())
    {
        model = 0;
    }
    else
    {
        error = "imp_sol must be 0, 1, 2 or a legacy Boolean value";
        return false;
    }
    return true;
}

bool parse_spreads(const std::vector<std::string>& words, std::vector<double>& spreads, std::string& error)
{
    error.clear();
    std::vector<double> values;
    for (const std::string& word : words)
    {
        const char* begin = word.c_str();
        char* end = nullptr;
        const double value = std::strtod(begin, &end);
        const bool complete = end != begin && *end == '\0';
        if (!complete || !std::isfinite(value))
        {
            error = "sccs_corespread values must be finite real numbers, got '" + word + "'";
            return false;
        }
        values.push_back(value);
    }
    spreads.swap(values);
    return true;
}

namespace
{
// Calculation settings that SCCS supports; empty when all are supported.
std::string context_error(const Input_para& input)
{
    std::string error;
    if (input.assume_isolated != "none" && input.assume_isolated != "pcc_0d"
        && input.assume_isolated != "pcc_2d")
    {
        error = "SCCS currently supports assume_isolated none, pcc_0d or pcc_2d";
    }
    else if (input.device != "cpu" || input.esolver_type != "ksdft")
    {
        error = "SCCS requires CPU KS-DFT";
    }
    else if (input.basis_type != "pw" && input.basis_type != "lcao")
    {
        error = "SCCS requires basis_type pw or lcao";
    }
    else if ((input.calculation != "scf" && input.calculation != "relax") || input.cal_stress)
    {
        error = "SCCS supports calculation scf or fixed-cell relax without stress";
    }
    else if (input.nspin != 1 && input.nspin != 2)
    {
        error = "SCCS requires nspin 1 or 2";
    }
    else if (input.efield_flag || input.gate_flag || input.dfthalf_type != 0
             || input.deepks_scf || input.deepks_out_labels || input.deepks_bandgap
             || input.deepks_v_delta || input.dm_to_rho)
    {
        error = "SCCS does not support electric/gate fields, DFT-1/2, DeePKS or dm_to_rho";
    }
    else if (!input.vl_in_h || !input.vion_in_h || !input.vh_in_h)
    {
        error = "SCCS requires the local ionic and Hartree potentials in the Hamiltonian";
    }
    return error;
}

// nan and inf pass the ordered comparisons below. The INPUT reader already
// rejects them; this also covers an Input_para filled in code. Signed values
// stay allowed for gamma/pressure and for the lowpass pair, where non-positive
// values turn the filter off.
std::string nonfinite_error(const Input_para& input)
{
    const std::vector<std::pair<const char*, double>> values = {
        {"sccs_epsilon", input.sccs_epsilon},
        {"sccs_rho_min", input.sccs_rho_min},
        {"sccs_rho_max", input.sccs_rho_max},
        {"sccs_gamma", input.sccs_gamma},
        {"sccs_pressure", input.sccs_pressure},
        {"sccs_tol_rms", input.sccs_tol_rms},
        {"sccs_tol_max", input.sccs_tol_max},
        {"sccs_surface_eta", input.sccs_surface_eta},
        {"sccs_lowpass_p1", input.sccs_lowpass_p1},
        {"sccs_lowpass_p2", input.sccs_lowpass_p2},
        {"sccs_start_drho", input.sccs_start_drho}};
    for (const std::pair<const char*, double>& value : values)
    {
        if (!std::isfinite(value.second))
        {
            return std::string(value.first) + " must be finite";
        }
    }
    return std::string();
}

// Cavity and surface model parameters, starting with the preset.
std::string model_error(const Input_para& input)
{
    const std::vector<std::string> presets = {"custom", "vacuum", "water-neutral", "water-cation", "water-anion"};
    if (std::find(presets.begin(), presets.end(), input.sccs_preset) == presets.end())
    {
        return "Unknown sccs_preset";
    }
    const bool p1_positive = input.sccs_lowpass_p1 > 0.0;
    const bool p2_positive = input.sccs_lowpass_p2 > 0.0;
    const bool pcc = input.assume_isolated == "pcc_0d" || input.assume_isolated == "pcc_2d";
    std::string error;
    if (input.sccs_epsilon < 1.0)
    {
        error = "sccs_epsilon must be at least 1";
    }
    else if (input.sccs_rho_min <= 0.0 || input.sccs_rho_max <= input.sccs_rho_min)
    {
        error = "sccs_rho_min and sccs_rho_max must satisfy 0 < sccs_rho_min < sccs_rho_max";
    }
    else if (input.sccs_surface_eta <= 0.0)
    {
        error = "sccs_surface_eta must be positive";
    }
    else if (p1_positive != p2_positive)
    {
        error = "sccs_lowpass_p1 and sccs_lowpass_p2 must be both positive or both non-positive";
    }
    else if (p1_positive && !pcc)
    {
        error = "sccs_lowpass_p1 and sccs_lowpass_p2 require assume_isolated pcc_0d or pcc_2d";
    }
    else if (input.sccs_solvent_mode != "electronic" && input.sccs_solvent_mode != "full")
    {
        error = "sccs_solvent_mode must be electronic or full";
    }
    else if (input.sccs_corespread.empty())
    {
        error = "sccs_corespread must contain one value or exactly nat values";
    }
    return error;
}

// sqrt-CG solver controls and the delayed start.
std::string solver_error(const Input_para& input)
{
    const bool delayed_start = input.sccs_start_drho > 0.0;
    std::string error;
    if (input.sccs_maxiter <= 0)
    {
        error = "sccs_maxiter must be positive";
    }
    else if (input.sccs_tol_rms <= 0.0)
    {
        error = "sccs_tol_rms must be positive";
    }
    else if (input.sccs_tol_max <= 0.0)
    {
        error = "sccs_tol_max must be positive";
    }
    else if (input.sccs_start_drho < 0.0)
    {
        error = "sccs_start_drho must be non-negative";
    }
    else if (delayed_start && input.sccs_start_drho <= input.scf_thr)
    {
        error = "sccs_start_drho must be zero or larger than scf_thr";
    }
    else if (input.sccs_start_nmax <= 0 || (delayed_start && input.sccs_start_nmax >= input.scf_nmax))
    {
        error = "sccs_start_nmax must be positive and smaller than scf_nmax when delayed start is enabled";
    }
    return error;
}
} // namespace

bool check_sccs_inp(const Input_para& input, std::string& error)
{
    error.clear();
    if (input.imp_sol < 0 || input.imp_sol > 2)
    {
        error = "imp_sol must be 0, 1 or 2";
        return false;
    }
    if (input.imp_sol != 2)
    {
        return true;
    }
    error = context_error(input);
    if (error.empty())
    {
        error = nonfinite_error(input);
    }
    if (error.empty())
    {
        error = model_error(input);
    }
    if (error.empty())
    {
        error = solver_error(input);
    }
    return error.empty();
}
} // namespace ModuleIO
