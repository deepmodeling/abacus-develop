#include "read_inp_sccs.h"
#include "read_input.h"
#include "read_input_tool.h"
#include "source_base/tool_quit.h"

#include <algorithm>
#include <cctype>
#include <cmath>

namespace ModuleIO
{
bool parse_solvation_model(const std::string& value, int& model, std::string& error)
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
    if (normalized == "2") { model = 2; }
    else if (std::find(yes.begin(), yes.end(), normalized) != yes.end()) { model = 1; }
    else if (std::find(no.begin(), no.end(), normalized) != no.end()) { model = 0; }
    else
    {
        error = "imp_sol must be 0, 1, 2 or a legacy Boolean value";
        return false;
    }
    return true;
}

bool validate_sccs_input(const Input_para& input, std::string& error)
{
    error.clear();
    if (input.imp_sol < 0 || input.imp_sol > 2)
    {
        error = "imp_sol must be 0, 1 or 2";
        return false;
    }
    if (input.imp_sol != 2) { return true; }
    if (input.assume_isolated != "none") { error = "SCCS currently requires assume_isolated none"; }
    else if (input.device != "cpu" || input.esolver_type != "ksdft") { error = "SCCS requires CPU KS-DFT"; }
    else if (input.basis_type != "pw" && input.basis_type != "lcao") { error = "SCCS requires basis_type pw or lcao"; }
    else if (input.calculation != "scf" || input.cal_stress)
    { error = "SCCS currently supports SCF energies and forces without stress"; }
    else if (input.nspin != 1 && input.nspin != 2) { error = "SCCS requires nspin 1 or 2"; }
    else if (input.efield_flag || input.gate_flag || input.dfthalf_type != 0
             || input.deepks_scf || input.deepks_out_labels || input.deepks_bandgap
             || input.deepks_v_delta || input.dm_to_rho)
    { error = "SCCS does not support electric/gate fields, DFT-1/2, DeePKS or dm_to_rho"; }
    else if (!input.vl_in_h || !input.vion_in_h || !input.vh_in_h)
    { error = "SCCS requires the local ionic and Hartree potentials in the Hamiltonian"; }
    if (!error.empty()) { return false; }
    const std::vector<std::string> presets = {"custom", "vacuum", "water-neutral", "water-cation", "water-anion"};
    if (std::find(presets.begin(), presets.end(), input.sccs_preset) == presets.end())
    {
        error = "Unknown sccs_preset";
        return false;
    }
    if (input.sccs_maxiter <= 0 || !std::isfinite(input.sccs_epsilon) || input.sccs_epsilon < 1.0
        || !std::isfinite(input.sccs_rho_min) || !std::isfinite(input.sccs_rho_max)
        || input.sccs_rho_min <= 0.0 || input.sccs_rho_max <= input.sccs_rho_min
        || !std::isfinite(input.sccs_gamma) || !std::isfinite(input.sccs_pressure)
        || !std::isfinite(input.sccs_tol_rms) || input.sccs_tol_rms <= 0.0
        || !std::isfinite(input.sccs_tol_max) || input.sccs_tol_max <= 0.0
        || !std::isfinite(input.sccs_surface_eta) || input.sccs_surface_eta <= 0.0)
    {
        error = "Invalid SCCS numerical parameters";
        return false;
    }
    return true;
}

void ReadInput::item_sccs()
{
    {
        Input_Item item("imp_sol");
        item.annotation = "implicit solvent model";
        item.category = "Implicit solvation model";
        item.type = "Integer";
        item.description = "Select 0 for vacuum, 1 for the original ABACUS implicit solvation model, or 2 for SCCS. Legacy Boolean values remain accepted as 0 or 1. SCCS currently supports neutral periodic CPU KS-DFT SCF calculations with basis_type pw or lcao and nspin 1 or 2, without stress, external fields or other correction models.";
        item.default_value = "0";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            std::string error;
            const bool valid = parse_solvation_model(item.str_values[0], para.input.imp_sol, error);
            if (!valid) { ModuleBase::WARNING_QUIT("ReadInput", error); }
        };
        sync_int(input.imp_sol);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            std::string error;
            const bool valid = validate_sccs_input(para.input, error);
            if (!valid) { ModuleBase::WARNING_QUIT("ReadInput", error); }
        };
        this->add_item(item);
    }
    {
        Input_Item item("sccs_preset");
        item.annotation = "SCCS parameter preset";
        item.category = "Implicit solvation model";
        item.type = "String";
        item.description = "SCCS parameter preset. Allowed values: custom, vacuum, water-neutral, water-cation, water-anion. Non-custom presets override sccs_epsilon, sccs_rho_min, sccs_rho_max, sccs_gamma and sccs_pressure; solver controls and surface regularization remain user-controlled.";
        item.default_value = "custom";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_string(input.sccs_preset);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_epsilon");
        item.annotation = "Bulk dielectric constant >= 1 for the custom preset";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Bulk dielectric constant >= 1 for the custom preset. Water presets use 78.3 and vacuum uses 1.";
        item.default_value = "78.3";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_epsilon);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_rho_min");
        item.annotation = "Lower electronic-density cavity threshold for the custom preset";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Lower electronic-density cavity threshold for the custom preset. Water-neutral/vacuum: 1e-4; water-cation: 2e-4; water-anion: 2.4e-3.";
        item.default_value = "1e-4";
        item.unit = "e/bohr^3";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_rho_min);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_rho_max");
        item.annotation = "Upper electronic-density cavity threshold for the custom preset";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Upper electronic-density cavity threshold for the custom preset. Water-neutral/vacuum: 5e-3; water-cation: 3.5e-3; water-anion: 1.55e-2.";
        item.default_value = "5e-3";
        item.unit = "e/bohr^3";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_rho_max);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_gamma");
        item.annotation = "Surface coefficient for the custom preset";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Surface coefficient for the custom preset. Water-neutral: 47.9; water-cation: 5; water-anion/vacuum: 0.";
        item.default_value = "0.0";
        item.unit = "dyn/cm";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_gamma);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_pressure");
        item.annotation = "Volume coefficient for the custom preset";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Volume coefficient for the custom preset. Water-neutral: -0.36; water-cation: 0.125; water-anion: 0.45; vacuum: 0.";
        item.default_value = "0.0";
        item.unit = "GPa";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_pressure);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_maxiter");
        item.annotation = "Positive maximum number of SCCS sqrt-CG iterations";
        item.category = "Implicit solvation model";
        item.type = "Integer";
        item.description = "Positive maximum number of SCCS sqrt-CG iterations. Failure to satisfy both residual tolerances terminates the calculation.";
        item.default_value = "200";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_int(input.sccs_maxiter);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_tol_rms");
        item.annotation = "Positive RMS charge-residual tolerance for the SCCS sqrt-CG solver";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Positive RMS charge-residual tolerance for the SCCS sqrt-CG solver.";
        item.default_value = "1e-10";
        item.unit = "e/bohr^3";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_tol_rms);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_tol_max");
        item.annotation = "Positive maximum charge-residual tolerance for the SCCS sqrt-CG solver";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Positive maximum charge-residual tolerance for the SCCS sqrt-CG solver.";
        item.default_value = "1e-8";
        item.unit = "e/bohr^3";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_tol_max);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_surface_eta");
        item.annotation = "Positive regularization of the SCCS surface gradient norm";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Positive regularization of the SCCS surface gradient norm.";
        item.default_value = "1e-8";
        item.unit = "bohr^-1";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_surface_eta);
        this->add_item(item);
    }
}
} // namespace ModuleIO
