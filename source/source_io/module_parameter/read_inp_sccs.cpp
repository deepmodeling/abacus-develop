#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

#include <cmath>

namespace ModuleIO
{
namespace
{
void check_sccs_preset(const Input_para& input)
{
    const std::vector<std::string> allowed
        = {"custom", "vacuum", "water-neutral", "water-cation", "water-anion"};
    if (std::find(allowed.begin(), allowed.end(), input.sccs_preset) == allowed.end())
    {
        ModuleBase::WARNING_QUIT("ReadInput", nofound_str(allowed, "sccs_preset"));
    }
}

void check_sccs_start_drho(const Input_para& input)
{
    if (!std::isfinite(input.sccs_start_drho)
        || input.sccs_start_drho < 0.0)
    {
        ModuleBase::WARNING_QUIT(
            "ReadInput",
            "sccs_start_drho must be finite and non-negative");
    }
}

void check_sccs_start_nmax(const Input_para& input)
{
    if (input.sccs_start_nmax <= 0
        || (input.sccs_start_drho > 0.0
            && input.sccs_start_nmax >= input.scf_nmax))
    {
        ModuleBase::WARNING_QUIT(
            "ReadInput",
            "sccs_start_nmax must be positive and smaller than scf_nmax when delayed start is enabled");
    }
}

void check_sccs_mixing_type(const Input_para& input)
{
    const std::vector<std::string> allowed = {"linear", "pulay", "anderson"};
    if (std::find(allowed.begin(), allowed.end(), input.sccs_mixing_type) == allowed.end())
    {
        ModuleBase::WARNING_QUIT("ReadInput", nofound_str(allowed, "sccs_mixing_type"));
    }
}

void check_sccs_mixing_parameters(const Input_para& input)
{
    if (input.sccs_mixing_ndim < 2 || !std::isfinite(input.sccs_mixing)
        || input.sccs_mixing <= 0.0 || input.sccs_mixing > 1.0
        || !std::isfinite(input.sccs_mixing_min)
        || !std::isfinite(input.sccs_mixing_max)
        || input.sccs_mixing_min <= 0.0 || input.sccs_mixing_max > 1.0
        || input.sccs_mixing_min > input.sccs_mixing_max
        || (input.sccs_mixing_adaptive
            && (input.sccs_mixing < input.sccs_mixing_min
                || input.sccs_mixing > input.sccs_mixing_max)))
    {
        ModuleBase::WARNING_QUIT("ReadInput", "invalid SCCS mixing parameters");
    }
}

void check_sccs_numerical_parameters(const Input_para& input)
{
    check_sccs_mixing_parameters(input);
    if (input.sccs_maxiter <= 0 || !std::isfinite(input.sccs_epsilon)
        || input.sccs_epsilon < 1.0 || !std::isfinite(input.sccs_rho_min)
        || !std::isfinite(input.sccs_rho_max) || input.sccs_rho_min <= 0.0
        || input.sccs_rho_max <= input.sccs_rho_min
        || !std::isfinite(input.sccs_gamma)
        || !std::isfinite(input.sccs_pressure)
        || !std::isfinite(input.sccs_tol_rms)
        || input.sccs_tol_rms <= 0.0 || !std::isfinite(input.sccs_tol_max)
        || input.sccs_tol_max <= 0.0
        || !std::isfinite(input.sccs_surface_eta)
        || input.sccs_surface_eta <= 0.0)
    {
        ModuleBase::WARNING_QUIT("ReadInput", "invalid SCCS numerical parameters");
    }
}
} // namespace

void ReadInput::item_sccs()
{
    {
        Input_Item item("sccs_preset");
        item.annotation = "SCCS parameter preset";
        item.category = "Implicit solvation model";
        item.type = "String";
        item.description = "Allowed values: custom, vacuum, water-neutral, water-cation, water-anion (use these exact lowercase names). custom uses sccs_epsilon, sccs_rho_min, sccs_rho_max, sccs_gamma and sccs_pressure from INPUT. Other presets override those five values; their individual parameter descriptions list all effective values. The preset is not chosen automatically from the net charge. Solver controls, surface regularization, delayed start and debug settings remain user-controlled for every preset. Vacuum has no dielectric or non-electrostatic solvent contribution, but assume_isolated can still enable PCC.";
        item.default_value = "custom";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_string(input.sccs_preset);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            check_sccs_preset(para.input);
        };
        this->add_item(item);
    }
#define ADD_SCCS_REAL_ITEM(NAME, MEMBER, DESCRIPTION, DEFAULT_VALUE, UNIT) \
    { \
        Input_Item item(NAME); \
        item.annotation = DESCRIPTION; \
        item.category = "Implicit solvation model"; \
        item.type = "Real"; \
        item.description = DESCRIPTION; \
        item.default_value = DEFAULT_VALUE; \
        item.unit = UNIT; \
        item.set_availability("imp_sol==2"); \
        read_sync_double(input.MEMBER); \
        this->add_item(item); \
    }
    ADD_SCCS_REAL_ITEM("sccs_epsilon", sccs_epsilon, "SCCS bulk relative permittivity (dimensionless, at least 1). For sccs_preset=custom, use this INPUT value (default 78.3). Effective value for vacuum: 1; water-neutral, water-cation and water-anion: 78.3. Non-custom presets override this INPUT value.", "78.3", "")
    ADD_SCCS_REAL_ITEM("sccs_rho_min", sccs_rho_min, "SCCS lower cavity-density threshold (positive). For sccs_preset=custom, use this INPUT value (default 1.0e-4 bohr^-3). Effective values: vacuum=1.0e-4, water-neutral=1.0e-4, water-cation=2.0e-4, water-anion=2.4e-3 bohr^-3. Non-custom presets override this INPUT value.", "1.0e-4", "bohr^-3")
    ADD_SCCS_REAL_ITEM("sccs_rho_max", sccs_rho_max, "SCCS upper cavity-density threshold (greater than sccs_rho_min). For sccs_preset=custom, use this INPUT value (default 5.0e-3 bohr^-3). Effective values: vacuum=5.0e-3, water-neutral=5.0e-3, water-cation=3.5e-3, water-anion=1.55e-2 bohr^-3. Non-custom presets override this INPUT value.", "5.0e-3", "bohr^-3")
    ADD_SCCS_REAL_ITEM("sccs_gamma", sccs_gamma, "SCCS effective surface coefficient. For sccs_preset=custom, use this INPUT value (default 0 dyn/cm). Effective values: vacuum=0, water-neutral=47.9, water-cation=5.0, water-anion=0 dyn/cm. Non-custom presets override this INPUT value.", "0.0", "dyn/cm")
    ADD_SCCS_REAL_ITEM("sccs_pressure", sccs_pressure, "SCCS effective volume coefficient. For sccs_preset=custom, use this INPUT value (default 0 GPa). Effective values: vacuum=0, water-neutral=-0.36, water-cation=0.125, water-anion=0.45 GPa. Non-custom presets override this INPUT value.", "0.0", "GPa")
    ADD_SCCS_REAL_ITEM("sccs_mixing", sccs_mixing, "Initial SCCS polarization damping factor for linear, pulay and anderson; range (0, 1]. User-controlled for every sccs_preset, default 0.5. With sccs_mixing_adaptive=1, it must lie within sccs_mixing_min and sccs_mixing_max. Reducing this value can help difficult inner iterations, at the cost of slower convergence.", "0.5", "")
    ADD_SCCS_REAL_ITEM("sccs_mixing_min",
                       sccs_mixing_min,
                       "Lower bound for adaptive SCCS polarization mixing; "
                       "must be positive and no greater than sccs_mixing_max. User-controlled "
                       "for every sccs_preset, default 0.1; inactive when sccs_mixing_adaptive=0",
                       "0.1",
                       "")
    ADD_SCCS_REAL_ITEM("sccs_mixing_max",
                       sccs_mixing_max,
                       "Upper bound for adaptive SCCS polarization mixing; must not exceed one. "
                       "User-controlled for every sccs_preset, default 0.8; inactive when sccs_mixing_adaptive=0",
                       "0.8",
                       "")
    ADD_SCCS_REAL_ITEM("sccs_tol_rms", sccs_tol_rms, "Positive SCCS polarization RMS residual tolerance; user-controlled for every sccs_preset, default 1.0e-10 e/bohr^3.", "1.0e-10", "e/bohr^3")
    ADD_SCCS_REAL_ITEM("sccs_tol_max", sccs_tol_max, "Positive SCCS polarization maximum residual tolerance; user-controlled for every sccs_preset, default 1.0e-8 e/bohr^3.", "1.0e-8", "e/bohr^3")
    ADD_SCCS_REAL_ITEM("sccs_surface_eta", sccs_surface_eta, "Positive SCCS surface regularization; user-controlled for every sccs_preset, default 1.0e-8 bohr^-1.", "1.0e-8", "bohr^-1")
#undef ADD_SCCS_REAL_ITEM
    {
        Input_Item item("sccs_start_drho");
        item.annotation = "SCCS delayed-start density threshold";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Delay SCCS on a cold start until DRHO is at or below this value. Zero starts SCCS immediately. Once activated, SCCS remains active for all later electronic and ionic steps. PCC remains active during the delay. User-controlled for every sccs_preset, default 0.";
        item.default_value = "0.0";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_start_drho);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            check_sccs_start_drho(para.input);
        };
        this->add_item(item);
    }
    {
        Input_Item item("sccs_start_nmax");
        item.annotation = "SCCS delayed-start iteration limit";
        item.category = "Implicit solvation model";
        item.type = "Integer";
        item.description = "Force delayed SCCS activation at this electronic iteration if the SCCS start DRHO threshold has not yet been reached. The value must be positive, and smaller than scf_nmax when delayed start is enabled. User-controlled for every sccs_preset, default 30; inactive when sccs_start_drho=0.";
        item.default_value = "30";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_int(input.sccs_start_nmax);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            check_sccs_start_nmax(para.input);
        };
        this->add_item(item);
    }
    {
        Input_Item item("sccs_debug");
        item.annotation = "detailed SCCS diagnostics";
        item.category = "Implicit solvation model";
        item.type = "Integer";
        item.description = "SCCS/PCC output level: 0 suppresses per-SCF summaries and diagnostics; 1 prints the iteration count, elapsed seconds and correction energy; 2 additionally prints all mixing, multipole and energy diagnostics. Applies to standalone PCC as well as SCCS.";
        item.default_value = "0";
        item.unit = "";
        read_sync_int(input.sccs_debug);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.sccs_debug < 0 || para.input.sccs_debug > 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sccs_debug must be 0, 1, or 2");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("sccs_maxiter");
        item.annotation = "SCCS polarization iteration limit";
        item.category = "Implicit solvation model";
        item.type = "Integer";
        item.description = "Positive maximum number of inner SCCS polarization iterations. User-controlled for every sccs_preset, default 200; failure to converge within this limit terminates the calculation.";
        item.default_value = "200";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_int(input.sccs_maxiter);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            check_sccs_numerical_parameters(para.input);
        };
        this->add_item(item);
    }
    {
        Input_Item item("sccs_mixing_adaptive");
        item.annotation = "adaptive SCCS polarization mixing";
        item.category = "Implicit solvation model";
        item.type = "Boolean";
        item.description
            = "Allowed values: 0 (fixed damping, default) or 1 (adaptive damping). "
              "User-controlled for every sccs_preset. Adapt the damping factor within sccs_mixing_min "
              "and sccs_mixing_max. The initial value is sccs_mixing. Three "
              "consecutive residual ratios below 0.7 increase the factor by "
              "10%. A ratio above 2.0, or two consecutive ratios above 1.1, "
              "halves the factor, clears the acceleration history, and forces "
              "one linear recovery step.";
        item.default_value = "0";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_bool(input.sccs_mixing_adaptive);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            check_sccs_numerical_parameters(para.input);
        };
        this->add_item(item);
    }
    {
        Input_Item item("sccs_mixing_type");
        item.annotation = "SCCS polarization mixing method";
        item.category = "Implicit solvation model";
        item.type = "String";
        item.description
            = "Allowed INPUT values are exactly linear, pulay and anderson. "
              "linear: damped fixed-point iteration; pulay: Pulay DIIS; "
              "anderson: Anderson acceleration using differences of iterates and residuals. "
              "These select the inner SCCS polarization solver, independently of the outer SCF mixing_type. "
              "User-controlled for every sccs_preset, default linear. "
              "broyden and andersonb are not accepted values. Accelerated methods fall back "
              "to a linear step when insufficient or unusable history is available.";
        item.default_value = "linear";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_string(input.sccs_mixing_type);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            check_sccs_mixing_type(para.input);
        };
        this->add_item(item);
    }
    {
        Input_Item item("sccs_mixing_ndim");
        item.annotation = "SCCS accelerated-mixing history length";
        item.category = "Implicit solvation model";
        item.type = "Integer";
        item.description
            = "Number of history vectors retained by sccs_mixing_type=pulay or anderson. "
              "Must be at least 2; user-controlled for every sccs_preset, default 8. "
              "Unused by linear mixing, but the value must still satisfy the input range.";
        item.default_value = "8";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_int(input.sccs_mixing_ndim);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            check_sccs_numerical_parameters(para.input);
        };
        this->add_item(item);
    }
}
} // namespace ModuleIO
