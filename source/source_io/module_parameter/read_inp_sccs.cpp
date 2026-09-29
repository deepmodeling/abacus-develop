#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

#include <cmath>
#include <iostream>

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
    // A threshold at or below scf_thr lets the SCF converge before SCCS starts.
    if (input.sccs_start_drho > 0.0 && input.sccs_start_drho <= input.scf_thr)
    {
        ModuleBase::WARNING_QUIT(
            "ReadInput",
            "sccs_start_drho must be zero or larger than scf_thr");
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

// The sqrt-preconditioned CG replaced the polarization fixed point, so the
// sccs_mixing* keys are read for old inputs but ignored.
const char* const deprecated_sccs_mixing_description
    = "Deprecated and ignored. SCCS is solved with the ENVIRON sqrt-preconditioned CG for "
      "every assume_isolated value (none, pcc_0d, pcc_2d), which uses no polarization "
      "mixing. The key is still read so that old INPUT files run; setting it only prints a "
      "warning. It is independent of the electronic SCF mixing_type and mixing_beta.";

// check_value runs on rank 0 before warning.log is opened, so report on stdout.
void warn_deprecated_sccs_mixing(const Input_Item& item)
{
    if (item.is_read())
    {
        std::cout << " WARNING: " << item.label
                  << " is deprecated and ignored; SCCS always uses the sqrt-preconditioned CG"
                  << std::endl;
    }
}

void check_sccs_lowpass(const Input_para& input)
{
    const bool p1_positive = input.sccs_lowpass_p1 > 0.0;
    const bool p2_positive = input.sccs_lowpass_p2 > 0.0;
    if (!std::isfinite(input.sccs_lowpass_p1) || !std::isfinite(input.sccs_lowpass_p2)
        || p1_positive != p2_positive)
    {
        ModuleBase::WARNING_QUIT("ReadInput",
                                 "sccs_lowpass_p1 and sccs_lowpass_p2 must be finite and both positive or both non-positive");
    }
    if (p1_positive && input.assume_isolated != "pcc_0d" && input.assume_isolated != "pcc_2d")
    {
        ModuleBase::WARNING_QUIT("ReadInput",
                                 "sccs_lowpass_p1 and sccs_lowpass_p2 require assume_isolated pcc_0d or pcc_2d");
    }
}

void check_sccs_solvent_mode(const Input_para& input)
{
    const std::vector<std::string> allowed = {"electronic", "full"};
    if (std::find(allowed.begin(), allowed.end(), input.sccs_solvent_mode) == allowed.end())
    {
        ModuleBase::WARNING_QUIT("ReadInput", nofound_str(allowed, "sccs_solvent_mode"));
    }
    if (!std::isfinite(input.sccs_corespread) || input.sccs_corespread <= 0.0)
    {
        ModuleBase::WARNING_QUIT("ReadInput", "sccs_corespread must be positive and finite");
    }
}

void check_sccs_numerical_parameters(const Input_para& input)
{
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
    check_sccs_lowpass(input);
    check_sccs_solvent_mode(input);
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
    ADD_SCCS_REAL_ITEM("sccs_tol_rms", sccs_tol_rms, "Positive RMS tolerance of the SCCS inner charge residual in e/bohr^3; user-controlled for every sccs_preset, default 1.0e-10. It applies to the ENVIRON sqrt-preconditioned CG solution of the generalized Poisson equation for every assume_isolated value; with pcc_0d or pcc_2d the preconditioner Poisson solve includes the analytic open-boundary correction. ENVIRON stops its CG when the unnormalized sum of squared residuals falls below its tol; the corresponding RMS is sqrt(tol/N) for N FFT grid points.", "1.0e-10", "e/bohr^3")
    ADD_SCCS_REAL_ITEM("sccs_tol_max", sccs_tol_max, "Positive maximum tolerance of the SCCS inner charge residual in e/bohr^3; user-controlled for every sccs_preset, default 1.0e-8. The sqrt-CG stops only when both sccs_tol_rms and sccs_tol_max are satisfied.", "1.0e-8", "e/bohr^3")
    ADD_SCCS_REAL_ITEM("sccs_surface_eta", sccs_surface_eta, "Positive SCCS surface regularization; user-controlled for every sccs_preset, default 1.0e-8 bohr^-1.", "1.0e-8", "bohr^-1")
    ADD_SCCS_REAL_ITEM("sccs_corespread",
                       sccs_corespread,
                       "Spread of the core-electron Gaussians of sccs_solvent_mode full, as Environ "
                       "corespread: exp(-r^2/spread^2), positive, default 0.5 bohr. Used only "
                       "with sccs_solvent_mode full.",
                       "0.5",
                       "bohr")
    ADD_SCCS_REAL_ITEM("sccs_lowpass_p1",
                       sccs_lowpass_p1,
                       "Low-pass filter of the SCCS switching-function derivatives, as Environ "
                       "deriv_lowpass_p1 with deriv_method fft: when sccs_lowpass_p1 and "
                       "sccs_lowpass_p2 are both positive, every Fourier derivative of the "
                       "switching function is multiplied by 0.5 erfc(p1 G^2/Gcut^2 - p2), Gcut^2 "
                       "being the ecutrho sphere, and the electronic potential becomes the exact "
                       "derivative of the discrete SCCS energy, so forces agree with energy "
                       "differences. Only with assume_isolated pcc_0d or pcc_2d. The default -1 "
                       "turns it off and reproduces Environ deriv_method fft (continuum cavity "
                       "potential). 10 with sccs_lowpass_p2 5 was validated at ecutrho 300-500 Ry; "
                       "the filter changes the model energy (about 10 meV for H3O+).",
                       "-1",
                       "")
    ADD_SCCS_REAL_ITEM("sccs_lowpass_p2",
                       sccs_lowpass_p2,
                       "Offset of the SCCS switching-function low-pass filter, as Environ "
                       "deriv_lowpass_p2; see sccs_lowpass_p1. Both must be positive or both "
                       "non-positive. Default -1 (off).",
                       "-1",
                       "")
#undef ADD_SCCS_REAL_ITEM
    {
        Input_Item item("sccs_start_drho");
        item.annotation = "SCCS delayed-start density threshold";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Delay SCCS on a cold start until DRHO is at or below this value. Zero starts SCCS immediately; a positive value must exceed scf_thr so that the SCF cannot converge before SCCS starts. Once activated, SCCS remains active for all later electronic and ionic steps. PCC remains active during the delay. User-controlled for every sccs_preset, default 0.";
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
        item.description = "SCCS/PCC output level: 0 suppresses per-SCF summaries and diagnostics; 1 prints the iteration count, elapsed seconds and correction energy; 2 additionally prints residual, warm-start, timing, Gauss-law (PCC), multipole and energy diagnostics, and verifies the sqrt-CG fixed point with one extra Poisson solve per SCCS evaluation. Applies to standalone PCC as well as SCCS.";
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
        item.annotation = "SCCS sqrt-CG iteration limit";
        item.category = "Implicit solvation model";
        item.type = "Integer";
        item.description = "Positive maximum inner iteration count; user-controlled for every sccs_preset, default 200. SCCS_ITER counts sqrt-preconditioned CG iterations for every assume_isolated value. Failure to converge terminates the calculation. No discrete adjoint is solved.";
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
        Input_Item item("sccs_solvent_mode");
        item.annotation = "density that defines the SCCS cavity";
        item.category = "Implicit solvation model";
        item.type = "String";
        item.description
            = "Allowed values: electronic (default) and full, as Environ solvent_mode. electronic "
              "builds the dielectric cavity from the valence electron density. full adds, on "
              "every atom except hydrogen, a Gaussian of the valence charge with spread "
              "sccs_corespread, as Environ does for the core electrons. Use full when the "
              "pseudo-valence density at a nucleus drops below sccs_rho_max (for example "
              "some S and Cl norm-conserving pseudopotentials with the water-anion or "
              "water-cation preset): electronic mode then puts dielectric inside the atom and "
              "the SCF diverges. The Gaussians shape only the cavity, not the solute charge; "
              "the ionic forces include their cavity term. The published SCCS presets were "
              "fitted with electronic mode.";
        item.default_value = "electronic";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_string(input.sccs_solvent_mode);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            check_sccs_solvent_mode(para.input);
        };
        this->add_item(item);
    }
#define ADD_DEPRECATED_SCCS_MIXING_ITEM(NAME, TYPE, DEFAULT_VALUE, READ) \
    { \
        Input_Item item(NAME); \
        item.annotation = "deprecated, ignored"; \
        item.category = "Implicit solvation model"; \
        item.type = TYPE; \
        item.description = deprecated_sccs_mixing_description; \
        item.default_value = DEFAULT_VALUE; \
        item.unit = ""; \
        item.set_availability("imp_sol==2"); \
        READ; \
        item.check_value = [](const Input_Item& self, const Parameter&) { \
            warn_deprecated_sccs_mixing(self); \
        }; \
        this->add_item(item); \
    }
    ADD_DEPRECATED_SCCS_MIXING_ITEM("sccs_mixing", "Real", "0.5", read_sync_double(input.sccs_mixing))
    ADD_DEPRECATED_SCCS_MIXING_ITEM("sccs_mixing_min", "Real", "0.1", read_sync_double(input.sccs_mixing_min))
    ADD_DEPRECATED_SCCS_MIXING_ITEM("sccs_mixing_max", "Real", "0.8", read_sync_double(input.sccs_mixing_max))
    ADD_DEPRECATED_SCCS_MIXING_ITEM("sccs_mixing_adaptive", "Boolean", "0", read_sync_bool(input.sccs_mixing_adaptive))
    ADD_DEPRECATED_SCCS_MIXING_ITEM("sccs_mixing_type", "String", "linear", read_sync_string(input.sccs_mixing_type))
    ADD_DEPRECATED_SCCS_MIXING_ITEM("sccs_mixing_ndim", "Integer", "8", read_sync_int(input.sccs_mixing_ndim))
#undef ADD_DEPRECATED_SCCS_MIXING_ITEM
}
} // namespace ModuleIO
