#include "read_inp_sccs.h"
#include "read_input.h"
#include "read_input_tool.h"
#include "source_base/tool_quit.h"

namespace ModuleIO
{
void ReadInput::item_sccs()
{
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
        item.description = "Bulk dielectric constant >= 1 for the custom preset. The default 1 (no dielectric) follows the Environ env_static_permittivity default for environ_type input. Water presets use 78.3 and vacuum uses 1.";
        item.default_value = "1.0";
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
    {
        Input_Item item("sccs_lowpass_p1");
        item.annotation = "SCCS switching-derivative low-pass slope";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Low-pass filter of the SCCS switching-function derivatives, as Environ deriv_lowpass_p1 with deriv_method fft. Only with assume_isolated pcc_0d or pcc_2d. The default -1 turns it off and reproduces Environ deriv_method fft (continuum cavity potential). With lowpass disabled (the default), analytical forces may differ from finite differences of the self-consistent energy. For geometry optimization with PCC, consider enabling lowpass and check force accuracy against finite differences. With lowpass disabled, the cavity potential uses the FFT gradient of the PCC-corrected potential, which oscillates around the potential step at the cell boundary half a cell from the system center; keep the dielectric transition region several bohr away from that boundary. Typical values are sccs_lowpass_p1 10 and sccs_lowpass_p2 5. The filter changes the model energy, so compare only results computed with the same values.";
        item.default_value = "-1";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_lowpass_p1);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_lowpass_p2");
        item.annotation = "SCCS switching-derivative low-pass offset";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Offset of the SCCS switching-function low-pass filter, as Environ deriv_lowpass_p2; see sccs_lowpass_p1. Both must be positive or both non-positive. Default -1 (off).";
        item.default_value = "-1";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_lowpass_p2);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_solvent_mode");
        item.annotation = "density that defines the SCCS cavity";
        item.category = "Implicit solvation model";
        item.type = "String";
        item.description = "Allowed values: electronic (default) and full, as Environ solvent_mode. electronic builds the dielectric cavity from the valence electron density. full adds valence-charge Gaussians selected by sccs_corespread.";
        item.default_value = "electronic";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_string(input.sccs_solvent_mode);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_corespread");
        item.annotation = "Widths of the full-cavity Gaussians";
        item.category = "Implicit solvation model";
        item.type = "Vector of Real";
        item.description = "Widths of the valence-charge Gaussians used only with sccs_solvent_mode full. One value applies to every atom except those whose known atomic number equals the pseudopotential valence charge (within 1e-8). More than one value requires exactly nat values in STRU atom order (grouped by atom type), and overrides this automatic exclusion. A value <= 0 disables the cavity Gaussian on that atom. A per-atom list sets symmetry 1 to 0, because the symmetry analysis does not distinguish equivalent atoms by their widths.";
        item.default_value = "0.5";
        item.unit = "bohr";
        item.set_availability("imp_sol==2");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            std::string error;
            const bool valid = parse_spreads(item.str_values, para.input.sccs_corespread, error);
            if (!valid)
            {
                ModuleBase::WARNING_QUIT("ReadInput", error);
            }
        };
        sync_doublevec(input.sccs_corespread, para.input.sccs_corespread.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_solvent_radius");
        item.annotation = "solvent radius of the solvent-aware SCCS cavity";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Solvent radius in bohr of the solvent-aware SCCS cavity, as Environ solvent_radius (Andreussi et al., J. Chem. Theory Comput. 15, 1996 (2019)). A positive value fills the parts of the cavity that a probe sphere of radius sccs_radial_scale * sccs_solvent_radius cannot reach, for example small voids between atoms. Unlike Environ, the probe sums all periodic images within its cutoff, so cells narrower than the probe are not truncated. Non-negative; default 0 (off).";
        item.default_value = "0.0";
        item.unit = "bohr";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_solvent_radius);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_radial_scale");
        item.annotation = "solvent-aware probe radius over sccs_solvent_radius";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Probe radius of the solvent-aware SCCS cavity in units of sccs_solvent_radius, as Environ radial_scale; at least 1, default 2. Used only when sccs_solvent_radius is positive.";
        item.default_value = "2.0";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_radial_scale);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_radial_spread");
        item.annotation = "erfc spread of the solvent-aware probe";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "erfc spread of the solvent-aware probe sphere, as Environ radial_spread; positive, default 0.5 bohr. Used only when sccs_solvent_radius is positive.";
        item.default_value = "0.5";
        item.unit = "bohr";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_radial_spread);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_filling_threshold");
        item.annotation = "probe solute fraction that fills a point";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Solute fraction of the solvent-aware probe sphere above which a point is filled, as Environ filling_threshold; between 0 and 1 (exclusive), default 0.825. Used only when sccs_solvent_radius is positive.";
        item.default_value = "0.825";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_filling_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_filling_spread");
        item.annotation = "erfc width of the solvent-aware filling";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Width of the erfc step of the solvent-aware filling in the probe solute fraction, as Environ filling_spread; positive, default 0.02. Used only when sccs_solvent_radius is positive.";
        item.default_value = "0.02";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_filling_spread);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_start_drho");
        item.annotation = "SCCS delayed-start density threshold";
        item.category = "Implicit solvation model";
        item.type = "Real";
        item.description = "Delays the start of SCCS until DRHO is at or below this value; until then the SCF runs without the solvent, while PCC, if set, stays on. 0 (default) starts SCCS immediately. A positive value must be larger than scf_thr, so that the SCF cannot converge before SCCS starts. When SCCS starts, the charge mixing restarts, that iteration does not count as converged, and the iteration is printed on the screen. SCCS then stays on for the rest of the run, including later ionic steps. Not set by sccs_preset.";
        item.default_value = "0.0";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_double(input.sccs_start_drho);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_start_nmax");
        item.annotation = "SCCS delayed-start iteration limit";
        item.category = "Implicit solvation model";
        item.type = "Integer";
        item.description = "With sccs_start_drho > 0, starts SCCS at this electronic iteration at the latest, if DRHO has not yet reached sccs_start_drho. Must be positive, and smaller than scf_nmax when sccs_start_drho > 0. Default 30. Not set by sccs_preset.";
        item.default_value = "30";
        item.unit = "";
        item.set_availability("imp_sol==2");
        read_sync_int(input.sccs_start_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("sccs_debug");
        item.annotation = "SCCS report level";
        item.category = "Implicit solvation model";
        item.type = "Integer";
        item.description = "Level of the SCCS reports, written to the running log (running_*.log in OUT.{suffix}): 0 (default) writes none; 1 writes the sqrt-CG iteration count and the correction energies of every SCF iteration; 2 also writes the residuals, warm start, cavity volume and surface (with the solvent-aware filled volume), the FFT counts of the Coulomb operator in the sqrt-CG solve and the count of reused fixed sources, and with PCC the Gauss-law check (far-field polarization charge and its expected value), the charge moments and the energy components.";
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
}
} // namespace ModuleIO
