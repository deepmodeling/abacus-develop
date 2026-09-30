#include "surchem.h"

#include <iomanip>
#include <cmath>
#include <ostream>

double surchem::Acav = 0;
double surchem::Ael = 0;
double surchem::Epcc = 0;

surchem::surchem()
{
    TOTN_real = nullptr;
    delta_phi = nullptr;
    epspot = nullptr;
    Vcav = ModuleBase::matrix();
    Vel = ModuleBase::matrix();
    qs = 0;
}

void surchem::set_parameters(const SurchemParameters& parameters)
{
    if (parameters.use_sccs || parameters.pcc_boundary != ModulePcc::Boundary::Periodic)
    {
        if (parameters.expected_electron_count < 0.0
            || parameters.expected_ionic_charge < 0.0
            || parameters.normalization_tolerance <= 0.0
            || parameters.pool_process_count <= 0
            || !std::isfinite(parameters.start_drho)
            || parameters.start_drho < 0.0
            || parameters.start_nmax <= 0)
        {
            throw std::invalid_argument("SCCS/PCC system and reduction parameters are invalid");
        }
    }
    if (parameters.debug < 0 || parameters.debug > 2
        || (parameters.use_legacy_solvent && (parameters.use_sccs
            || parameters.pcc_boundary != ModulePcc::Boundary::Periodic))
        || (parameters.use_sccs && parameters.sccs_config.boundary != parameters.pcc_boundary))
    {
        throw std::invalid_argument("inconsistent solvent/PCC configuration or debug level");
    }
    this->parameters_ = parameters;
    this->fixed_source_cache_ = FixedSourceCache();
    this->parameters_set_ = true;
    this->sccs_active_ = parameters.use_sccs && parameters.start_drho <= 0.0;
    this->sccs_state_ = ModuleSccs::SccsState();
    this->sccs_result_ = ModuleSccs::SccsResult();
    this->pcc_result_valid_ = false;
    this->electrostatic_correction_ry_.clear();
}

bool surchem::uses_sccs() const
{
    return this->parameters_set_ && this->parameters_.use_sccs;
}

bool surchem::uses_pcc() const
{
    return this->parameters_set_
           && this->parameters_.pcc_boundary != ModulePcc::Boundary::Periodic;
}

bool surchem::sccs_is_active() const
{
    return this->uses_sccs() && this->sccs_active_;
}

bool surchem::try_activate_sccs(const int electronic_iteration, const double drho)
{
    if (!this->uses_sccs() || this->sccs_active_)
    {
        return false;
    }
    if (drho > this->parameters_.start_drho
        && electronic_iteration < this->parameters_.start_nmax)
    {
        return false;
    }
    this->sccs_active_ = true;
    this->fixed_source_cache_ = FixedSourceCache();
    this->sccs_state_ = ModuleSccs::SccsState();
    this->sccs_result_ = ModuleSccs::SccsResult();
    this->pcc_result_valid_ = false;
    return true;
}

const std::vector<double>& surchem::electrostatic_correction() const
{
    return this->electrostatic_correction_ry_;
}

const ModuleSccs::SccsResult& surchem::sccs_result() const
{
    if (!this->uses_sccs())
    {
        throw std::logic_error("SCCS result requested while the legacy solvent backend is active");
    }
    return this->sccs_result_;
}

void surchem::write_iteration(std::ostream& output, const double drho) const
{
    if (this->parameters_.debug == 0)
    {
        return;
    }
    if (this->sccs_is_active() || this->uses_pcc())
    {
        this->write_sccs_iteration(output);
    }
    else if (this->uses_sccs())
    {
        output << " SCCS_DEFERRED DRHO " << drho
               << " START_DRHO " << this->parameters_.start_drho
               << " START_NMAX " << this->parameters_.start_nmax << std::endl;
    }
}

void surchem::write_sccs_iteration(std::ostream& output) const
{
    if (this->parameters_.debug == 0)
    {
        return;
    }
    if (!this->sccs_is_active())
    {
        if (!this->uses_pcc() || !this->pcc_result_valid_)
        {
            throw std::logic_error("correction summary requires a current SCCS or PCC result");
        }
        const std::streamsize precision = output.precision();
        const std::ios_base::fmtflags flags = output.flags();
        output << " E_PCC/Ry " << std::setprecision(12) << this->pcc_energy_rydberg_ << '\n';
        if (this->parameters_.debug >= 2)
        {
            if (this->parameters_.pcc_boundary == ModulePcc::Boundary::Pcc0d)
            {
                output << " PCC0D_ORIGIN X/Bohr " << this->pcc_geometry_.origin.x
                       << " Y/Bohr " << this->pcc_geometry_.origin.y
                       << " Z/Bohr " << this->pcc_geometry_.origin.z << '\n'
                       << " PCC0D_MOMENTS POINT Q/e " << this->pcc_moments_.charge
                       << " PX/eBohr " << this->pcc_moments_.dipole.x
                       << " PY/eBohr " << this->pcc_moments_.dipole.y
                       << " PZ/eBohr " << this->pcc_moments_.dipole.z
                       << " QRR/eBohr2 " << this->pcc_moments_.quadrupole_trace << '\n';
            }
            else
            {
                output << " PCC2D_ORIGIN Y/Bohr " << this->pcc_2d_geometry_.origin_y << '\n'
                       << " PCC2D_MOMENTS Q_POINT/e " << this->pcc_2d_moments_.charge
                       << " PY_POINT/eBohr " << this->pcc_2d_moments_.dipole_y
                       << " QYY_POINT/eBohr2 " << this->pcc_2d_moments_.quadrupole_yy << '\n';
            }
            output << " PCC_ENERGY PCC_POINT/Ha " << 0.5 * this->pcc_energy_rydberg_
                   << " PCC_USED/Ry " << this->pcc_energy_rydberg_ << '\n';
        }
        output.flags(flags);
        output.precision(precision);
        return;
    }
    const ModuleSccs::SccsResult& result = this->sccs_result();
    const double solvation_energy_rydberg
        = 2.0 * (result.electrostatic.reaction_energy
                 + result.non_electrostatic.surface_energy
                 + result.non_electrostatic.volume_energy);
    const std::streamsize previous_precision = output.precision();
    const std::ios_base::fmtflags previous_flags = output.flags();
    output << " SCCS_ITER " << result.response.polarization.iterations
           << " E_SOL/Ry " << std::setprecision(8) << solvation_energy_rydberg;
    if (this->uses_pcc())
    {
        output << " E_PCC/Ry " << 2.0 * result.vacuum_pcc_energy;
    }
    output << '\n';

    if (this->parameters_.debug >= 2)
    {
        const ModuleSccs::PolarizationResult& polarization = result.response.polarization;
        output << " SCCS_RESIDUAL RMS " << polarization.residual_rms
               << " MAX " << polarization.residual_max
               << " WARM_START " << polarization.warm_started << '\n';
        if (polarization.fixed_point_checked)
        {
            output << " SCCS_CG_FIXED_POINT_DEFECT RMS " << polarization.fixed_point_defect_rms
                   << " MAX " << polarization.fixed_point_defect_max << '\n';
        }
        if (this->parameters_.sccs_config.boundary != ModulePcc::Boundary::Periodic)
        {
            // Gauss's law: FAR_FIELD is the solution's net screening charge,
            // DENSITY the dielectric_of_potential integral with its grid error.
            const bool slab = this->parameters_.sccs_config.boundary == ModulePcc::Boundary::Pcc2d;
            const double solute_charge
                = slab ? result.solute_moments_2d.charge : result.solute_moments.charge;
            const double density_charge
                = slab ? result.polarization_moments_2d.charge : result.polarization_moments.charge;
            const double epsilon_bulk = this->parameters_.sccs_config.cavity.epsilon_bulk;
            const double expected_charge = -(1.0 - 1.0 / epsilon_bulk) * solute_charge;
            output << " SCCS_GAUSS Q_POL_FAR_FIELD/e "
                   << result.response.far_field_polarization_charge
                   << " Q_POL_DENSITY/e " << density_charge
                   << " Q_POL_EXPECTED/e " << expected_charge << '\n';
        }
        // Local output-rank FFT counts of the sqrt-CG solve; timings are in
        // the ModuleBase::timer summary.
        const ModuleSccs::CoulombTransformCounts& counts = result.forward_transforms;
        output << " SCCS_FFT R2G_CALLS " << counts.forward_calls
               << " G2R_CALLS " << counts.inverse_calls
               << " CACHED_SOURCES " << result.reused_fixed_sources << '\n';
    }

    if (this->parameters_.debug >= 2
        && this->parameters_.sccs_config.boundary == ModulePcc::Boundary::Pcc0d)
    {
        const ModuleBase::Vector3<double>& origin = this->sccs_state_.pcc_geometry.origin;
        output << std::setprecision(12)
               << " PCC0D_ORIGIN X/Bohr " << origin.x
               << " Y/Bohr " << origin.y
               << " Z/Bohr " << origin.z << '\n';
        const ModulePcc::MultipoleMoments* moments[] = {
            &result.solute_moments,
            &result.point_solute_moments,
            &result.polarization_moments,
            &result.screened_moments};
        const char* labels[] = {"SMOOTH", "POINT", "POLARIZATION", "SCREENED"};
        for (int index = 0; index < 4; ++index)
        {
            output << " PCC0D_MOMENTS " << labels[index]
                   << " Q/e " << moments[index]->charge
                   << " PX/eBohr " << moments[index]->dipole.x
                   << " PY/eBohr " << moments[index]->dipole.y
                   << " PZ/eBohr " << moments[index]->dipole.z
                   << " QRR/eBohr2 " << moments[index]->quadrupole_trace
                   << '\n';
        }
        output << " PCC0D_ENERGY"
               << " REACTION/Ha " << result.electrostatic.reaction_energy
               << " PCC_SMOOTH/Ha " << result.smooth_vacuum_pcc_energy
               << " PCC_POINT/Ha " << result.vacuum_pcc_energy
               << " PCC_USED/Ry " << 2.0 * result.vacuum_pcc_energy
               << '\n';
    }

    if (this->parameters_.debug >= 2
        && this->parameters_.sccs_config.boundary == ModulePcc::Boundary::Pcc2d)
    {
        output << std::setprecision(12)
               << " PCC2D_MOMENTS"
               << " Q_SMOOTH/e " << result.solute_moments_2d.charge
               << " PY_SMOOTH/eBohr " << result.solute_moments_2d.dipole_y
               << " QYY_SMOOTH/eBohr2 " << result.solute_moments_2d.quadrupole_yy
               << " Q_POINT/e " << result.point_solute_moments_2d.charge
               << " PY_POINT/eBohr " << result.point_solute_moments_2d.dipole_y
               << " QYY_POINT/eBohr2 " << result.point_solute_moments_2d.quadrupole_yy
               << '\n';
        output << " PCC2D_ENERGY"
               << " REACTION/Ha " << result.electrostatic.reaction_energy
               << " PCC_SMOOTH/Ha " << result.smooth_vacuum_pcc_energy
               << " PCC_POINT/Ha " << result.vacuum_pcc_energy
               << " PCC_USED/Ry " << 2.0 * result.vacuum_pcc_energy
               << '\n';
    }
    output.flags(previous_flags);
    output.precision(previous_precision);
}

void surchem::write_sccs_diagnostics(std::ostream& output) const
{
    if (this->parameters_.debug < 2)
    {
        return;
    }
    const ModuleSccs::SccsResult& result = this->sccs_result();
    const std::streamsize previous_precision = output.precision();
    output << std::setprecision(16);
    output << " SCCS_DIAGNOSTIC reaction_energy_hartree "
           << result.electrostatic.reaction_energy << '\n';
    output << " SCCS_DIAGNOSTIC smooth_vacuum_pcc_energy_hartree "
           << result.smooth_vacuum_pcc_energy << '\n';
    output << " SCCS_DIAGNOSTIC point_vacuum_pcc_energy_hartree "
           << result.vacuum_pcc_energy << '\n';
    output << " SCCS_DIAGNOSTIC electrostatic_energy_rydberg "
           << 2.0 * (result.electrostatic.reaction_energy
                     + result.vacuum_pcc_energy)
           << '\n';

    if (this->parameters_.sccs_config.boundary == ModulePcc::Boundary::Pcc2d)
    {
        output << " SCCS_DIAGNOSTIC smooth_solute_charge "
               << result.solute_moments_2d.charge << '\n';
        output << " SCCS_DIAGNOSTIC smooth_solute_dipole_y "
               << result.solute_moments_2d.dipole_y << '\n';
        output << " SCCS_DIAGNOSTIC smooth_solute_quadrupole_yy "
               << result.solute_moments_2d.quadrupole_yy << '\n';
        output << " SCCS_DIAGNOSTIC point_solute_charge "
               << result.point_solute_moments_2d.charge << '\n';
        output << " SCCS_DIAGNOSTIC point_solute_dipole_y "
               << result.point_solute_moments_2d.dipole_y << '\n';
        output << " SCCS_DIAGNOSTIC point_solute_quadrupole_yy "
               << result.point_solute_moments_2d.quadrupole_yy << '\n';
        output << " SCCS_DIAGNOSTIC polarization_charge "
               << result.polarization_moments_2d.charge << '\n';
        output << " SCCS_DIAGNOSTIC far_field_polarization_charge "
               << result.response.far_field_polarization_charge << '\n';
        output << " SCCS_DIAGNOSTIC polarization_dipole_y "
               << result.polarization_moments_2d.dipole_y << '\n';
        output << " SCCS_DIAGNOSTIC polarization_quadrupole_yy "
               << result.polarization_moments_2d.quadrupole_yy << '\n';
        output << " SCCS_DIAGNOSTIC screened_charge "
               << result.screened_moments_2d.charge << '\n';
        output << " SCCS_DIAGNOSTIC screened_dipole_y "
               << result.screened_moments_2d.dipole_y << '\n';
        output << " SCCS_DIAGNOSTIC screened_quadrupole_yy "
               << result.screened_moments_2d.quadrupole_yy << '\n';
    }
    output.precision(previous_precision);
}

void surchem::allocate(const int &nrxx, const int &nspin)
{
    assert(nrxx >= 0);
    assert(nspin > 0);

    delete[] TOTN_real;
    delete[] delta_phi;
    delete[] epspot;
    if (nrxx > 0)
    {
        TOTN_real = new double[nrxx];
        delta_phi = new double[nrxx];
        epspot = new double[nrxx];
    }
    else
    {
        TOTN_real = nullptr;
        delta_phi = nullptr;
        epspot = nullptr;
    }
    Vcav.create(nspin, nrxx);
    Vel.create(nspin, nrxx);

    ModuleBase::GlobalFunc::ZEROS(delta_phi, nrxx);
    ModuleBase::GlobalFunc::ZEROS(TOTN_real, nrxx);
    ModuleBase::GlobalFunc::ZEROS(epspot, nrxx);
    return;
}

void surchem::clear()
{
    this->pcc_result_valid_ = false;
    delete[] TOTN_real;
    delete[] delta_phi;
    delete[] epspot;
    this->TOTN_real = nullptr;
    this->delta_phi = nullptr;
    this->epspot = nullptr;

    this->Vcav.create(0, 0); 
    this->Vel.create(0, 0);
    this->sccs_state_ = ModuleSccs::SccsState();
    this->sccs_result_ = ModuleSccs::SccsResult();
    this->fixed_source_cache_ = FixedSourceCache();
    this->electrostatic_correction_ry_.clear();
}

surchem::~surchem()
{
    this->clear();
}
