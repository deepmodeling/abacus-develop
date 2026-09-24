#include "surchem.h"

#include <iomanip>
#include <cmath>
#include <ostream>

double surchem::Acav = 0;
double surchem::Ael = 0;

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
    if (parameters.use_sccs)
    {
        if (parameters.expected_electron_count < 0.0
            || parameters.expected_ionic_charge < 0.0
            || parameters.normalization_tolerance <= 0.0
            || parameters.pool_process_count <= 0
            || !std::isfinite(parameters.start_drho)
            || parameters.start_drho < 0.0
            || parameters.start_nmax <= 0)
        {
            throw std::invalid_argument("SCCS system and reduction parameters are invalid");
        }
    }
    this->parameters_ = parameters;
    this->parameters_set_ = true;
    this->sccs_active_ = parameters.use_sccs && parameters.start_drho <= 0.0;
    this->sccs_state_ = ModuleSccs::SccsState();
    this->sccs_result_ = ModuleSccs::SccsResult();
    this->sccs_elapsed_seconds_ = 0.0;
}

bool surchem::uses_sccs() const
{
    return this->parameters_set_ && this->parameters_.use_sccs;
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
    this->sccs_state_ = ModuleSccs::SccsState();
    this->sccs_result_ = ModuleSccs::SccsResult();
    this->sccs_elapsed_seconds_ = 0.0;
    return true;
}

const ModuleSccs::SccsResult& surchem::sccs_result() const
{
    if (!this->uses_sccs())
    {
        throw std::logic_error("SCCS result requested while the legacy solvent backend is active");
    }
    return this->sccs_result_;
}

void surchem::write_sccs_iteration(std::ostream& output) const
{
    if (!this->sccs_is_active())
    {
        throw std::logic_error("SCCS iteration requested before delayed activation");
    }
    const ModuleSccs::SccsResult& result = this->sccs_result();
    const double solvation_energy_rydberg
        = 2.0 * (result.electrostatic.reaction_energy
                 + result.vacuum_pcc_energy
                 + result.ionic_shape_pcc_energy
                 + result.non_electrostatic.surface_energy
                 + result.non_electrostatic.volume_energy);
    const std::streamsize previous_precision = output.precision();
    const std::ios_base::fmtflags previous_flags = output.flags();
    output << " SCCS_ITER " << result.response.polarization.iterations
           << " SCCS_TIME/s " << std::fixed << std::setprecision(2)
           << this->sccs_elapsed_seconds_
           << " E_SOL/Ry " << std::defaultfloat << std::setprecision(8)
           << solvation_energy_rydberg << '\n';

    if (this->parameters_.debug)
    {
        output << " SCCS_MIXING VALUE "
               << result.response.polarization.final_mixing
               << " RESTARTS "
               << result.response.polarization.mixing_restarts << '\n';
    }

    if (this->parameters_.debug
        && this->parameters_.sccs_config.boundary == ModuleSccs::Boundary::Pcc2d)
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
               << " PCC_ION_SHAPE/Ha " << result.ionic_shape_pcc_energy
               << " PCC_USED/Ry "
               << 2.0 * (result.vacuum_pcc_energy + result.ionic_shape_pcc_energy)
               << '\n';
    }
    output.flags(previous_flags);
    output.precision(previous_precision);
}

void surchem::write_sccs_diagnostics(std::ostream& output) const
{
    if (!this->parameters_.debug)
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
    output << " SCCS_DIAGNOSTIC ionic_shape_pcc_energy_hartree "
           << result.ionic_shape_pcc_energy << '\n';
    output << " SCCS_DIAGNOSTIC electrostatic_energy_rydberg "
           << 2.0 * (result.electrostatic.reaction_energy
                     + result.vacuum_pcc_energy
                     + result.ionic_shape_pcc_energy)
           << '\n';

    if (this->parameters_.sccs_config.boundary == ModuleSccs::Boundary::Pcc2d)
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
    this->sccs_elapsed_seconds_ = 0.0;
}

surchem::~surchem()
{
    this->clear();
}
