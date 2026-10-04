#include "pot_pcc.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_grid_geometry.h"
#include "source_cell/cell_tools.h"

#include <cmath>
#include <algorithm>

namespace
{
void require_valid_on_pool(const bool valid, const std::string& error, const std::string& function)
{
    double invalid = valid ? 0.0 : 1.0;
#ifdef __MPI
    Parallel_Reduce::reduce_pool(invalid);
#endif
    if (invalid != 0)
    {
        std::string message = error;
        if (message.empty())
        {
            message = "invalid PCC input or result on another process in this pool";
        }
        ModuleBase::WARNING_QUIT(function, message);
    }
}

bool fractional_translation(const ModuleBase::Matrix3& rotation,
                             const ModuleBase::Vector3<double>& translation)
{
    const double tolerance = 1.0e-6;
    const double identity_error = std::abs(rotation.e11 - 1.0) + std::abs(rotation.e22 - 1.0)
                                  + std::abs(rotation.e33 - 1.0) + std::abs(rotation.e12)
                                  + std::abs(rotation.e13) + std::abs(rotation.e21)
                                  + std::abs(rotation.e23) + std::abs(rotation.e31)
                                  + std::abs(rotation.e32);
    if (identity_error > tolerance)
    {
        return false;
    }
    const double dx = translation.x - std::round(translation.x);
    const double dy = translation.y - std::round(translation.y);
    const double dz = translation.z - std::round(translation.z);
    return std::abs(dx) > tolerance || std::abs(dy) > tolerance || std::abs(dz) > tolerance;
}

bool primitive_symmetry(const ModuleSymmetry::Symmetry& symmetry)
{
    for (int operation = 0; operation < symmetry.nrotk; ++operation)
    {
        const bool forbidden = fractional_translation(symmetry.gmatrix[operation], symmetry.gtrans[operation]);
        if (forbidden)
        {
            return false;
        }
    }
    for (int operation = 0; operation < symmetry.nrotk_anti; ++operation)
    {
        const bool forbidden = fractional_translation(symmetry.gmatrix_anti[operation],
                                                       symmetry.gtrans_anti[operation]);
        if (forbidden)
        {
            return false;
        }
    }
    const ModuleBase::Matrix3 identity;
    for (const auto& translation : symmetry.ptrans)
    {
        const bool forbidden = fractional_translation(identity, translation);
        if (forbidden)
        {
            return false;
        }
    }
    return true;
}
bool slab_operation_valid(const ModuleBase::Matrix3& rotation,
                           const ModuleBase::Vector3<double>& translation,
                           const int axis)
{
    const double entries[3][3] = {{rotation.e11, rotation.e12, rotation.e13},
                                  {rotation.e21, rotation.e22, rotation.e23},
                                  {rotation.e31, rotation.e32, rotation.e33}};
    const double translations[3] = {translation.x, translation.y, translation.z};
    for (int other = 0; other < 3; ++other)
    {
        if (other != axis && (std::abs(entries[axis][other]) > 1.0e-6
                              || std::abs(entries[other][axis]) > 1.0e-6))
        {
            return false;
        }
    }
    const double fraction = translations[axis] - std::round(translations[axis]);
    return entries[axis][axis] < 0.0 || std::abs(fraction) <= 1.0e-6;
}

bool slab_symmetry_valid(const ModuleSymmetry::Symmetry& symmetry, const int axis)
{
    for (int operation = 0; operation < symmetry.nrotk; ++operation)
    {
        const bool valid = slab_operation_valid(symmetry.gmatrix[operation], symmetry.gtrans[operation], axis);
        if (!valid) { return false; }
    }
    for (int operation = 0; operation < symmetry.nrotk_anti; ++operation)
    {
        const bool valid = slab_operation_valid(symmetry.gmatrix_anti[operation], symmetry.gtrans_anti[operation], axis);
        if (!valid) { return false; }
    }
    const ModuleBase::Matrix3 identity;
    for (const auto& translation : symmetry.ptrans)
    {
        const bool valid = slab_operation_valid(identity, translation, axis);
        if (!valid) { return false; }
    }
    return true;
}
} // namespace

namespace elecstate
{

PotPcc::PotPcc(const ModulePW::PW_Basis* basis)
    : dimension_(Dimension::molecule), open_axis_(2)
{
    this->rho_basis_ = basis;
    this->dynamic_mode = true;
}

PotPcc::PotPcc(const ModulePW::PW_Basis* basis, const Dimension dimension, const int open_axis)
    : dimension_(dimension), open_axis_(open_axis)
{
    this->rho_basis_ = basis;
    this->dynamic_mode = true;
    if (open_axis < 0 || open_axis > 2)
    {
        ModuleBase::WARNING_QUIT("PotPcc", "pcc_2d_axis must be 0, 1 or 2");
    }
}

void PotPcc::validate_kpoints(const std::vector<ModuleBase::Vector3<double>>& points,
                              const int count,
                              const int open_axis)
{
    double maximum = 0.0;
    if (open_axis < 0 || open_axis > 2 || count < 0 || count > static_cast<int>(points.size()))
    {
        maximum = 1.0;
    }
    else
    {
        for (int ik = 0; ik < count; ++ik)
        {
            const double components[3] = {points[ik].x, points[ik].y, points[ik].z};
            const double coordinate = components[open_axis];
            const double absolute = std::abs(coordinate);
            if (!std::isfinite(coordinate)) { maximum = 1.0; }
            else { maximum = std::max(maximum, absolute); }
        }
    }
#ifdef __MPI
    // All pools must reject a forbidden open-direction k point together.
    Parallel_Reduce::reduce_max(maximum);
#endif
    if (maximum > 1.0e-12)
    {
        ModuleBase::WARNING_QUIT("PotPcc::validate_kpoints", "pcc_2d requires Gamma-only sampling along pcc_2d_axis");
    }
}

ModuleBase::Vector3<double> PotPcc::relative_position(const ModuleBase::Vector3<double>& position) const
{
    if (dimension_ == Dimension::molecule)
    {
        return unitcell::relative_position(position, geometry_);
    }
    const double coordinate = unitcell::relative_coordinate(position, slab_);
    return slab_.normal * coordinate;
}

double PotPcc::correction_energy() const
{
    if (dimension_ == Dimension::molecule) { return pcc_0d_energy(moments_, parameters_); }
    return pcc_2d_energy(moments_, slab_parameters_);
}

double PotPcc::correction_potential(const ModuleBase::Vector3<double>& position) const
{
    if (dimension_ == Dimension::molecule) { return pcc_0d_potential(moments_, position, parameters_); }
    const double coordinate = position * slab_.normal;
    return pcc_2d_potential(moments_, coordinate, slab_.normal, slab_parameters_);
}

ModuleBase::Vector3<double> PotPcc::correction_force(const double charge,
                                                    const ModuleBase::Vector3<double>& position) const
{
    if (dimension_ == Dimension::molecule) { return pcc_0d_force(moments_, charge, position, parameters_); }
    const double coordinate = position * slab_.normal;
    return pcc_2d_force(moments_, charge, coordinate, slab_.normal, slab_parameters_);
}

void PotPcc::prepare_ions(const UnitCell& cell, ChargeMoments& ionic_moments)
{
    ModuleBase::timer::start("PotPcc", "prepare_ions");
    std::string error;
    if (dimension_ == Dimension::molecule)
    {
        const bool geometry_valid = unitcell::make_orthogonal_cell(cell.latvec, cell.lat0, 1.0e-10, geometry_, error);
        require_valid_on_pool(geometry_valid, error, "PotPcc::prepare_ions");
        const bool parameters_valid = make_pcc_0d_parameters(geometry_, 1.0e-10, parameters_, error);
        require_valid_on_pool(parameters_valid, error, "PotPcc::prepare_ions");
        const bool symmetry_valid = primitive_symmetry(cell.symm);
        require_valid_on_pool(symmetry_valid, "PCC 0D requires a primitive cell without fractional-translation symmetry", "PotPcc::prepare_ions");
    }
    else
    {
        const bool geometry_valid = unitcell::make_slab_cell(cell.latvec, cell.lat0, open_axis_, 1.0e-6, slab_, error);
        require_valid_on_pool(geometry_valid, error, "PotPcc::prepare_ions");
        const bool parameters_valid = make_pcc_2d_parameters(slab_, slab_parameters_, error);
        require_valid_on_pool(parameters_valid, error, "PotPcc::prepare_ions");
        const bool symmetry_valid = slab_symmetry_valid(cell.symm, open_axis_);
        require_valid_on_pool(symmetry_valid, "PCC 2D symmetry must preserve the open axis without fractional translations along it", "PotPcc::prepare_ions");
    }

    const std::vector<unitcell::AtomData> atoms = unitcell::get_atom_data(cell.atoms, cell.ntype, cell.lat0);
    const int atom_count = static_cast<int>(atoms.size());
    const bool atom_count_valid = atom_count == cell.nat;
    require_valid_on_pool(atom_count_valid, "PCC atom count does not match UnitCell", "PotPcc::prepare_ions");
    std::vector<double> masses(atom_count);
    ionic_positions_.resize(atom_count);
    ionic_charges_.resize(atom_count);
    for (int atom = 0; atom < atom_count; ++atom)
    {
        ionic_positions_[atom] = atoms[atom].position;
        ionic_charges_[atom] = atoms[atom].valence_charge;
        masses[atom] = atoms[atom].mass;
    }
    if (dimension_ == Dimension::molecule)
    {
        ModuleBase::Vector3<double> center;
        const bool center_valid = unitcell::weighted_center(ionic_positions_, masses, geometry_, center, error);
        require_valid_on_pool(center_valid, error, "PotPcc::prepare_ions");
        geometry_.origin = center;
    }
    else
    {
        double center = 0.0;
        const bool center_valid = unitcell::weighted_center(ionic_positions_, masses, slab_, center, error);
        require_valid_on_pool(center_valid, error, "PotPcc::prepare_ions");
        slab_.origin = center;
    }
    std::vector<ModuleBase::Vector3<double>> relative_ions(atom_count);
    for (int atom = 0; atom < atom_count; ++atom)
    {
        relative_ions[atom] = this->relative_position(ionic_positions_[atom]);
    }
    const double* ionic_charge_data = ionic_charges_.data();
    const ModuleBase::Vector3<double>* ionic_position_data = relative_ions.data();
    const bool ions_valid = charge_moments(ionic_charge_data,
                                            ionic_position_data,
                                            atom_count,
                                            1.0,
                                            ionic_moments,
                                            error);
    require_valid_on_pool(ions_valid, error, "PotPcc::prepare_ions");

    ModuleBase::timer::end("PotPcc", "prepare_ions");
}

ChargeMoments PotPcc::collect_electrons(
    const Charge& charge,
    const UnitCell& cell,
    std::vector<ModuleBase::Vector3<double>>& positions) const
{
    ModuleBase::timer::start("PotPcc", "collect_electrons");
    const ModulePW::PW_Basis& basis = *this->rho_basis_;
    std::string error;
    std::vector<double> electronic_charge(basis.nrxx, 0.0);
    const bool positions_valid = ModulePW::grid_positions(basis, cell.latvec, cell.lat0, positions, error);
    require_valid_on_pool(positions_valid, error, "PotPcc::collect_electrons");
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        for (int spin = 0; spin < charge.nspin; ++spin)
        {
            electronic_charge[ir] -= charge.rho[spin][ir];
        }
        positions[ir] = this->relative_position(positions[ir]);
    }
    ChargeMoments electronic_moments;
    const double volume_element = cell.omega / basis.nxyz;
    const double* electronic_data = electronic_charge.data();
    const ModuleBase::Vector3<double>* position_data = positions.data();
    const bool density_valid = charge_moments(electronic_data,
                                               position_data,
                                               basis.nrxx,
                                               volume_element,
                                               electronic_moments,
                                               error);
    require_valid_on_pool(density_valid, error, "PotPcc::collect_electrons");
    double reduced[5] = {electronic_moments.charge, electronic_moments.dipole.x,
                         electronic_moments.dipole.y, electronic_moments.dipole.z,
                         electronic_moments.second_moment};
#ifdef __MPI
    Parallel_Reduce::reduce_pool(reduced, 5);
#endif
    electronic_moments.charge = reduced[0];
    electronic_moments.dipole = ModuleBase::Vector3<double>(reduced[1], reduced[2], reduced[3]);
    electronic_moments.second_moment = reduced[4];
    ModuleBase::timer::end("PotPcc", "collect_electrons");
    return electronic_moments;
}

void PotPcc::cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential)
{
    ModuleBase::timer::start("PotPcc", "cal_v_eff");
    result_valid_ = false;
    const bool storage_valid = charge != nullptr && cell != nullptr && this->rho_basis_ != nullptr;
    require_valid_on_pool(storage_valid, "PCC requires charge, cell and PW basis storage", "PotPcc::cal_v_eff");
    const ModulePW::PW_Basis& basis = *this->rho_basis_;
    const bool grid_valid = basis.nx > 0 && basis.ny > 0 && basis.nz > 0 && basis.nxyz > 0
                            && basis.nplane >= 0 && basis.nrxx == basis.nx * basis.ny * basis.nplane
                            && (charge->nspin == 1 || charge->nspin == 2)
                            && potential.nr == charge->nspin && potential.nc == basis.nrxx
                            && std::isfinite(cell->omega) && cell->omega > 0.0;
    require_valid_on_pool(grid_valid, "PCC requires an initialized grid and nspin=1/2 potential", "PotPcc::cal_v_eff");
    bool density_storage_valid = charge->rho != nullptr;
    if (density_storage_valid && basis.nrxx > 0)
    {
        for (int spin = 0; spin < charge->nspin; ++spin)
        {
            density_storage_valid = density_storage_valid && charge->rho[spin] != nullptr;
        }
    }
    require_valid_on_pool(density_storage_valid, "PCC charge density is not available", "PotPcc::cal_v_eff");
    ChargeMoments ionic_moments;
    this->prepare_ions(*cell, ionic_moments);
    std::vector<ModuleBase::Vector3<double>> positions;
    const ChargeMoments electronic_moments = this->collect_electrons(*charge, *cell, positions);
    // Ions are replicated on every rank and must be added after the reduction.
    moments_ = add_charge_moments(ionic_moments, electronic_moments);
    energy_rydberg_ = 2.0 * this->correction_energy();
    electron_potential_.resize(basis.nrxx);
    bool result_finite = std::isfinite(energy_rydberg_);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double positive_potential = this->correction_potential(positions[ir]);
        electron_potential_[ir] = -2.0 * positive_potential;
        result_finite = result_finite && std::isfinite(electron_potential_[ir]);
    }
    require_valid_on_pool(result_finite, "PCC energy or potential is nonfinite", "PotPcc::cal_v_eff");
    for (int spin = 0; spin < charge->nspin; ++spin)
    {
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            potential(spin, ir) += electron_potential_[ir];
        }
    }
    result_valid_ = true;
    ModuleBase::timer::end("PotPcc", "cal_v_eff");
}

double PotPcc::get_energy() const
{
    if (!result_valid_)
    {
        ModuleBase::WARNING_QUIT("PotPcc::get_energy", "PCC energy requires an updated potential");
    }
    return energy_rydberg_;
}

const std::vector<double>& PotPcc::electron_potential() const
{
    if (!result_valid_)
    {
        ModuleBase::WARNING_QUIT("PotPcc::electron_potential", "PCC output requires an updated potential");
    }
    return electron_potential_;
}

void PotPcc::add_force(const UnitCell& cell, ModuleBase::matrix& force) const
{
    ModuleBase::timer::start("PotPcc", "add_force");
    if (!result_valid_ || force.nr != cell.nat || force.nc != 3)
    {
        ModuleBase::WARNING_QUIT("PotPcc::add_force", "PCC force requires a current potential and nat-by-3 output");
    }
    const std::vector<unitcell::AtomData> atoms = unitcell::get_atom_data(cell.atoms, cell.ntype, cell.lat0);
    if (atoms.size() != ionic_positions_.size() || atoms.size() != static_cast<std::size_t>(cell.nat))
    {
        ModuleBase::WARNING_QUIT("PotPcc::add_force", "PCC atom count changed since the potential update");
    }
    for (int atom = 0; atom < cell.nat; ++atom)
    {
        const ModuleBase::Vector3<double> difference = atoms[atom].position - ionic_positions_[atom];
        if (difference.norm2() > 1.0e-24 || atoms[atom].valence_charge != ionic_charges_[atom])
        {
            ModuleBase::WARNING_QUIT("PotPcc::add_force", "PCC atoms changed since the potential update");
        }
        const ModuleBase::Vector3<double> relative = this->relative_position(ionic_positions_[atom]);
        const ModuleBase::Vector3<double> correction = this->correction_force(ionic_charges_[atom], relative);
        force(atom, 0) += 2.0 * correction.x;
        force(atom, 1) += 2.0 * correction.y;
        force(atom, 2) += 2.0 * correction.z;
    }
    ModuleBase::timer::end("PotPcc", "add_force");
}

} // namespace elecstate
