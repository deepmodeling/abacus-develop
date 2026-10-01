#include "surchem_input.h"
#include "surchem.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_cell/klist.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace ModuleSurchem
{
namespace
{

const double symmetry_tolerance = 1.0e-6;

double matrix_element(const ModuleBase::Matrix3& matrix, const int row, const int column)
{
    const double elements[3][3] = {{matrix.e11, matrix.e12, matrix.e13},
                                   {matrix.e21, matrix.e22, matrix.e23},
                                   {matrix.e31, matrix.e32, matrix.e33}};
    return elements[row][column];
}

double vector_component(const ModuleBase::Vector3<double>& vector, const int index)
{
    return index == 0 ? vector.x : (index == 1 ? vector.y : vector.z);
}

bool is_lattice_translation(const double fraction)
{
    const double nearest_integer = std::round(fraction);
    const double translation_error = fraction - nearest_integer;
    return std::abs(translation_error) < symmetry_tolerance;
}

bool is_identity(const ModuleBase::Matrix3& rotation)
{
    for (int row = 0; row < 3; ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            const double expected = row == column ? 1.0 : 0.0;
            const double element = matrix_element(rotation, row, column);
            const double identity_error = element - expected;
            if (std::abs(identity_error) > symmetry_tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

bool is_pure_translation(const ModuleBase::Matrix3& rotation, const ModuleBase::Vector3<double>& translation)
{
    const bool lattice_translation = is_lattice_translation(translation.x)
                                     && is_lattice_translation(translation.y)
                                     && is_lattice_translation(translation.z);
    return is_identity(rotation) && !lattice_translation;
}

// Operations act on direct coordinates as r' = r G + t. The open axis maps
// onto +-itself exactly when row and column `axis` of G vanish off the
// diagonal; with G(axis, axis) = +1 a fractional t along it would copy the
// slab within the cell.
std::string pcc_2d_operation_violation(const ModuleBase::Matrix3& rotation,
                                       const ModuleBase::Vector3<double>& translation,
                                       const int axis)
{
    for (int other = 0; other < 3; ++other)
    {
        if (other == axis)
        {
            continue;
        }
        const double row_element = matrix_element(rotation, axis, other);
        const double column_element = matrix_element(rotation, other, axis);
        if (std::abs(row_element) > symmetry_tolerance || std::abs(column_element) > symmetry_tolerance)
        {
            return "a symmetry operation mixes the pcc_2d open lattice vector with a periodic one";
        }
    }
    const double diagonal = matrix_element(rotation, axis, axis);
    const double translation_along_axis = vector_component(translation, axis);
    if (diagonal > 0.0 && !is_lattice_translation(translation_along_axis))
    {
        return "a symmetry operation translates the structure by a fraction of the pcc_2d open lattice vector";
    }
    return std::string();
}

std::string operation_violation(const SurchemParameters& parameters,
                                const ModulePcc::Boundary boundary,
                                const ModuleBase::Matrix3& rotation,
                                const ModuleBase::Vector3<double>& translation)
{
    if (boundary == ModulePcc::Boundary::Pcc2d)
    {
        return pcc_2d_operation_violation(rotation, translation, parameters.pcc_2d_axis);
    }
    if (is_pure_translation(rotation, translation))
    {
        return "a symmetry operation is a fractional translation, so the pcc_0d cell is not primitive";
    }
    return std::string();
}

} // namespace

SurchemParameters make_parameters(const Input_para& inp,
                                  const UnitCell& ucell,
                                  const double electron_count,
                                  const int pool_process_count)
{
    SurchemParameters parameters;
    parameters.eb_k = inp.eb_k;
    parameters.tau = inp.tau;
    parameters.sigma_k = inp.sigma_k;
    parameters.nc_k = inp.nc_k;
    parameters.use_sccs = inp.imp_sol == 2;
    parameters.use_legacy_solvent = inp.imp_sol == 1;
    if (inp.assume_isolated == "pcc_0d" || inp.assume_isolated == "pcc_2d")
    {
        parameters.pcc_boundary = ModulePcc::parse_boundary(inp.assume_isolated);
        parameters.pcc_2d_axis = inp.pcc_2d_axis;
    }
    if (!parameters.use_sccs && parameters.pcc_boundary == ModulePcc::Boundary::Periodic)
    {
        return parameters;
    }
    parameters.pool_process_count = pool_process_count;
    parameters.debug = inp.sccs_debug;
    parameters.expected_electron_count = electron_count;
    for (int atom_type = 0; atom_type < ucell.ntype; ++atom_type)
    {
        parameters.expected_ionic_charge
            += ucell.atoms[atom_type].ncpp.zv * ucell.atoms[atom_type].na;
    }
    if (!parameters.use_sccs)
    {
        return parameters;
    }
    const ModuleSccs::Preset preset = ModuleSccs::parse_preset(inp.sccs_preset);
    if (preset == ModuleSccs::Preset::Custom)
    {
        parameters.sccs_config.cavity.density_min = inp.sccs_rho_min;
        parameters.sccs_config.cavity.density_max = inp.sccs_rho_max;
        parameters.sccs_config.cavity.epsilon_bulk = inp.sccs_epsilon;
        parameters.sccs_config.surface_tension
            = ModuleSccs::dyn_per_cm_to_hartree_per_bohr2(inp.sccs_gamma);
        parameters.sccs_config.pressure
            = ModuleSccs::gpa_to_hartree_per_bohr3(inp.sccs_pressure);
    }
    else if (preset == ModuleSccs::Preset::Vacuum)
    {
        parameters.sccs_config = ModuleSccs::vacuum_preset();
    }
    else
    {
        parameters.sccs_config = ModuleSccs::water_preset(preset);
    }
    parameters.sccs_config.boundary = parameters.pcc_boundary;
    parameters.sccs_config.core_electrons = inp.sccs_solvent_mode == "full";
    parameters.sccs_config.core_spread = inp.sccs_corespread;
    parameters.sccs_config.cavity.lowpass_p1 = inp.sccs_lowpass_p1;
    parameters.sccs_config.cavity.lowpass_p2 = inp.sccs_lowpass_p2;
    parameters.sccs_config.max_iterations = inp.sccs_maxiter;
    parameters.sccs_config.tolerance_rms = inp.sccs_tol_rms;
    parameters.sccs_config.tolerance_max = inp.sccs_tol_max;
    parameters.sccs_config.surface_regularization = inp.sccs_surface_eta;
    parameters.sccs_config.check_fixed_point = inp.sccs_debug >= 2;
    parameters.start_drho = inp.sccs_start_drho;
    parameters.start_nmax = inp.sccs_start_nmax;
    ModuleSccs::validate_config(parameters.sccs_config);

    const double net_charge
        = parameters.expected_ionic_charge - parameters.expected_electron_count;
    // The periodic Poisson solver drops the G = 0 component of the net charge.
    // In an inhomogeneous dielectric this leaves a cell-size error with no
    // closed-form correction. Allowed, but warned. Vacuum (epsilon 1) is the
    // usual periodic charged calculation.
    if (parameters.sccs_config.boundary == ModulePcc::Boundary::Periodic
        && parameters.sccs_config.cavity.epsilon_bulk > 1.0
        && std::abs(net_charge) > parameters.normalization_tolerance)
    {
        ModuleBase::WARNING("surchem",
                            "charged SCCS with periodic boundaries: the periodic Poisson solver "
                            "drops the G = 0 component of the net charge, so the energy depends "
                            "on the cell size. Use assume_isolated pcc_0d or pcc_2d for converged "
                            "charged energies.");
    }
    return parameters;
}

void validate_kpoints(const SurchemParameters& parameters,
                           const K_Vectors& kv)
{
    const ModulePcc::Boundary boundary
        = parameters.use_sccs ? parameters.sccs_config.boundary : parameters.pcc_boundary;
    if (boundary != ModulePcc::Boundary::Pcc2d)
    {
        return;
    }

    const int axis = parameters.pcc_2d_axis;
    double maximum_open_direction_k = 0.0;
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        const ModuleBase::Vector3<double>& k_direct = kv.kvec_d[ik];
        const double open_component = axis == 0 ? k_direct.x : (axis == 1 ? k_direct.y : k_direct.z);
        const double absolute_open_direction_k = std::abs(open_component);
        maximum_open_direction_k = std::max(maximum_open_direction_k, absolute_open_direction_k);
    }
    Parallel_Reduce::reduce_max(maximum_open_direction_k);
    if (maximum_open_direction_k > 1.0e-12)
    {
        ModuleBase::WARNING_QUIT(
            "surchem",
            "pcc_2d requires Gamma-only sampling along the open lattice vector (pcc_2d_axis)");
    }
}

std::string pcc_symmetry_violation(const SurchemParameters& parameters,
                                   const ModuleSymmetry::Symmetry& symmetry)
{
    const ModulePcc::Boundary boundary
        = parameters.use_sccs ? parameters.sccs_config.boundary : parameters.pcc_boundary;
    if (boundary == ModulePcc::Boundary::Periodic)
    {
        return std::string();
    }
    for (int operation = 0; operation < symmetry.nrotk; ++operation)
    {
        const std::string violation
            = operation_violation(parameters, boundary, symmetry.gmatrix[operation], symmetry.gtrans[operation]);
        if (!violation.empty())
        {
            return violation;
        }
    }
    for (int operation = 0; operation < symmetry.nrotk_anti; ++operation)
    {
        const std::string violation = operation_violation(parameters,
                                                          boundary,
                                                          symmetry.gmatrix_anti[operation],
                                                          symmetry.gtrans_anti[operation]);
        if (!violation.empty())
        {
            return violation;
        }
    }
    // rhog_symmetry also imposes the primitive-cell translations when
    // pricell_loop is set; they must keep the PCC frame as well.
    if (ModuleSymmetry::Symmetry::pricell_loop)
    {
        const ModuleBase::Matrix3 identity;
        for (std::size_t cell = 0; cell < symmetry.ptrans.size(); ++cell)
        {
            const std::string violation = operation_violation(parameters, boundary, identity, symmetry.ptrans[cell]);
            if (!violation.empty())
            {
                return violation;
            }
        }
    }
    return std::string();
}

void validate_symmetry(const SurchemParameters& parameters, const ModuleSymmetry::Symmetry& symmetry)
{
    const std::string violation = pcc_symmetry_violation(parameters, symmetry);
    if (violation.empty())
    {
        return;
    }
    std::ostringstream message;
    message << violation << "; the PCC correction does not have this symmetry. "
            << "Use a primitive cell with the vacuum only along the open direction, or set symmetry 0 or -1";
    ModuleBase::WARNING_QUIT("surchem", message.str());
}

} // namespace ModuleSurchem
