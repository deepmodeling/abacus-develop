#include "socket_frame.h"

#include "source_relax/socket_frame_utils.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace FrameUtils;

namespace SocketFrame
{
namespace
{
CellValidation make_failed_cell_validation()
{
    CellValidation result;
    result.ok = false;
    result.message.clear();
    result.determinant_bohr3 = 0.0;
    result.condition_number_2 = std::numeric_limits<double>::infinity();
    result.inverse_residual = std::numeric_limits<double>::infinity();
    result.computed_inverse_wire_bohr_inv.fill(0.0);
    return result;
}

bool validate_cell_entries(const Matrix9& cell_wire,
                           const Matrix9& inverse_wire,
                           std::string& message)
{
    if (!is_finite_matrix(cell_wire) || !is_finite_matrix(inverse_wire))
    {
        message = "cell and received inverse entries must be finite";
        return false;
    }
    return true;
}

bool validate_cell_tolerances(const double max_condition_number,
                              const double inverse_absolute_tolerance,
                              const double inverse_relative_tolerance,
                              std::string& message)
{
    if (!std::isfinite(max_condition_number) || max_condition_number <= 0.0
        || !std::isfinite(inverse_absolute_tolerance) || inverse_absolute_tolerance < 0.0
        || !std::isfinite(inverse_relative_tolerance) || inverse_relative_tolerance < 0.0)
    {
        message = "cell validation tolerances must be finite and nonnegative";
        return false;
    }
    return true;
}

bool compute_cell_scale(const Matrix9& cell_wire, double& scale, Matrix9& scaled_cell,
                        std::string& message)
{
    scale = 0.0;
    for (std::size_t index = 0; index < cell_wire.size(); ++index)
    {
        scale = std::max(scale, std::fabs(cell_wire[index]));
    }
    if (scale == 0.0)
    {
        message = "cell determinant must be positive";
        return false;
    }
    for (std::size_t index = 0; index < cell_wire.size(); ++index)
    {
        scaled_cell[index] = cell_wire[index] / scale;
    }
    return true;
}

bool compute_cell_determinant(const double scale,
                              const Matrix9& scaled_cell,
                              double& determinant_bohr3,
                              std::string& message)
{
    const long double determinant_scaled = scaled_determinant(scaled_cell);
    if (determinant_scaled <= 0.0L)
    {
        message = "cell determinant must be positive";
        return false;
    }
    const long double scale_long = scale;
    const long double determinant
        = determinant_scaled * scale_long * scale_long * scale_long;
    if (!std::isfinite(determinant)
        || determinant > static_cast<long double>(std::numeric_limits<double>::max()))
    {
        message = "cell determinant is not representable as a finite double";
        return false;
    }
    determinant_bohr3 = static_cast<double>(determinant);
    if (!std::isfinite(determinant_bohr3) || determinant_bohr3 <= 0.0)
    {
        message = "cell determinant is not representable as a positive finite double";
        return false;
    }
    return true;
}

bool compute_cell_svd(const Matrix9& scaled_cell,
                      double singular_values[kMatrixDimension],
                      Matrix9& orthogonal_columns,
                      Matrix9& right_vectors,
                      std::string& message)
{
    orthogonal_columns = scaled_cell;
    if (!one_sided_jacobi(orthogonal_columns, right_vectors))
    {
        message = "cell singular-value iteration did not converge";
        return false;
    }
    for (int column = 0; column < kMatrixDimension; ++column)
    {
        singular_values[column] = std::sqrt(column_norm_squared(orthogonal_columns, column));
    }
    return true;
}

bool compute_condition_number(const double singular_values[kMatrixDimension],
                              const double max_condition_number,
                              double& condition_number_2,
                              std::string& message)
{
    double largest_singular = 0.0;
    double smallest_singular = std::numeric_limits<double>::infinity();
    for (int column = 0; column < kMatrixDimension; ++column)
    {
        largest_singular = std::max(largest_singular, singular_values[column]);
        smallest_singular = std::min(smallest_singular, singular_values[column]);
    }
    if (smallest_singular == 0.0 || !std::isfinite(smallest_singular))
    {
        message = "cell is singular";
        return false;
    }
    condition_number_2 = largest_singular / smallest_singular;
    if (!std::isfinite(condition_number_2)
        || condition_number_2 >= max_condition_number)
    {
        message = "cell condition number is not below the configured maximum";
        return false;
    }
    return true;
}

void compute_cell_inverse(const double scale,
                          const double singular_values[kMatrixDimension],
                          const Matrix9& orthogonal_columns,
                          const Matrix9& right_vectors,
                          Matrix9& computed_inverse)
{
    for (int row = 0; row < kMatrixDimension; ++row)
    {
        for (int column = 0; column < kMatrixDimension; ++column)
        {
            long double inverse_value = 0.0L;
            for (int singular = 0; singular < kMatrixDimension; ++singular)
            {
                const long double sigma = singular_values[singular];
                inverse_value
                    += static_cast<long double>(right_vectors[row * kMatrixDimension + singular])
                       * orthogonal_columns[column * kMatrixDimension + singular]
                       / (static_cast<long double>(scale) * sigma * sigma);
            }
            computed_inverse[row * kMatrixDimension + column]
                = static_cast<double>(inverse_value);
        }
    }
}

bool check_inverse_consistency(const Matrix9& cell_wire,
                               const Matrix9& inverse_wire,
                               const double condition_number_2,
                               const double inverse_absolute_tolerance,
                               const double inverse_relative_tolerance,
                               double& inverse_residual,
                               std::string& message)
{
    const double direct_inverse_residual
        = received_inverse_residual(cell_wire, inverse_wire, false);
    const double transposed_inverse_residual
        = received_inverse_residual(cell_wire, inverse_wire, true);
    inverse_residual = std::min(direct_inverse_residual, transposed_inverse_residual);
    const double residual_limit
        = inverse_absolute_tolerance
          + inverse_relative_tolerance * condition_number_2
                * std::numeric_limits<double>::epsilon();
    if (!std::isfinite(inverse_residual) || inverse_residual > residual_limit)
    {
        message = "received cell inverse is inconsistent with the cell";
        return false;
    }
    return true;
}
} // namespace

Matrix9 transpose_matrix9(const Matrix9& values)
{
    return {{values[0], values[3], values[6],
             values[1], values[4], values[7],
             values[2], values[5], values[8]}};
}

CellValidation validate_ipi_cell(const Matrix9& cell_wire,
                                 const Matrix9& inverse_wire,
                                 double max_condition_number,
                                 double inverse_absolute_tolerance,
                                 double inverse_relative_tolerance)
{
    CellValidation result = make_failed_cell_validation();

    if (!validate_cell_entries(cell_wire, inverse_wire, result.message))
    {
        return result;
    }
    if (!validate_cell_tolerances(max_condition_number,
                                  inverse_absolute_tolerance,
                                  inverse_relative_tolerance,
                                  result.message))
    {
        return result;
    }

    double scale = 0.0;
    Matrix9 scaled_cell;
    if (!compute_cell_scale(cell_wire, scale, scaled_cell, result.message))
    {
        return result;
    }
    if (!compute_cell_determinant(scale, scaled_cell, result.determinant_bohr3, result.message))
    {
        return result;
    }

    double singular_values[kMatrixDimension];
    Matrix9 orthogonal_columns;
    Matrix9 right_vectors;
    if (!compute_cell_svd(scaled_cell, singular_values, orthogonal_columns, right_vectors,
                          result.message))
    {
        return result;
    }
    if (!compute_condition_number(singular_values, max_condition_number,
                                  result.condition_number_2, result.message))
    {
        return result;
    }

    compute_cell_inverse(scale, singular_values, orthogonal_columns, right_vectors,
                         result.computed_inverse_wire_bohr_inv);
    if (!check_inverse_consistency(cell_wire, inverse_wire, result.condition_number_2,
                                   inverse_absolute_tolerance, inverse_relative_tolerance,
                                   result.inverse_residual, result.message))
    {
        return result;
    }

    result.ok = true;
    return result;
}

bool validate_positions(const std::vector<double>& positions_bohr,
                        std::size_t coordinate_count,
                        std::string& message)
{
    if (positions_bohr.size() != coordinate_count)
    {
        message = "position coordinate count does not match the validated atom count";
        return false;
    }
    for (std::size_t index = 0; index < positions_bohr.size(); ++index)
    {
        if (!std::isfinite(positions_bohr[index]))
        {
            message = "position coordinates must be finite";
            return false;
        }
    }
    message.clear();
    return true;
}

bool checked_position_count(std::int32_t nat_socket,
                            int nat_expected,
                            std::size_t& coordinate_count,
                            std::string& message)
{
    if (nat_socket != nat_expected)
    {
        message = "socket atom count does not match the expected atom count";
        return false;
    }
    if (nat_socket < 0)
    {
        message = "socket atom count must not be negative";
        return false;
    }
    const std::size_t atom_count = static_cast<std::size_t>(nat_socket);
    if (atom_count > std::numeric_limits<std::size_t>::max() / 3)
    {
        message = "socket position coordinate count is not representable";
        return false;
    }
    coordinate_count = 3 * atom_count;
    message.clear();
    return true;
}

namespace
{
VirialConversion make_failed_virial_conversion()
{
    VirialConversion result;
    result.ok = false;
    result.message.clear();
    result.wire_virial_hartree.fill(0.0);
    result.max_antisymmetric_component = 0.0;
    return result;
}

bool validate_virial_inputs(const Matrix9& stress_ry_per_bohr3,
                            const double volume_bohr3,
                            const double antisymmetric_absolute_tolerance,
                            const double antisymmetric_relative_tolerance,
                            std::string& message)
{
    if (!is_finite_matrix(stress_ry_per_bohr3))
    {
        message = "stress entries must be finite";
        return false;
    }
    if (!std::isfinite(volume_bohr3) || volume_bohr3 <= 0.0)
    {
        message = "cell volume must be finite and positive";
        return false;
    }
    if (!std::isfinite(antisymmetric_absolute_tolerance)
        || antisymmetric_absolute_tolerance < 0.0
        || !std::isfinite(antisymmetric_relative_tolerance)
        || antisymmetric_relative_tolerance < 0.0)
    {
        message = "stress symmetry tolerances must be finite and nonnegative";
        return false;
    }
    return true;
}

bool check_stress_symmetry(const Matrix9& stress_ry_per_bohr3,
                           const double antisymmetric_absolute_tolerance,
                           const double antisymmetric_relative_tolerance,
                           double& max_antisymmetric_component,
                           std::string& message)
{
    double maximum_stress = 0.0;
    for (std::size_t index = 0; index < stress_ry_per_bohr3.size(); ++index)
    {
        maximum_stress = std::max(maximum_stress, std::fabs(stress_ry_per_bohr3[index]));
    }
    for (int row = 0; row < kMatrixDimension; ++row)
    {
        for (int column = row + 1; column < kMatrixDimension; ++column)
        {
            const double difference
                = std::fabs(stress_ry_per_bohr3[row * kMatrixDimension + column]
                            - stress_ry_per_bohr3[column * kMatrixDimension + row]);
            max_antisymmetric_component
                = std::max(max_antisymmetric_component, difference);
        }
    }
    const double symmetry_limit
        = antisymmetric_absolute_tolerance + antisymmetric_relative_tolerance * maximum_stress;
    if (!std::isfinite(max_antisymmetric_component)
        || max_antisymmetric_component > symmetry_limit)
    {
        message = "stress tensor is not symmetric within tolerance";
        return false;
    }
    return true;
}

bool compute_symmetric_virial(const Matrix9& stress_ry_per_bohr3,
                              const double volume_bohr3,
                              Matrix9& wire_virial_hartree,
                              std::string& message)
{
    Matrix9 virial;
    for (int row = 0; row < kMatrixDimension; ++row)
    {
        for (int column = 0; column < kMatrixDimension; ++column)
        {
            const long double symmetric_stress
                = 0.5L
                  * (static_cast<long double>(stress_ry_per_bohr3[row * kMatrixDimension + column])
                     + stress_ry_per_bohr3[column * kMatrixDimension + row]);
            const long double converted = 0.5L * volume_bohr3 * symmetric_stress;
            if (!std::isfinite(converted)
                || std::fabs(converted)
                       > static_cast<long double>(std::numeric_limits<double>::max()))
            {
                message = "converted virial is not representable as finite doubles";
                return false;
            }
            virial[row * kMatrixDimension + column] = static_cast<double>(converted);
        }
    }
    wire_virial_hartree = transpose_matrix9(virial);
    return true;
}
} // namespace

VirialConversion make_ipi_virial(const Matrix9& stress_ry_per_bohr3,
                                  double volume_bohr3,
                                  double antisymmetric_absolute_tolerance,
                                  double antisymmetric_relative_tolerance)
{
    VirialConversion result = make_failed_virial_conversion();

    if (!validate_virial_inputs(stress_ry_per_bohr3,
                                volume_bohr3,
                                antisymmetric_absolute_tolerance,
                                antisymmetric_relative_tolerance,
                                result.message))
    {
        return result;
    }
    if (!check_stress_symmetry(stress_ry_per_bohr3,
                               antisymmetric_absolute_tolerance,
                               antisymmetric_relative_tolerance,
                               result.max_antisymmetric_component,
                               result.message))
    {
        return result;
    }
    if (!compute_symmetric_virial(stress_ry_per_bohr3,
                                  volume_bohr3,
                                  result.wire_virial_hartree,
                                  result.message))
    {
        return result;
    }

    result.ok = true;
    return result;
}
} // namespace SocketFrame
