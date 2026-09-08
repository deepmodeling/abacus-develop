#include "socket_frame.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
const int MATRIX_DIMENSION = 3;
const int MAX_JACOBI_SWEEPS = 32;

bool is_finite_matrix(const SocketFrame::Matrix9& values)
{
    for (std::size_t index = 0; index < values.size(); ++index)
    {
        if (!std::isfinite(values[index]))
        {
            return false;
        }
    }
    return true;
}

double column_norm_squared(const SocketFrame::Matrix9& values, int column)
{
    double norm_squared = 0.0;
    for (int row = 0; row < MATRIX_DIMENSION; ++row)
    {
        const double value = values[row * MATRIX_DIMENSION + column];
        norm_squared += value * value;
    }
    return norm_squared;
}

double column_dot(const SocketFrame::Matrix9& values, int first, int second)
{
    double dot = 0.0;
    for (int row = 0; row < MATRIX_DIMENSION; ++row)
    {
        dot += values[row * MATRIX_DIMENSION + first] * values[row * MATRIX_DIMENSION + second];
    }
    return dot;
}

bool columns_are_orthogonal(const SocketFrame::Matrix9& values)
{
    const double multiplier = 32.0 * std::numeric_limits<double>::epsilon();
    const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    for (int pair = 0; pair < 3; ++pair)
    {
        const int first = pairs[pair][0];
        const int second = pairs[pair][1];
        const double first_norm = column_norm_squared(values, first);
        const double second_norm = column_norm_squared(values, second);
        const double tolerance = multiplier * std::sqrt(first_norm * second_norm);
        if (std::fabs(column_dot(values, first, second)) > tolerance)
        {
            return false;
        }
    }
    return true;
}

void rotate_columns(SocketFrame::Matrix9& values, int first, int second, double cosine, double sine)
{
    for (int row = 0; row < MATRIX_DIMENSION; ++row)
    {
        const int first_index = row * MATRIX_DIMENSION + first;
        const int second_index = row * MATRIX_DIMENSION + second;
        const double first_value = values[first_index];
        const double second_value = values[second_index];
        values[first_index] = cosine * first_value - sine * second_value;
        values[second_index] = sine * first_value + cosine * second_value;
    }
}

bool one_sided_jacobi(SocketFrame::Matrix9& columns, SocketFrame::Matrix9& right_vectors)
{
    right_vectors = {{1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0}};
    const double multiplier = 32.0 * std::numeric_limits<double>::epsilon();
    const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};

    for (int sweep = 0; sweep < MAX_JACOBI_SWEEPS; ++sweep)
    {
        for (int pair = 0; pair < 3; ++pair)
        {
            const int first = pairs[pair][0];
            const int second = pairs[pair][1];
            const double first_norm = column_norm_squared(columns, first);
            const double second_norm = column_norm_squared(columns, second);
            const double dot = column_dot(columns, first, second);
            const double tolerance = multiplier * std::sqrt(first_norm * second_norm);
            if (std::fabs(dot) <= tolerance)
            {
                continue;
            }

            const double tau = (second_norm - first_norm) / (2.0 * dot);
            const double tangent
                = std::copysign(1.0 / (std::fabs(tau) + std::hypot(1.0, tau)), tau);
            const double cosine = 1.0 / std::sqrt(1.0 + tangent * tangent);
            const double sine = tangent * cosine;
            rotate_columns(columns, first, second, cosine, sine);
            rotate_columns(right_vectors, first, second, cosine, sine);
        }

        if (columns_are_orthogonal(columns))
        {
            return true;
        }
    }
    return false;
}

long double scaled_determinant(const SocketFrame::Matrix9& values)
{
    const long double a00 = values[0];
    const long double a01 = values[1];
    const long double a02 = values[2];
    const long double a10 = values[3];
    const long double a11 = values[4];
    const long double a12 = values[5];
    const long double a20 = values[6];
    const long double a21 = values[7];
    const long double a22 = values[8];
    return a00 * (a11 * a22 - a12 * a21)
           - a01 * (a10 * a22 - a12 * a20)
           + a02 * (a10 * a21 - a11 * a20);
}

double received_inverse_residual(const SocketFrame::Matrix9& cell,
                                 const SocketFrame::Matrix9& inverse,
                                 bool transpose_inverse)
{
    long double maximum = 0.0L;
    for (int row = 0; row < MATRIX_DIMENSION; ++row)
    {
        for (int column = 0; column < MATRIX_DIMENSION; ++column)
        {
            long double product = 0.0L;
            for (int inner = 0; inner < MATRIX_DIMENSION; ++inner)
            {
                const int inverse_index = transpose_inverse
                                              ? column * MATRIX_DIMENSION + inner
                                              : inner * MATRIX_DIMENSION + column;
                product += static_cast<long double>(cell[row * MATRIX_DIMENSION + inner])
                           * inverse[inverse_index];
            }
            const long double expected = row == column ? 1.0L : 0.0L;
            maximum = std::max(maximum, std::fabs(product - expected));
        }
    }
    return static_cast<double>(maximum);
}
} // namespace

namespace SocketFrame
{
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
    CellValidation result;
    result.ok = false;
    result.message.clear();
    result.determinant_bohr3 = 0.0;
    result.condition_number_2 = std::numeric_limits<double>::infinity();
    result.inverse_residual = std::numeric_limits<double>::infinity();
    result.computed_inverse_wire_bohr_inv.fill(0.0);

    if (!is_finite_matrix(cell_wire) || !is_finite_matrix(inverse_wire))
    {
        result.message = "cell and received inverse entries must be finite";
        return result;
    }
    if (!std::isfinite(max_condition_number) || max_condition_number <= 0.0
        || !std::isfinite(inverse_absolute_tolerance) || inverse_absolute_tolerance < 0.0
        || !std::isfinite(inverse_relative_tolerance) || inverse_relative_tolerance < 0.0)
    {
        result.message = "cell validation tolerances must be finite and nonnegative";
        return result;
    }

    double scale = 0.0;
    for (std::size_t index = 0; index < cell_wire.size(); ++index)
    {
        scale = std::max(scale, std::fabs(cell_wire[index]));
    }
    if (scale == 0.0)
    {
        result.message = "cell determinant must be positive";
        return result;
    }

    Matrix9 scaled_cell;
    for (std::size_t index = 0; index < cell_wire.size(); ++index)
    {
        scaled_cell[index] = cell_wire[index] / scale;
    }
    const long double determinant_scaled = scaled_determinant(scaled_cell);
    if (determinant_scaled <= 0.0L)
    {
        result.message = "cell determinant must be positive";
        return result;
    }
    const long double scale_long = scale;
    const long double determinant
        = determinant_scaled * scale_long * scale_long * scale_long;
    if (!std::isfinite(determinant)
        || determinant > static_cast<long double>(std::numeric_limits<double>::max()))
    {
        result.message = "cell determinant is not representable as a finite double";
        return result;
    }
    result.determinant_bohr3 = static_cast<double>(determinant);
    if (!std::isfinite(result.determinant_bohr3) || result.determinant_bohr3 <= 0.0)
    {
        result.message = "cell determinant is not representable as a positive finite double";
        return result;
    }

    Matrix9 orthogonal_columns = scaled_cell;
    Matrix9 right_vectors;
    if (!one_sided_jacobi(orthogonal_columns, right_vectors))
    {
        result.message = "cell singular-value iteration did not converge";
        return result;
    }

    double singular_values[MATRIX_DIMENSION];
    double largest_singular = 0.0;
    double smallest_singular = std::numeric_limits<double>::infinity();
    for (int column = 0; column < MATRIX_DIMENSION; ++column)
    {
        singular_values[column] = std::sqrt(column_norm_squared(orthogonal_columns, column));
        largest_singular = std::max(largest_singular, singular_values[column]);
        smallest_singular = std::min(smallest_singular, singular_values[column]);
    }
    if (smallest_singular == 0.0 || !std::isfinite(smallest_singular))
    {
        result.message = "cell is singular";
        return result;
    }
    result.condition_number_2 = largest_singular / smallest_singular;
    if (!std::isfinite(result.condition_number_2)
        || result.condition_number_2 >= max_condition_number)
    {
        result.message = "cell condition number is not below the configured maximum";
        return result;
    }

    for (int row = 0; row < MATRIX_DIMENSION; ++row)
    {
        for (int column = 0; column < MATRIX_DIMENSION; ++column)
        {
            long double inverse_value = 0.0L;
            for (int singular = 0; singular < MATRIX_DIMENSION; ++singular)
            {
                const long double sigma = singular_values[singular];
                inverse_value
                    += static_cast<long double>(right_vectors[row * MATRIX_DIMENSION + singular])
                       * orthogonal_columns[column * MATRIX_DIMENSION + singular]
                       / (static_cast<long double>(scale) * sigma * sigma);
            }
            result.computed_inverse_wire_bohr_inv[row * MATRIX_DIMENSION + column]
                = static_cast<double>(inverse_value);
        }
    }

    const double direct_inverse_residual
        = received_inverse_residual(cell_wire, inverse_wire, false);
    const double transposed_inverse_residual
        = received_inverse_residual(cell_wire, inverse_wire, true);
    result.inverse_residual = std::min(direct_inverse_residual, transposed_inverse_residual);
    const double residual_limit
        = inverse_absolute_tolerance
          + inverse_relative_tolerance * result.condition_number_2
                * std::numeric_limits<double>::epsilon();
    if (!std::isfinite(result.inverse_residual) || result.inverse_residual > residual_limit)
    {
        result.message = "received cell inverse is inconsistent with the cell";
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

VirialConversion make_ipi_virial(const Matrix9& stress_ry_per_bohr3,
                                  double volume_bohr3,
                                  double antisymmetric_absolute_tolerance,
                                  double antisymmetric_relative_tolerance)
{
    VirialConversion result;
    result.ok = false;
    result.message.clear();
    result.wire_virial_hartree.fill(0.0);
    result.max_antisymmetric_component = 0.0;

    if (!is_finite_matrix(stress_ry_per_bohr3))
    {
        result.message = "stress entries must be finite";
        return result;
    }
    if (!std::isfinite(volume_bohr3) || volume_bohr3 <= 0.0)
    {
        result.message = "cell volume must be finite and positive";
        return result;
    }
    if (!std::isfinite(antisymmetric_absolute_tolerance)
        || antisymmetric_absolute_tolerance < 0.0
        || !std::isfinite(antisymmetric_relative_tolerance)
        || antisymmetric_relative_tolerance < 0.0)
    {
        result.message = "stress symmetry tolerances must be finite and nonnegative";
        return result;
    }

    double maximum_stress = 0.0;
    for (std::size_t index = 0; index < stress_ry_per_bohr3.size(); ++index)
    {
        maximum_stress = std::max(maximum_stress, std::fabs(stress_ry_per_bohr3[index]));
    }
    for (int row = 0; row < MATRIX_DIMENSION; ++row)
    {
        for (int column = row + 1; column < MATRIX_DIMENSION; ++column)
        {
            const double difference
                = std::fabs(stress_ry_per_bohr3[row * MATRIX_DIMENSION + column]
                            - stress_ry_per_bohr3[column * MATRIX_DIMENSION + row]);
            result.max_antisymmetric_component
                = std::max(result.max_antisymmetric_component, difference);
        }
    }
    const double symmetry_limit
        = antisymmetric_absolute_tolerance + antisymmetric_relative_tolerance * maximum_stress;
    if (!std::isfinite(result.max_antisymmetric_component)
        || result.max_antisymmetric_component > symmetry_limit)
    {
        result.message = "stress tensor is not symmetric within tolerance";
        return result;
    }

    Matrix9 virial;
    for (int row = 0; row < MATRIX_DIMENSION; ++row)
    {
        for (int column = 0; column < MATRIX_DIMENSION; ++column)
        {
            const long double symmetric_stress
                = 0.5L
                  * (static_cast<long double>(stress_ry_per_bohr3[row * MATRIX_DIMENSION + column])
                     + stress_ry_per_bohr3[column * MATRIX_DIMENSION + row]);
            const long double converted = 0.5L * volume_bohr3 * symmetric_stress;
            if (!std::isfinite(converted)
                || std::fabs(converted)
                       > static_cast<long double>(std::numeric_limits<double>::max()))
            {
                result.message = "converted virial is not representable as finite doubles";
                return result;
            }
            virial[row * MATRIX_DIMENSION + column] = static_cast<double>(converted);
        }
    }
    result.wire_virial_hartree = transpose_matrix9(virial);
    result.ok = true;
    return result;
}
} // namespace SocketFrame
