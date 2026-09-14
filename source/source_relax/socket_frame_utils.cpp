#include "socket_frame_utils.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace FrameUtils
{
const int kMaxJacobiSweeps = 32;

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
    for (int row = 0; row < kMatrixDimension; ++row)
    {
        const double value = values[row * kMatrixDimension + column];
        norm_squared += value * value;
    }
    return norm_squared;
}

double column_dot(const SocketFrame::Matrix9& values, int first, int second)
{
    double dot = 0.0;
    for (int row = 0; row < kMatrixDimension; ++row)
    {
        dot += values[row * kMatrixDimension + first] * values[row * kMatrixDimension + second];
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
    for (int row = 0; row < kMatrixDimension; ++row)
    {
        const int first_index = row * kMatrixDimension + first;
        const int second_index = row * kMatrixDimension + second;
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

    for (int sweep = 0; sweep < kMaxJacobiSweeps; ++sweep)
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
    for (int row = 0; row < kMatrixDimension; ++row)
    {
        for (int column = 0; column < kMatrixDimension; ++column)
        {
            long double product = 0.0L;
            for (int inner = 0; inner < kMatrixDimension; ++inner)
            {
                const int inverse_index = transpose_inverse
                                              ? column * kMatrixDimension + inner
                                              : inner * kMatrixDimension + column;
                product += static_cast<long double>(cell[row * kMatrixDimension + inner])
                           * inverse[inverse_index];
            }
            const long double expected = row == column ? 1.0L : 0.0L;
            maximum = std::max(maximum, std::fabs(product - expected));
        }
    }
    return static_cast<double>(maximum);
}
} // namespace FrameUtils
