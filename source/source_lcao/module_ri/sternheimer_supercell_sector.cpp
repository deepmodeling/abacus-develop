#include "source_lcao/module_ri/sternheimer_supercell_sector.h"

#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ModuleRI
{
namespace
{

using Complex = std::complex<double>;

std::size_t column_major_index(const int row, const int column, const int rows)
{
    return static_cast<std::size_t>(row)
           + static_cast<std::size_t>(rows) * static_cast<std::size_t>(column);
}

int validated_cell_count(const std::array<int, 3>& repeats,
                         const SternheimerReducedKPoint& primitive_kpoint)
{
    int cell_count = 1;
    for (int direction = 0; direction != 3; ++direction)
    {
        const int repeat = repeats[static_cast<std::size_t>(direction)];
        const double kpoint = primitive_kpoint[static_cast<std::size_t>(direction)];
        if (repeat <= 0 || !std::isfinite(kpoint)
            || std::abs(kpoint * repeat - std::round(kpoint * repeat)) > 1.0e-10)
        {
            throw std::invalid_argument(
                "Sternheimer supercell sector requires positive repeats and a commensurate finite k point.");
        }
        if (cell_count > std::numeric_limits<int>::max() / repeat)
        {
            throw std::overflow_error("Sternheimer supercell sector cell count overflowed.");
        }
        cell_count *= repeat;
    }
    return cell_count;
}

std::vector<Complex> translation_sector_basis(const int primitive_basis_size,
                                              const std::array<int, 3>& repeats,
                                              const SternheimerReducedKPoint& primitive_kpoint)
{
    const int cell_count = repeats[0] * repeats[1] * repeats[2];
    const int full_basis_size = cell_count * primitive_basis_size;
    std::vector<Complex> basis(
        static_cast<std::size_t>(full_basis_size) * static_cast<std::size_t>(primitive_basis_size),
        Complex(0.0, 0.0));
    const double normalization = 1.0 / std::sqrt(static_cast<double>(cell_count));
    const double two_pi = 2.0 * std::acos(-1.0);
    for (int ix = 0; ix != repeats[0]; ++ix)
    {
        for (int iy = 0; iy != repeats[1]; ++iy)
        {
            for (int iz = 0; iz != repeats[2]; ++iz)
            {
                const int cell_index = (ix * repeats[1] + iy) * repeats[2] + iz;
                const double angle = two_pi
                                     * (primitive_kpoint[0] * ix + primitive_kpoint[1] * iy
                                        + primitive_kpoint[2] * iz);
                const Complex phase = normalization * std::exp(Complex(0.0, angle));
                for (int orbital = 0; orbital != primitive_basis_size; ++orbital)
                {
                    const int row = cell_index * primitive_basis_size + orbital;
                    basis[column_major_index(row, orbital, full_basis_size)] = phase;
                }
            }
        }
    }
    return basis;
}

SternheimerSupercellSector recover_sector_from_eigenbasis_coefficients(
    const std::vector<double>& eigenvalues,
    const std::vector<Complex>& eigenbasis_sector_coefficients,
    const std::array<int, 3>& repeats,
    const SternheimerReducedKPoint& primitive_kpoint,
    const int primitive_basis_size,
    const int full_basis_size)
{
    const std::size_t expected_coefficients
        = static_cast<std::size_t>(full_basis_size)
          * static_cast<std::size_t>(primitive_basis_size);
    if (eigenbasis_sector_coefficients.size() != expected_coefficients)
    {
        throw std::invalid_argument("Recovered Sternheimer sector coefficient dimensions are inconsistent.");
    }

    const std::size_t sector_matrix_size
        = static_cast<std::size_t>(primitive_basis_size)
          * static_cast<std::size_t>(primitive_basis_size);
    std::vector<Complex> overlap(sector_matrix_size, Complex(0.0, 0.0));
    std::vector<Complex> hamiltonian(sector_matrix_size, Complex(0.0, 0.0));
    for (int column = 0; column != primitive_basis_size; ++column)
    {
        for (int row = 0; row != primitive_basis_size; ++row)
        {
            Complex overlap_element(0.0, 0.0);
            Complex hamiltonian_element(0.0, 0.0);
            for (int band = 0; band != full_basis_size; ++band)
            {
                const Complex bra = std::conj(eigenbasis_sector_coefficients[
                    column_major_index(band, row, full_basis_size)]);
                const Complex ket = eigenbasis_sector_coefficients[
                    column_major_index(band, column, full_basis_size)];
                overlap_element += bra * ket;
                hamiltonian_element
                    += bra * eigenvalues[static_cast<std::size_t>(band)] * ket;
            }
            const std::size_t index
                = column_major_index(row, column, primitive_basis_size);
            overlap[index] = overlap_element;
            hamiltonian[index] = hamiltonian_element;
        }
    }
    const std::vector<Complex> overlap_original = overlap;

    constexpr int problem_type = 1;
    constexpr char jobz = 'V';
    constexpr char uplo = 'U';
    std::vector<double> sector_eigenvalues(
        static_cast<std::size_t>(primitive_basis_size), 0.0);
    int info = 0;
    int lwork = -1;
    Complex work_query(0.0, 0.0);
    std::vector<double> rwork(
        static_cast<std::size_t>(std::max(1, 3 * primitive_basis_size - 2)), 0.0);
    zhegv_(&problem_type,
           &jobz,
           &uplo,
           &primitive_basis_size,
           hamiltonian.data(),
           &primitive_basis_size,
           overlap.data(),
           &primitive_basis_size,
           sector_eigenvalues.data(),
           &work_query,
           &lwork,
           rwork.data(),
           &info);
    if (info != 0)
    {
        throw std::runtime_error("Sternheimer supercell sector LAPACK workspace query failed.");
    }
    lwork = std::max(1, static_cast<int>(std::real(work_query)));
    std::vector<Complex> work(static_cast<std::size_t>(lwork), Complex(0.0, 0.0));
    zhegv_(&problem_type,
           &jobz,
           &uplo,
           &primitive_basis_size,
           hamiltonian.data(),
           &primitive_basis_size,
           overlap.data(),
           &primitive_basis_size,
           sector_eigenvalues.data(),
           work.data(),
           &lwork,
           rwork.data(),
           &info);
    if (info != 0)
    {
        std::ostringstream message;
        message << "Sternheimer supercell sector generalized eigensolve failed; zhegv info="
                << info << '.';
        throw std::runtime_error(message.str());
    }

    const std::vector<Complex> sector_basis
        = translation_sector_basis(primitive_basis_size, repeats, primitive_kpoint);
    SternheimerSupercellSector result;
    result.eigenvalues = std::move(sector_eigenvalues);
    result.coefficients.assign(
        static_cast<std::size_t>(primitive_basis_size),
        std::vector<Complex>(static_cast<std::size_t>(full_basis_size),
                             Complex(0.0, 0.0)));
    for (int state = 0; state != primitive_basis_size; ++state)
    {
        for (int primitive_orbital = 0; primitive_orbital != primitive_basis_size;
             ++primitive_orbital)
        {
            const Complex coefficient
                = hamiltonian[column_major_index(primitive_orbital,
                                                 state,
                                                 primitive_basis_size)];
            for (int full_orbital = 0; full_orbital != full_basis_size; ++full_orbital)
            {
                result.coefficients[static_cast<std::size_t>(state)]
                                   [static_cast<std::size_t>(full_orbital)]
                    += sector_basis[column_major_index(full_orbital,
                                                       primitive_orbital,
                                                       full_basis_size)]
                       * coefficient;
            }
        }
    }

    for (int column = 0; column != primitive_basis_size; ++column)
    {
        for (int row = 0; row != primitive_basis_size; ++row)
        {
            Complex metric_element(0.0, 0.0);
            for (int left = 0; left != primitive_basis_size; ++left)
            {
                for (int right = 0; right != primitive_basis_size; ++right)
                {
                    metric_element
                        += std::conj(hamiltonian[column_major_index(left,
                                                                  row,
                                                                  primitive_basis_size)])
                           * overlap_original[column_major_index(left,
                                                                 right,
                                                                 primitive_basis_size)]
                           * hamiltonian[column_major_index(right,
                                                            column,
                                                            primitive_basis_size)];
                }
            }
            const Complex expected
                = row == column ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
            result.max_orthonormality_error
                = std::max(result.max_orthonormality_error,
                           std::abs(metric_element - expected));
        }
    }

    for (int state = 0; state != primitive_basis_size; ++state)
    {
        double residual_squared = 0.0;
        double norm_squared = 0.0;
        for (int full_band = 0; full_band != full_basis_size; ++full_band)
        {
            Complex eigenbasis_coefficient(0.0, 0.0);
            for (int primitive_orbital = 0; primitive_orbital != primitive_basis_size;
                 ++primitive_orbital)
            {
                eigenbasis_coefficient
                    += eigenbasis_sector_coefficients[column_major_index(
                           full_band, primitive_orbital, full_basis_size)]
                       * hamiltonian[column_major_index(primitive_orbital,
                                                        state,
                                                        primitive_basis_size)];
            }
            const double coefficient_norm = std::norm(eigenbasis_coefficient);
            const double difference = eigenvalues[static_cast<std::size_t>(full_band)]
                                      - result.eigenvalues[static_cast<std::size_t>(state)];
            norm_squared += coefficient_norm;
            residual_squared += difference * difference * coefficient_norm;
        }
        const double relative_residual
            = norm_squared > 0.0 ? std::sqrt(residual_squared / norm_squared)
                                 : std::sqrt(residual_squared);
        result.max_full_space_residual
            = std::max(result.max_full_space_residual, relative_residual);
    }
    return result;
}

} // namespace

SternheimerSupercellSector recover_sternheimer_supercell_sector(
    const std::vector<double>& eigenvalues,
    const std::vector<std::vector<std::complex<double>>>& eigenvectors,
    const std::array<int, 3>& repeats,
    const SternheimerReducedKPoint& primitive_kpoint)
{
    const int cell_count = validated_cell_count(repeats, primitive_kpoint);
    const int full_basis_size = static_cast<int>(eigenvectors.size());
    if (full_basis_size <= 0 || eigenvalues.size() != eigenvectors.size()
        || full_basis_size % cell_count != 0)
    {
        throw std::invalid_argument(
            "Sternheimer supercell sector requires a complete square eigensystem and translation-major AO blocks.");
    }
    for (const auto& eigenvector: eigenvectors)
    {
        if (eigenvector.size() != static_cast<std::size_t>(full_basis_size))
        {
            throw std::invalid_argument("Sternheimer supercell sector eigensystem is not square.");
        }
    }

    const int primitive_basis_size = full_basis_size / cell_count;
    std::vector<Complex> eigenvector_matrix(
        static_cast<std::size_t>(full_basis_size) * static_cast<std::size_t>(full_basis_size));
    for (int band = 0; band != full_basis_size; ++band)
    {
        for (int orbital = 0; orbital != full_basis_size; ++orbital)
        {
            eigenvector_matrix[column_major_index(orbital, band, full_basis_size)]
                = eigenvectors[static_cast<std::size_t>(band)][static_cast<std::size_t>(orbital)];
        }
    }

    std::vector<Complex> eigenbasis_sector_coefficients
        = translation_sector_basis(primitive_basis_size, repeats, primitive_kpoint);
    std::vector<int> pivots(static_cast<std::size_t>(full_basis_size), 0);
    int info = 0;
    zgesv_(&full_basis_size,
           &primitive_basis_size,
           eigenvector_matrix.data(),
           &full_basis_size,
           pivots.data(),
           eigenbasis_sector_coefficients.data(),
           &full_basis_size,
           &info);
    if (info != 0)
    {
        std::ostringstream message;
        message << "Sternheimer supercell sector could not invert the complete LCAO eigenvector matrix; zgesv info="
                << info << '.';
        throw std::runtime_error(message.str());
    }

    const std::size_t sector_matrix_size
        = static_cast<std::size_t>(primitive_basis_size) * static_cast<std::size_t>(primitive_basis_size);
    std::vector<Complex> overlap(sector_matrix_size, Complex(0.0, 0.0));
    std::vector<Complex> hamiltonian(sector_matrix_size, Complex(0.0, 0.0));
    for (int column = 0; column != primitive_basis_size; ++column)
    {
        for (int row = 0; row != primitive_basis_size; ++row)
        {
            Complex overlap_element(0.0, 0.0);
            Complex hamiltonian_element(0.0, 0.0);
            for (int band = 0; band != full_basis_size; ++band)
            {
                const Complex bra = std::conj(eigenbasis_sector_coefficients[
                    column_major_index(band, row, full_basis_size)]);
                const Complex ket = eigenbasis_sector_coefficients[
                    column_major_index(band, column, full_basis_size)];
                overlap_element += bra * ket;
                hamiltonian_element += bra * eigenvalues[static_cast<std::size_t>(band)] * ket;
            }
            const std::size_t index = column_major_index(row, column, primitive_basis_size);
            overlap[index] = overlap_element;
            hamiltonian[index] = hamiltonian_element;
        }
    }
    const std::vector<Complex> overlap_original = overlap;

    constexpr int problem_type = 1;
    constexpr char jobz = 'V';
    constexpr char uplo = 'U';
    std::vector<double> sector_eigenvalues(static_cast<std::size_t>(primitive_basis_size), 0.0);
    int lwork = -1;
    Complex work_query(0.0, 0.0);
    std::vector<double> rwork(static_cast<std::size_t>(std::max(1, 3 * primitive_basis_size - 2)), 0.0);
    zhegv_(&problem_type,
           &jobz,
           &uplo,
           &primitive_basis_size,
           hamiltonian.data(),
           &primitive_basis_size,
           overlap.data(),
           &primitive_basis_size,
           sector_eigenvalues.data(),
           &work_query,
           &lwork,
           rwork.data(),
           &info);
    if (info != 0)
    {
        throw std::runtime_error("Sternheimer supercell sector LAPACK workspace query failed.");
    }
    lwork = std::max(1, static_cast<int>(std::real(work_query)));
    std::vector<Complex> work(static_cast<std::size_t>(lwork), Complex(0.0, 0.0));
    zhegv_(&problem_type,
           &jobz,
           &uplo,
           &primitive_basis_size,
           hamiltonian.data(),
           &primitive_basis_size,
           overlap.data(),
           &primitive_basis_size,
           sector_eigenvalues.data(),
           work.data(),
           &lwork,
           rwork.data(),
           &info);
    if (info != 0)
    {
        std::ostringstream message;
        message << "Sternheimer supercell sector generalized eigensolve failed; zhegv info=" << info << '.';
        throw std::runtime_error(message.str());
    }

    const std::vector<Complex> sector_basis
        = translation_sector_basis(primitive_basis_size, repeats, primitive_kpoint);
    SternheimerSupercellSector result;
    result.eigenvalues = std::move(sector_eigenvalues);
    result.coefficients.assign(static_cast<std::size_t>(primitive_basis_size),
                               std::vector<Complex>(static_cast<std::size_t>(full_basis_size),
                                                    Complex(0.0, 0.0)));
    for (int state = 0; state != primitive_basis_size; ++state)
    {
        for (int primitive_orbital = 0; primitive_orbital != primitive_basis_size; ++primitive_orbital)
        {
            const Complex coefficient
                = hamiltonian[column_major_index(primitive_orbital, state, primitive_basis_size)];
            for (int full_orbital = 0; full_orbital != full_basis_size; ++full_orbital)
            {
                result.coefficients[static_cast<std::size_t>(state)][static_cast<std::size_t>(full_orbital)]
                    += sector_basis[column_major_index(full_orbital,
                                                       primitive_orbital,
                                                       full_basis_size)]
                       * coefficient;
            }
        }
    }

    for (int column = 0; column != primitive_basis_size; ++column)
    {
        for (int row = 0; row != primitive_basis_size; ++row)
        {
            Complex metric_element(0.0, 0.0);
            for (int left = 0; left != primitive_basis_size; ++left)
            {
                for (int right = 0; right != primitive_basis_size; ++right)
                {
                    metric_element
                        += std::conj(hamiltonian[column_major_index(left, row, primitive_basis_size)])
                           * overlap_original[column_major_index(left, right, primitive_basis_size)]
                           * hamiltonian[column_major_index(right, column, primitive_basis_size)];
                }
            }
            const Complex expected = row == column ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
            result.max_orthonormality_error
                = std::max(result.max_orthonormality_error, std::abs(metric_element - expected));
        }
    }

    for (int state = 0; state != primitive_basis_size; ++state)
    {
        double residual_squared = 0.0;
        double norm_squared = 0.0;
        for (int full_band = 0; full_band != full_basis_size; ++full_band)
        {
            Complex eigenbasis_coefficient(0.0, 0.0);
            for (int primitive_orbital = 0; primitive_orbital != primitive_basis_size;
                 ++primitive_orbital)
            {
                eigenbasis_coefficient
                    += eigenbasis_sector_coefficients[column_major_index(full_band,
                                                                         primitive_orbital,
                                                                         full_basis_size)]
                       * hamiltonian[column_major_index(primitive_orbital,
                                                        state,
                                                        primitive_basis_size)];
            }
            const double coefficient_norm = std::norm(eigenbasis_coefficient);
            const double difference = eigenvalues[static_cast<std::size_t>(full_band)]
                                      - result.eigenvalues[static_cast<std::size_t>(state)];
            norm_squared += coefficient_norm;
            residual_squared += difference * difference * coefficient_norm;
        }
        const double relative_residual
            = norm_squared > 0.0 ? std::sqrt(residual_squared / norm_squared) : std::sqrt(residual_squared);
        result.max_full_space_residual = std::max(result.max_full_space_residual, relative_residual);
    }
    return result;
}

std::vector<SternheimerSupercellKPointSector> recover_all_sternheimer_supercell_sectors(
    const std::vector<double>& eigenvalues,
    const std::vector<std::vector<std::complex<double>>>& eigenvectors,
    const std::array<int, 3>& repeats)
{
    const int cell_count = validated_cell_count(repeats, {0.0, 0.0, 0.0});
    const int full_basis_size = static_cast<int>(eigenvectors.size());
    if (full_basis_size <= 0 || eigenvalues.size() != eigenvectors.size()
        || full_basis_size % cell_count != 0)
    {
        throw std::invalid_argument(
            "Sternheimer supercell sectors require a complete square eigensystem and translation-major AO blocks.");
    }
    for (const auto& eigenvector: eigenvectors)
    {
        if (eigenvector.size() != static_cast<std::size_t>(full_basis_size))
        {
            throw std::invalid_argument("Sternheimer supercell sector eigensystem is not square.");
        }
    }
    const int primitive_basis_size = full_basis_size / cell_count;

    std::vector<SternheimerSupercellKPointSector> sectors;
    sectors.reserve(static_cast<std::size_t>(cell_count));
    std::vector<Complex> all_sector_right_hand_sides(
        static_cast<std::size_t>(full_basis_size)
            * static_cast<std::size_t>(full_basis_size),
        Complex(0.0, 0.0));
    for (int iz = 0; iz != repeats[2]; ++iz)
    {
        for (int iy = 0; iy != repeats[1]; ++iy)
        {
            for (int ix = 0; ix != repeats[0]; ++ix)
            {
                SternheimerSupercellKPointSector record;
                record.kpoint = {static_cast<double>(ix) / repeats[0],
                                 static_cast<double>(iy) / repeats[1],
                                 static_cast<double>(iz) / repeats[2]};
                const std::vector<Complex> sector_basis
                    = translation_sector_basis(primitive_basis_size,
                                               repeats,
                                               record.kpoint);
                const std::size_t column_offset
                    = sectors.size() * static_cast<std::size_t>(primitive_basis_size)
                      * static_cast<std::size_t>(full_basis_size);
                std::copy(sector_basis.begin(),
                          sector_basis.end(),
                          all_sector_right_hand_sides.begin()
                              + static_cast<std::ptrdiff_t>(column_offset));
                sectors.push_back(std::move(record));
            }
        }
    }

    std::vector<Complex> eigenvector_matrix(
        static_cast<std::size_t>(full_basis_size)
        * static_cast<std::size_t>(full_basis_size));
    for (int band = 0; band != full_basis_size; ++band)
    {
        for (int orbital = 0; orbital != full_basis_size; ++orbital)
        {
            eigenvector_matrix[column_major_index(orbital, band, full_basis_size)]
                = eigenvectors[static_cast<std::size_t>(band)]
                              [static_cast<std::size_t>(orbital)];
        }
    }
    std::vector<int> pivots(static_cast<std::size_t>(full_basis_size), 0);
    int info = 0;
    zgesv_(&full_basis_size,
           &full_basis_size,
           eigenvector_matrix.data(),
           &full_basis_size,
           pivots.data(),
           all_sector_right_hand_sides.data(),
           &full_basis_size,
           &info);
    if (info != 0)
    {
        std::ostringstream message;
        message << "Sternheimer supercell sectors could not invert the complete LCAO eigenvector matrix; zgesv info="
                << info << '.';
        throw std::runtime_error(message.str());
    }

    const std::size_t sector_coefficient_count
        = static_cast<std::size_t>(full_basis_size)
          * static_cast<std::size_t>(primitive_basis_size);
    for (std::size_t ik = 0; ik != sectors.size(); ++ik)
    {
        const auto begin = all_sector_right_hand_sides.begin()
                           + static_cast<std::ptrdiff_t>(ik * sector_coefficient_count);
        const std::vector<Complex> sector_coefficients(
            begin, begin + static_cast<std::ptrdiff_t>(sector_coefficient_count));
        sectors[ik].sector = recover_sector_from_eigenbasis_coefficients(
            eigenvalues,
            sector_coefficients,
            repeats,
            sectors[ik].kpoint,
            primitive_basis_size,
            full_basis_size);
    }
    return sectors;
}

} // namespace ModuleRI
