#include "dense_text_body.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>

namespace
{
/// Components with a magnitude below this threshold are flushed to zero in the
/// complex body so that numerical noise does not clutter the output.
constexpr double k_dense_text_zero_thr = 1.0e-16;

/// @brief Flat index of element (i, j) in the row-major file layout.
std::size_t row_major_index(const int i, const int j, const int nlocal)
{
    const std::size_t dim = static_cast<std::size_t>(nlocal);
    return static_cast<std::size_t>(i) * dim + static_cast<std::size_t>(j);
}

/// @brief Flat index of element (i, j) in the gathered buffer.
/// @details A column-major buffer stores (i, j) at j * nlocal + i.
std::size_t source_index(const int i, const int j, const int nlocal, const bool column_major)
{
    return ModuleIO::dense_buffer_index(i, j, nlocal, column_major);
}

/// @brief Validate the arguments shared by both overloads.
void check_arguments(const std::size_t mat_size, const int nlocal, const int precision)
{
    assert(nlocal > 0);
    assert(precision > 0);
    const std::size_t expected = static_cast<std::size_t>(nlocal) * static_cast<std::size_t>(nlocal);
    assert(mat_size >= expected);
    (void)mat_size;
    (void)expected;
}

/// @brief Apply the common stream formatting of the flat text body.
void set_body_format(std::ofstream& ofs, const int precision)
{
    ofs << std::fixed;
    ofs << std::scientific;
    ofs << std::setprecision(precision);
    ofs << std::right;
}
} // namespace

std::size_t ModuleIO::dense_buffer_index(const int i, const int j, const int nlocal, const bool column_major)
{
    assert(i >= 0);
    assert(j >= 0);
    assert(nlocal > 0);
    if (column_major)
    {
        return row_major_index(j, i, nlocal);
    }
    return row_major_index(i, j, nlocal);
}

void ModuleIO::write_dense_text_body(std::ofstream& ofs,
                                     const std::vector<double>& mat,
                                     const int nlocal,
                                     const int precision,
                                     const bool column_major)
{
    check_arguments(mat.size(), nlocal, precision);
    set_body_format(ofs, precision);
    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nlocal; ++j)
        {
            if (j % 8 == 0)
            {
                ofs << "\n";
            }
            ofs << " " << mat[source_index(i, j, nlocal, column_major)];
        }
    }
    // Terminate the record so that an appended record starts on a fresh line.
    ofs << std::endl;
}

void ModuleIO::write_dense_text_body(std::ofstream& ofs,
                                     const std::vector<std::complex<double>>& mat,
                                     const int nlocal,
                                     const int precision,
                                     const bool column_major)
{
    check_arguments(mat.size(), nlocal, precision);
    set_body_format(ofs, precision);
    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nlocal; ++j)
        {
            if (j % 4 == 0)
            {
                ofs << "\n";
            }
            const std::complex<double>& value = mat[source_index(i, j, nlocal, column_major)];
            double real_v = std::real(value);
            if (std::abs(real_v) < k_dense_text_zero_thr)
            {
                real_v = 0.0;
            }
            double imag_v = std::imag(value);
            if (std::abs(imag_v) < k_dense_text_zero_thr)
            {
                imag_v = 0.0;
            }
            ofs << " (" << real_v << "," << imag_v << ")";
        }
    }
    // Terminate the record so that an appended record starts on a fresh line.
    ofs << std::endl;
}
