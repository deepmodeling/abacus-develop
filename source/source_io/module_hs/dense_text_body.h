#ifndef DENSE_TEXT_BODY_H
#define DENSE_TEXT_BODY_H

#include <complex>
#include <cstddef>
#include <iosfwd>
#include <vector>

namespace ModuleIO
{

/// @brief Flat index of matrix element (i, j) in a gathered buffer.
/// @details A row-major buffer stores (i, j) at i * nlocal + j, while a
///          column-major buffer stores it at j * nlocal + i. Writers and
///          readers of the flat NAO text layout must share this mapping so that
///          the file always stays row-major.
/// @param[in] i            row index, must be non-negative
/// @param[in] j            column index, must be non-negative
/// @param[in] nlocal       matrix dimension, must be positive
/// @param[in] column_major whether the buffer is stored column-major
std::size_t dense_buffer_index(const int i, const int j, const int nlocal, const bool column_major);

/// @brief Write a dense square matrix body in the flat NAO text layout shared
///        by the DM(k), H(k)/S(k) and T(k) writers.
///
/// The documented file layout is row-major: element (i, j) is emitted at
/// position (i * nlocal + j), 8 real or 4 complex values per line, and the
/// record is terminated by a single newline so that an appended record starts
/// on a fresh line.
///
/// The gathered buffer keeps the storage order of the eigensolver. For the
/// column-major solvers (genelpa, elpa, scalapack_gvx, cusolver, cusolvermp,
/// pexsi and lapack) the local block is stored as mat[i + j * nrloc], so the
/// row-major destination of Cpxgemr2d keeps a column-major buffer. Reading such
/// a buffer as row-major emits the transpose of the matrix; for a Hermitian
/// matrix that is the complex conjugate, which flips the sign of every
/// imaginary off-diagonal element. The caller selects the traversal order with
/// column_major so the file always stays row-major.
///
/// @param[in] ofs          output stream, already open by the writing rank
/// @param[in] mat          gathered nlocal-by-nlocal matrix on the writing rank
/// @param[in] nlocal       matrix dimension, must be positive
/// @param[in] precision    number of significant digits, must be positive
/// @param[in] column_major true when the gathered buffer is stored column-major
void write_dense_text_body(std::ofstream& ofs,
                           const std::vector<double>& mat,
                           const int nlocal,
                           const int precision,
                           const bool column_major);

/// @brief Complex overload of write_dense_text_body.
///
/// Real and imaginary components whose magnitude is below an internal
/// threshold are flushed to zero so that numerical noise does not clutter the
/// output.
///
/// @param[in] ofs          output stream, already open by the writing rank
/// @param[in] mat          gathered nlocal-by-nlocal matrix on the writing rank
/// @param[in] nlocal       matrix dimension, must be positive
/// @param[in] precision    number of significant digits, must be positive
/// @param[in] column_major true when the gathered buffer is stored column-major
void write_dense_text_body(std::ofstream& ofs,
                           const std::vector<std::complex<double>>& mat,
                           const int nlocal,
                           const int precision,
                           const bool column_major);

} // namespace ModuleIO

#endif // DENSE_TEXT_BODY_H
