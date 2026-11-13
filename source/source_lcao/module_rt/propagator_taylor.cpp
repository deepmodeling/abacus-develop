#include "propagator.h"

#include "source_base/matrix2d_linalg.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <vector>

namespace module_rt
{
#ifdef __MPI
void Propagator::compute_propagator_taylor(
    Matrix2dView<const std::complex<double>> Stmp,
    Matrix2dView<const std::complex<double>> Htmp,
    Matrix2dView<std::complex<double>> U_operator,
    std::ofstream& ofs_running,
    const int print_matrix,
    const int tag) const
{
    assert(Htmp.size() > 0);
    const int nrloc = Htmp.nrloc();
    const int ncloc = Htmp.ncloc();
    const std::size_t matrix_size = Htmp.size();
    const int* descriptor = Htmp.desc();

    if (print_matrix)
    {
        ofs_running << std::endl;
        ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < nrloc; ++i)
        {
            for (int j = 0; j < ncloc; ++j)
            {
                const std::complex<double>& value = Stmp.at_local(i, j);
                ofs_running << value.real() << "+" << value.imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < nrloc; ++i)
        {
            for (int j = 0; j < ncloc; ++j)
            {
                const std::complex<double>& value = Htmp.at_local(i, j);
                ofs_running << value.real() << "+" << value.imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    std::vector<std::complex<double>> rank0(matrix_size);
    std::vector<std::complex<double>> rank2(matrix_size);
    std::vector<std::complex<double>> rank3(matrix_size);
    std::vector<std::complex<double>> rank4(matrix_size);
    std::vector<std::complex<double>> sinv(matrix_size);

    Matrix2dView<std::complex<double>> rank0_view(
        rank0.data(), nrloc, ncloc, descriptor);
    Matrix2dView<std::complex<double>> rank2_view(
        rank2.data(), nrloc, ncloc, descriptor);
    Matrix2dView<std::complex<double>> rank3_view(
        rank3.data(), nrloc, ncloc, descriptor);
    Matrix2dView<std::complex<double>> rank4_view(
        rank4.data(), nrloc, ncloc, descriptor);
    Matrix2dView<std::complex<double>> sinv_view(
        sinv.data(), nrloc, ncloc, descriptor);

    std::copy(Stmp.data(), Stmp.data() + matrix_size, sinv.data());
    const std::complex<double> beta(0.0, -0.5 * this->dt / tag);
    std::copy(Htmp.data(), Htmp.data() + matrix_size, U_operator.data());
    matrix2d_linalg::scale(U_operator, beta);
    matrix2d_linalg::solve(sinv_view, U_operator);
    matrix2d_linalg::set_identity(rank0_view);
    matrix2d_linalg::matmul(rank2_view, U_operator, U_operator);
    matrix2d_linalg::matmul(rank3_view, U_operator, rank2_view);
    matrix2d_linalg::matmul(rank4_view, U_operator, rank3_view);

    const std::complex<double> p1(1.0, 0.0);
    const std::complex<double> p2(1.0 / 2.0, 0.0);
    const std::complex<double> p3(1.0 / 6.0, 0.0);
    const std::complex<double> p4(1.0 / 24.0, 0.0);
    matrix2d_linalg::linear_combination(
        U_operator, p1, rank0_view, p1, U_operator);
    matrix2d_linalg::linear_combination(
        U_operator, p1, U_operator, p2, rank2_view);
    matrix2d_linalg::linear_combination(
        U_operator, p1, U_operator, p3, rank3_view);
    matrix2d_linalg::linear_combination(
        U_operator, p1, U_operator, p4, rank4_view);

    if (print_matrix)
    {
        ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < nrloc; ++i)
        {
            for (int j = 0; j < ncloc; ++j)
            {
                const std::complex<double>& value = U_operator.at_local(i, j);
                const double real_part = std::abs(value.real()) < 1.0e-8
                    ? 0.0
                    : value.real();
                const double imaginary_part = std::abs(value.imag()) < 1.0e-8
                    ? 0.0
                    : value.imag();
                ofs_running << real_part << "+" << imaginary_part << "i ";
            }
            ofs_running << std::endl;
        }
    }
}
#endif // __MPI
} // namespace module_rt
