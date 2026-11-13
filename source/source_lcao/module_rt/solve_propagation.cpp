#include "solve_propagation.h"

#include "source_base/constants.h"
#include "source_base/matrix2d_linalg.h"

#include <algorithm>
#include <complex>
#include <vector>

namespace module_rt
{
#ifdef __MPI
namespace
{
void solve_propagation_impl(
    const double dt,
    Matrix2dView<const std::complex<double>> overlap,
    Matrix2dView<const std::complex<double>> hamiltonian,
    const Matrix2dView<const std::complex<double>>* gauge,
    Matrix2dView<const std::complex<double>> psi_last,
    Matrix2dView<std::complex<double>> psi)
{
    std::vector<std::complex<double>> operator_A(hamiltonian.size());
    std::vector<std::complex<double>> operator_B(hamiltonian.size());
    Matrix2dView<const std::complex<double>> effective_hamiltonian = hamiltonian;
    std::vector<std::complex<double>> effective_hamiltonian_storage;
    if (gauge != nullptr)
    {
        effective_hamiltonian_storage.resize(hamiltonian.size());
        effective_hamiltonian = Matrix2dView<const std::complex<double>>(
            effective_hamiltonian_storage.data(),
            hamiltonian.nrloc(),
            hamiltonian.ncloc(),
            hamiltonian.desc());
        Matrix2dView<std::complex<double>> effective_output(
            effective_hamiltonian_storage.data(),
            hamiltonian.nrloc(),
            hamiltonian.ncloc(),
            hamiltonian.desc());
        const std::complex<double> one(1.0, 0.0);
        matrix2d_linalg::linear_combination(
            effective_output, one, hamiltonian, one, *gauge);
    }

    const double dt_au = dt / ModuleBase::AU_to_FS;
    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double> beta1(0.0, -0.25 * dt_au);
    const std::complex<double> beta2(0.0, 0.25 * dt_au);

    Matrix2dView<std::complex<double>> operator_A_view(
        operator_A.data(),
        hamiltonian.nrloc(),
        hamiltonian.ncloc(),
        hamiltonian.desc());
    Matrix2dView<std::complex<double>> operator_B_view(
        operator_B.data(),
        hamiltonian.nrloc(),
        hamiltonian.ncloc(),
        hamiltonian.desc());
    matrix2d_linalg::linear_combination(
        operator_A_view, alpha, overlap, beta2, effective_hamiltonian);
    matrix2d_linalg::linear_combination(
        operator_B_view, alpha, overlap, beta1, effective_hamiltonian);

    std::vector<std::complex<double>> rhs(psi_last.size());
    Matrix2dView<std::complex<double>> rhs_view(
        rhs.data(), psi_last.nrloc(), psi_last.ncloc(), psi_last.desc());
    matrix2d_linalg::matmul(rhs_view, operator_B_view, psi_last);
    matrix2d_linalg::solve(operator_A_view, rhs_view);
    std::copy(rhs.begin(), rhs.end(), psi.data());
}
} // namespace

void solve_propagation(
    const double dt,
    Matrix2dView<const std::complex<double>> overlap,
    Matrix2dView<const std::complex<double>> hamiltonian,
    Matrix2dView<const std::complex<double>> psi_last,
    Matrix2dView<std::complex<double>> psi)
{
    solve_propagation_impl(
        dt, overlap, hamiltonian, nullptr, psi_last, psi);
}

void solve_propagation(
    const double dt,
    Matrix2dView<const std::complex<double>> overlap,
    Matrix2dView<const std::complex<double>> hamiltonian,
    Matrix2dView<const std::complex<double>> gauge,
    Matrix2dView<const std::complex<double>> psi_last,
    Matrix2dView<std::complex<double>> psi)
{
    solve_propagation_impl(
        dt, overlap, hamiltonian, &gauge, psi_last, psi);
}
#endif // __MPI
} // namespace module_rt
