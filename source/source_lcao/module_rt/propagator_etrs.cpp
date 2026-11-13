#include "propagator.h"

#include "source_base/matrix2d_linalg.h"

#include <complex>
#include <vector>

namespace module_rt
{
#ifdef __MPI
void Propagator::compute_propagator_etrs(
    Matrix2dView<const std::complex<double>> Stmp,
    Matrix2dView<const std::complex<double>> Htmp,
    Matrix2dView<const std::complex<double>> H_laststep,
    Matrix2dView<std::complex<double>> U_operator,
    std::ofstream& ofs_running,
    const int print_matrix) const
{
    std::vector<std::complex<double>> first_storage(Htmp.size());
    std::vector<std::complex<double>> second_storage(Htmp.size());
    Matrix2dView<std::complex<double>> first(
        first_storage.data(), Htmp.nrloc(), Htmp.ncloc(), Htmp.desc());
    Matrix2dView<std::complex<double>> second(
        second_storage.data(), Htmp.nrloc(), Htmp.ncloc(), Htmp.desc());
    const int tag = 2;
    compute_propagator_taylor(
        Stmp, Htmp, first, ofs_running, print_matrix, tag);
    compute_propagator_taylor(
        Stmp, H_laststep, second, ofs_running, print_matrix, tag);
    matrix2d_linalg::matmul(U_operator, first, second);
}
#endif // __MPI
} // namespace module_rt
