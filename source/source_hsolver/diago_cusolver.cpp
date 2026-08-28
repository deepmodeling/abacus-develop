#include "diago_cusolver.h"

#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/blacs_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"

#include <memory>
#include <type_traits>
#include <vector>

using complex = std::complex<double>;

// Namespace for the diagonalization solver
namespace hsolver
{
// Initialize the DecomposedState variable for real and complex numbers
template <typename T>
int DiagoCusolver<T>::DecomposedState = 0;

template <typename T>
DiagoCusolver<T>::DiagoCusolver(const int nlocal_in, const int nbands_in) : nlocal(nlocal_in), nbands(nbands_in)
{
}

template <typename T>
DiagoCusolver<T>::~DiagoCusolver()
{
}

// Diagonalization function
template <typename T>
void DiagoCusolver<T>::diag(
    hamilt::MatrixBlock<T>& h_mat,
    hamilt::MatrixBlock<T>& s_mat,
    psi::Psi<T>& psi,
    Real* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoCusolver", "diag");
    ModuleBase::timer::start("DiagoCusolver", "cusolver");
    // Allocate memory for eigenvalues
    std::vector<double> eigen(this->nlocal, 0.0);
    std::vector<T> eigenvectors(h_mat.row * h_mat.col);
    this->dc.Dngvd(h_mat.row, h_mat.col, h_mat.p, s_mat.p, eigen.data(), eigenvectors.data());
    const int size = psi.get_nbands() * psi.get_nbasis();
    BlasConnector::copy(size, eigenvectors.data(), 1, psi.get_pointer(), 1);
    const int inc = 1;
    BlasConnector::copy(this->nbands, eigen.data(), inc, eigenvalue_in, inc);
    ModuleBase::timer::end("DiagoCusolver", "cusolver");
}

// Explicit instantiation of the DiagoCusolver class for real and complex numbers
template class DiagoCusolver<double>;
template class DiagoCusolver<complex>;

} // namespace hsolver
