#include "diago_cusolver.h"

#include "source_base/tool_title.h"
#include "source_base/timer.h"

#include <cstddef>
#include <cstdio>
#include <cstdlib>

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

namespace
{
void solve(
    Diag_Cusolver_gvd &solver,
    const int n,
    const int nbands,
    const double *h,
    const double *s,
    double *w,
    double *v)
{
    if (Diag_Cusolver_gvd::use_subset(n, nbands))
    {
        if (solver.solve_subset(n, nbands, h, s, w, v))
        {
            return;
        }
        std::fprintf(
            stderr, "cuSOLVER Dsygvdx could not return %d eigenpairs (n=%d); retrying with Dsygvd\n", nbands, n);
    }
    solver.solve_full(n, nbands, h, s, w, v);
}

void solve(
    Diag_Cusolver_gvd &solver,
    const int n,
    const int nbands,
    const complex *h,
    const complex *s,
    double *w,
    complex *v)
{
    solver.solve_full(n, nbands, h, s, w, v);
}
} // namespace

template <typename T>
void DiagoCusolver<T>::diag(
    ModuleBase::MatrixBlock<T>& h_mat,
    ModuleBase::MatrixBlock<T>& s_mat,
    psi::Psi<T>& psi,
    Real* eigenvalue_in)
{
    ModuleBase::timer::start("DiagoCusolver", "cusolver");
    ModuleBase::TITLE("DiagoCusolver", "diag");

    const std::size_t dimension = static_cast<std::size_t>(this->nlocal);
    if (this->nlocal <= 0 || this->nbands <= 0 || this->nbands > this->nlocal || h_mat.row != dimension ||
        h_mat.col != dimension || s_mat.row != dimension || s_mat.col != dimension ||
        psi.get_nbasis() != this->nlocal || psi.get_nbands() < this->nbands || !psi.get_k_first() ||
        psi.get_current_b() != 0 || !h_mat.p || !s_mat.p || !eigenvalue_in || !psi.get_pointer())
    {
        std::fprintf(
            stderr, "cuSOLVER invalid H/S or wavefunction layout (n=%d, nbands=%d)\n", this->nlocal, this->nbands);
        std::exit(EXIT_FAILURE);
    }

    T *eigenvectors = psi.get_pointer();
    solve(this->dc, this->nlocal, this->nbands, h_mat.p, s_mat.p, eigenvalue_in, eigenvectors);

    ModuleBase::timer::end("DiagoCusolver", "cusolver");
}

// Explicit instantiation of the DiagoCusolver class for real and complex numbers
template class DiagoCusolver<double>;
template class DiagoCusolver<complex>;

} // namespace hsolver
