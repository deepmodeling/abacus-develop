#ifndef HSOLVER_HS_OPERATOR_H
#define HSOLVER_HS_OPERATOR_H

#include "source_base/module_device/types.h"

namespace hsolver
{

/**
 * @brief What an iterative eigensolver needs from the generalized eigenproblem
 *        H x = e S x: the ability to apply H and S to a block of vectors, and
 *        nothing else.
 *
 * The solvers never see how H is built. The code that owns the Hamiltonian
 * implements this interface (hamilt::HamiltHSOperator in
 * source_hamilt/hamilt_hs_adapter.h) and the solvers only call hpsi()/spsi().
 * Unit tests implement it with a dense matrix.
 *
 * Block vectors are column major: vector i occupies x[i*ld, i*ld + ld). Only
 * the first npw rows carry data, ld is the leading dimension.
 */
template <typename T, typename Device = base_device::DEVICE_CPU>
class HSOperator
{
  public:
    virtual ~HSOperator() = default;

    /// switch H and S to k point ik; must precede hpsi()/spsi() for that k
    virtual void update_k(const int ik) = 0;

    /// hx[:, 0:nvec) = H * x[:, 0:nvec), both with leading dimension ld
    virtual void hpsi(const T* x, T* hx, const int ld, const int nvec) const = 0;

    /// sx[:, 0:nvec) = S * x[:, 0:nvec), both with leading dimension ld
    virtual void spsi(const T* x, T* sx, const int ld, const int nvec) const = 0;

    /// Hook used by DiagoIterAssist::diag_subspace_init: hcc is the n*n
    /// subspace Hamiltonian (column major, ld n). A Hamiltonian carrying a
    /// term that hpsi() does not cover (EXX in lcao_in_pw) adds it here.
    virtual void add_to_subspace_h(T* hcc, const int n) const
    {
    }

    /// Hook used by DiagoIterAssist::diag_subspace_init: vcc is the n*nband
    /// matrix of subspace eigenvectors, handed out right after the subspace
    /// diagonalization.
    virtual void export_subspace_vec(const T* vcc, const int n, const int nband) const
    {
    }
};

} // namespace hsolver

#endif // HSOLVER_HS_OPERATOR_H
