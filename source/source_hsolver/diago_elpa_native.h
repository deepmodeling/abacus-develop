#ifndef DIAGOELPANATIVE_H
#define DIAGOELPANATIVE_H

#include "source_base/macros.h"   // GetRealType
#include "source_hamilt/hamilt.h"
#include "source_basis/module_ao/parallel_orbitals.h"

namespace hsolver
{

template <typename T>
class DiagoElpaNative
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    /// @param nlocal_in global dimension of the NAO Hamiltonian
    /// @param nbands_in number of lowest eigenpairs to compute
    /// @param use_gpu_in offload to the NVIDIA-GPU ELPA kernels when ELPA was built with GPU support
    DiagoElpaNative(const int nlocal_in, const int nbands_in, const bool use_gpu_in)
        : nlocal(nlocal_in), nbands(nbands_in), use_gpu(use_gpu_in) {};

    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in);
#ifdef __MPI
    // diagnolization used in parallel-k case
    void diag_pool(hamilt::MatrixBlock<T>& h_mat, hamilt::MatrixBlock<T>& s_mat, psi::Psi<T>& psi, Real* eigenvalue_in, MPI_Comm& comm);
    MPI_Comm setmpicomm(); // set mpi comm;
    static int elpa_num_thread;  // need to set mpi_comm or not,-1 not,else the number of mpi needed
    static int lastmpinum; // last using mpi;

#endif

    static int DecomposedState;

  private:
    const int nlocal;
    const int nbands;
    const bool use_gpu;
};

template <typename T>
int DiagoElpaNative<T>::DecomposedState = 0;
#ifdef __MPI

template <typename T>
int DiagoElpaNative<T>::lastmpinum = -1;
template <typename T>
int DiagoElpaNative<T>::elpa_num_thread = -1;
#endif

} // namespace hsolver

#endif
