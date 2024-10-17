#ifndef DIAGOELPA_H
#define DIAGOELPA_H

#include "diagh.h"
#include "module_basis/module_ao/parallel_orbitals.h"

namespace hsolver
{

template <typename T>
class DiagoElpa
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in);
#ifdef __MPI
    // diagnolization used in parallel-k case
    void diag_pool(hamilt::MatrixBlock<T>& h_mat, hamilt::MatrixBlock<T>& s_mat, psi::Psi<T>& psi, Real* eigenvalue_in, MPI_Comm& comm);
    MPI_Comm setmpicomm(); // set mpi comm;
    static int elpa_num_thread;  // need to set mpi_comm or not,-1 not,else the number of mpi needed
#endif

    static int DecomposedState;

  private:
#ifdef __MPI
    bool ifElpaHandle(const bool& newIteration, const bool& ifNSCF);
    static int lastmpinum; // last using mpi;
#endif
};

template <typename T>
int DiagoElpa<T>::DecomposedState = 0;
template <typename T>
int DiagoElpa<T>::lastmpinum = -1;
template <typename T>
int DiagoElpa<T>::elpa_num_thread = -1;
} // namespace hsolver

#endif
