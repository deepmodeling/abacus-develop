#ifndef DIAGH_H
#define DIAGH_H

#include "module_base/macros.h"
#include "module_hamilt_general/hamilt.h"
#include "module_psi/psi.h"
#include "module_hamilt_general/matrixblock.h"

#include <string>

#ifdef __MPI
#include "mpi.h"
#endif

namespace hsolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class DiagH
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    virtual ~DiagH(){};
    // virtual void init()=0;
    std::string method = "none";

    virtual void diag(hamilt::Hamilt<T, Device>* phm_in, psi::Psi<T, Device>& psi, Real* eigenvalue_in)
    {
        ModuleBase::WARNING_QUIT("diagh", "diag method not implemented for the base class!");
    };
#ifdef __MPI
    // diagnolization used in parallel-k case
    virtual void diag_pool(hamilt::MatrixBlock<T>& h_mat, hamilt::MatrixBlock<T>& s_mat, psi::Psi<T>& psi, Real* eigenvalue_in, MPI_Comm& comm)
    {
        ModuleBase::WARNING_QUIT("diagh", "diag_pool method not implemented for the base class!");
    };
#endif
};

} // namespace hsolver

#endif