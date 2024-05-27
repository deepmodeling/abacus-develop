#ifndef DIAGH_H
#define DIAGH_H

#include <string>

#include "module_base/macros.h"

#include "module_hamilt_general/hamilt.h"
#include "module_psi/psi.h"

#ifdef __MPI
#include "mpi.h"
#endif

template<typename T> struct consts
{
    consts();
    T zero;
    T one;
    T neg_one;
};

namespace hsolver
{


struct diag_comm_info
{

    const int rank;
    const int nproc;

#ifndef __MPI
    diag_comm_info(const int rank_in, const int nproc_in) : rank(rank_in), nproc(nproc_in) {}
#else
    const MPI_Comm comm;
    diag_comm_info(const MPI_Comm &comm_in, const int rank_in, const int nproc_in) : comm(comm_in), rank(rank_in), nproc(nproc_in) {}
#endif
};


template <typename T, typename Device = base_device::DEVICE_CPU>
class DiagH
{
  private:
    using Real = typename GetTypeReal<T>::type;
  public:
    virtual ~DiagH(){};
    // virtual void init()=0;
    std::string method = "none";

    virtual void diag(hamilt::Hamilt<T, Device> *phm_in, psi::Psi<T, Device> &psi, Real *eigenvalue_in) {
      ModuleBase::WARNING_QUIT("diagh", "diag method not implemented for the base class!");
    };

};

} // namespace hsolver

#endif