#ifndef DIAGOELPA_H
#define DIAGOELPA_H

#include "diagh.h"
#include "module_basis/module_ao/parallel_orbitals.h"

namespace hsolver
{

    template<typename T>
    class DiagoElpa : public DiagH<T>
    {
    private:
        using Real = typename GetTypeReal<T>::type;

    public:
        void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in) override;
        #ifdef __MPI
        MPI_Comm setmpicomm();//set mpi comm;
        void set_comm_num(int cnum);
        #endif

    static int DecomposedState;
    

  private:
#ifdef __MPI
    bool ifElpaHandle(const bool& newIteration, const bool& ifNSCF);
    static int ifsetcomm;//need to set mpi_comm or not,-1 not,else the number of mpi needed
    static int lastmpinum;//last using mpi;
#endif
};

template <typename T>
int DiagoElpa<T>::DecomposedState = 0;
template <typename T>
int DiagoElpa<T>::lastmpinum=-1;
template <typename T>
//int DiagoElpa<T>::ifsetcomm=-1;
int DiagoElpa<T>::ifsetcomm=5;//modify only for check
} // namespace hsolver

#endif
