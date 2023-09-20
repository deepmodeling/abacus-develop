#ifndef MODULEHAMILT_H
#define MODULEHAMILT_H

#include "matrixblock.h"
#include "module_psi/psi.h"
#include "operator.h"

#include <complex>
#include <vector>

namespace hamilt
{

template<typename T, typename Device = psi::DEVICE_CPU>
class Hamilt
{
  public:
    virtual ~Hamilt(){};

    // for target K point, update consequence of hPsi() and matrix()
    virtual void updateHk(const int ik){return;}

    // core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const T* psi_in, T* hpsi, const size_t size) const{return;}
    virtual void sPsi(const T* psi_in, T* spsi, const size_t size) const{return;}

    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    virtual void matrix(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &sk_in){return;}
    virtual void matrix(MatrixBlock<double> &hk_in, MatrixBlock<double> &sk_in){return;}

    std::string classname = "none";

    int non_first_scf=0;

    // first node operator, add operations from each operators
    Operator<T, Device>* ops = nullptr;
    Operator<double, Device>* opsd = nullptr;

};

} // namespace hamilt

#endif