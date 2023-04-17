#ifndef MODULEHAMILT_H
#define MODULEHAMILT_H

#include "matrixblock.h"
#include "module_psi/psi.h"
#include "operator.h"

#include <complex>
#include <vector>

namespace hamilt
{

template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class Hamilt
{
  public:
    virtual ~Hamilt(){};

    // for target K point, update consequence of hPsi() and matrix()
    virtual void updateHk(const int ik){return;}

    // core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const std::complex<FPTYPE> *psi_in, std::complex<FPTYPE> *hpsi, const size_t size) const{return;}
    virtual void sPsi(const std::complex<FPTYPE> *psi_in, std::complex<FPTYPE> *spsi, const size_t size) const{return;}


    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    virtual void matrix(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &sk_in){return;}
    virtual void matrix(MatrixBlock<double> &hk_in, MatrixBlock<double> &sk_in){return;}
    // used for tddft evolve, added by zhaoht
    virtual void matrix_l(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &hk_last_in, MatrixBlock<std::complex<double>> &sk_in){return;}
    virtual void matrix_l(MatrixBlock<double> &hk_in, MatrixBlock<double> &hk_last_in, MatrixBlock<double> &sk_in){return;}

    std::string classname = "none";

    int non_first_scf=0;

    // first node operator, add operations from each operators
    Operator<std::complex<FPTYPE>, Device>* ops = nullptr;
    Operator<double, Device>* opsd = nullptr;

};

} // namespace hamilt

#endif