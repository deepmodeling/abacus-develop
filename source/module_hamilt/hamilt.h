#ifndef MODULEHAMILT_H
#define MODULEHAMILT_H

#include "matrixblock.h"
#include "module_psi/psi.h"
#include "operator.h"

#include <complex>
#include <vector>

namespace hamilt
{

class Hamilt
{
  public:
    virtual ~Hamilt(){};
    // construct Hamiltonian matrix with inputed electonic density
    virtual void constructHamilt() {return;}

    // for target K point, update consequence of hPsi() and matrix()
    virtual void updateHk(const int ik){return;}

    // core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const{return;}
    virtual void sPsi(const std::complex<double> *psi_in, std::complex<double> *spsi, const size_t size) const{return;}

    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    virtual void matrix(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &sk_in){return;}
    virtual void matrix(MatrixBlock<double> &hk_in, MatrixBlock<double> &sk_in){return;}

    std::string classname = "none";

    int non_first_scf=0;

  protected:
    // array, save operations from each operators
    std::vector<Operator*> ops;
};

} // namespace hamilt

#endif