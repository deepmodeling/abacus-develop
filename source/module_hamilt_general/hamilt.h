#ifndef MODULEHAMILT_H
#define MODULEHAMILT_H

#include <complex>
#include <vector>

#include "matrixblock.h"
#include "module_psi/psi.h"
#include "operator.h"

namespace hamilt
{

template<typename T, typename Device = psi::DEVICE_CPU>
class Hamilt
{
  public:
    virtual ~Hamilt(){};

    /// for target K point, update consequence of hPsi() and matrix()
    virtual void updateHk(const int ik){return;}

    /// refresh status of Hamiltonian, for example, refresh H(R) and S(R) in LCAO case
    virtual void refresh(){return;}

    /// core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const T* psi_in, T* hpsi, const size_t size) const{return;}
    virtual void sPsi(const T* psi_in, // psi
                      T* spsi,         // spsi
                      const int nrow,  // dimension of spsi: nbands * nrow
                      const int npw,   // number of plane waves
                      const int nbands // number of bands
    ) const
    {
        syncmem_op()(this->ctx, this->ctx, spsi, psi_in, static_cast<size_t>(nbands * nrow));
    }

    /// core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    virtual void matrix(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &sk_in){return;}
    virtual void matrix(MatrixBlock<double> &hk_in, MatrixBlock<double> &sk_in){return;}

    std::string classname = "none";

    int non_first_scf=0;

    /// first node operator, add operations from each operators
    Operator<T, Device>* ops = nullptr;
protected:
    Device* ctx = {};
    using syncmem_op = psi::memory::synchronize_memory_op<T, Device, Device>;
};

} // namespace hamilt

#endif