#ifndef HAMILTPW_H
#define HAMILTPW_H

#include "hamilt.h"

namespace hamilt
{

class HamiltPW : public Hamilt
{
  public:
    HamiltPW();
    ~HamiltPW();

    // for target K point, update consequence of hPsi() and matrix()
    void updateHk(const int ik) override;

    // core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const override;
    virtual void sPsi(const std::complex<double> *psi_in, std::complex<double> *spsi, const size_t size) const override;

  private:
};

} // namespace hamilt

#endif