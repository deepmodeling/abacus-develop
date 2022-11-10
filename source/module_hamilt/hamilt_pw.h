#ifndef HAMILTPW_H
#define HAMILTPW_H

#include "hamilt.h"
#include "module_elecstate/potentials/potential_new.h"

namespace hamilt
{

class HamiltPW : public Hamilt
{
  public:
    HamiltPW(elecstate::Potential* pot_in);
    ~HamiltPW();

    // for target K point, update consequence of hPsi() and matrix()
    void updateHk(const int ik) override;

    // core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void sPsi(const std::complex<double> *psi_in, std::complex<double> *spsi, const size_t size) const override;

  private:
};

} // namespace hamilt

#endif