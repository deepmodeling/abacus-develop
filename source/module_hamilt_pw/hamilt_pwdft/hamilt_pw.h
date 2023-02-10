#ifndef HAMILTPW_H
#define HAMILTPW_H

#include "module_hamilt_general/hamilt.h"
#include "module_elecstate/potentials/potential_new.h"

namespace hamilt
{

template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class HamiltPW : public Hamilt<FPTYPE, Device>
{
  public:
    HamiltPW(elecstate::Potential* pot_in);
    template<typename T_in, typename Device_in = Device>
    explicit HamiltPW(const HamiltPW<T_in, Device_in>* hamilt);
    ~HamiltPW();

    // for target K point, update consequence of hPsi() and matrix()
    void updateHk(const int ik) override;

    // core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void sPsi(const std::complex<FPTYPE> *psi_in, std::complex<FPTYPE> *spsi, const size_t size) const override;

  private:

    Device *ctx = {};
    using syncmem_complex_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, Device>;
};

} // namespace hamilt

#endif