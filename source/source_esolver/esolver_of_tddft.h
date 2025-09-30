#ifndef ESOLVER_OF_TDDFT_H
#define ESOLVER_OF_TDDFT_H

#include "esolver_of.h"
#include "source_pw/module_ofdft/evolve_ofdft.h"

namespace ModuleESolver
{
class ESolver_OF_TDDFT : public ESolver_OF
{
  public:
    ESolver_OF_TDDFT();
    ~ESolver_OF_TDDFT();

    virtual void runner(UnitCell& ucell, const int istep) override;

  protected:
    std::vector<std::complex<double>> phi_td;                     // pphi[i] = ppsi.get_pointer(i), which will be freed in ~Psi().
    Evolve_OFDFT* evolve_ofdft=nullptr;
};
} // namespace ModuleESolver

#endif
