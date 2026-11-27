#ifndef ABACUS_SOURCE_HAMILT_MODULE_GINT_GINT_ENV_GAMMA_H
#define ABACUS_SOURCE_HAMILT_MODULE_GINT_GINT_ENV_GAMMA_H

#include "gint.h"
#include "gint_info.h"
#include <memory>
#include <vector>

namespace ModuleGint
{

class Gint_env_gamma : public Gint
{
  public:
    Gint_env_gamma(
        const double *psid,
        const Parallel_Orbitals *pv,
        const int nbands,
        const int nlocal,
        double *rho);

    void cal_env_band(const int iband);

  private:
    // output
    double *rho_ = nullptr;

    // intermediate variable
    std::vector<double> wfc_gint_;
};

} // namespace ModuleGint

#endif // ABACUS_SOURCE_HAMILT_MODULE_GINT_GINT_ENV_GAMMA_H
