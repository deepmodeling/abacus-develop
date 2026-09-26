#ifndef ABACUS_XC_RVV10_PW_H
#define ABACUS_XC_RVV10_PW_H

#include "xc_rvv10.h"

#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace Rvv10
{
struct Evaluation
{
    double energy;
    double vtxc;
    std::vector<double> potential;
};

// CPU/double, full-complex FFT only. The existing PW transforms may distribute
// the density grid and G coefficients over an MPI pool; scalar outputs are
// reduced by the evaluator. PotRvv10 supplies the INPUT/SCF integration.
class Evaluator
{
  public:
    Evaluator(double b, double c);
    Evaluation evaluate(const ModulePW::PW_Basis& pw,
                        const std::vector<double>& total_density,
                        const std::vector<double>& valence_density) const;

  private:
    LocalModel model_;
    SplineBasis basis_;
    KernelTable kernel_;
};

} // namespace Rvv10
#endif
