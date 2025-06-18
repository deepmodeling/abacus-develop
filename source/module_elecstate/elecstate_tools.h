#ifndef ELECSTATE_TOOLS_H
#define ELECSTATE_TOOLS_H
#include "elecstate.h"
#include "source_base/matrix.h"

namespace elecstate
{
void calEBand(const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, fenergy& f_en);

void calculate_weights(const ModuleBase::matrix& ekb,
                       ModuleBase::matrix& wg,
                       const K_Vectors* klist,
                       efermi& eferm,
                       fenergy& f_en,
                       std::vector<double>& nelec_spin,
                       const bool skip_weights);

void fixed_weights(const std::vector<double>& ocp_kb,
                   const int& nbands,
                   const double& nelec,
                   const K_Vectors* klist,
                   ModuleBase::matrix& wg,
                   bool& skip_weights);
} // namespace elecstate

#endif