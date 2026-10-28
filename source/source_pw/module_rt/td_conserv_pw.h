#ifndef TD_CONSERV_PW_H
#define TD_CONSERV_PW_H

#include "source_estate/module_charge/charge.h"

namespace pw
{
/** @brief Globally reduced electron counts and orthogonality errors. */
struct TDConserv
{
    double npsi = 0.0;
    double nrho = 0.0;
    double orth_before = 0.0;
    double orth_after = 0.0;
};

/** @brief Count each k-point pool once and integrate the charge channels. */
TDConserv cal_td_conserv(double local_npsi,
                         const Charge& charge,
                         double omega,
                         int pool_rank,
                         double orth_before,
                         double orth_after,
                         bool full_gram);
} // namespace pw
#endif
