#ifndef EXX_STRESS_OP_H
#define EXX_STRESS_OP_H

#include "source_base/macros.h"
#include "source_base/module_device/types.h"

#include <complex>

namespace hamilt
{
template <typename T, typename Device>
struct exx_stress_accumulate_op
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const T* density_recip,
                    const Real* pot,
                    const Real* pot_stress,
                    const Real* gcar,
                    Real dkx,
                    Real dky,
                    Real dkz,
                    Real tpiba,
                    Real scalar,
                    int npw,
                    Real* sigma);
};
} // namespace hamilt

#endif // EXX_STRESS_OP_H
