#include "source_base/macros.h"

#ifndef CAL_DENSITY_REAL_OP_H
#define CAL_DENSITY_REAL_OP_H
namespace hamilt
{
template <typename T, typename Device>
struct cal_density_real_op
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const T *psi1, const T* psi2, T *out, double omega, int nrxx);
};
}
#endif //CAL_DENSITY_REAL_OP_H
