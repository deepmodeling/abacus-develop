#ifndef CAL_VEC_NORM_OP_H
#define CAL_VEC_NORM_OP_H
#include "source_base/macros.h"
namespace hamilt{
template <typename T, typename Device>
struct exx_cal_energy_op
{

    using FPTYPE = typename GetTypeReal<T>::type;
    FPTYPE operator()(const T *den, const FPTYPE *pot, FPTYPE scala, int npw);
};

} // namespace hamilt
#endif //CAL_VEC_NORM_OP_H
