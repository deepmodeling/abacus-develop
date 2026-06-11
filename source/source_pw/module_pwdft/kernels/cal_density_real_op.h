#include "source_base/macros.h"

#ifndef CAL_DENSITY_REAL_OP_H
#define CAL_DENSITY_REAL_OP_H
namespace hamilt
{
template <typename T, typename Device>
struct cal_density_real_op
{
    using Real = typename GetTypeReal<T>::type;

    // Single element-wise operation
    void operator()(const T *psi1, const T* psi2, T *out, double omega, int nrxx);

    // Batched element-wise operation (GPU only - CPU falls back to loop)
    void operator_batch(
        const T *psi_nk,           // Constant input (nrxx)
        const T *psi_mq_batch,     // Batch input (batch_size × nrxx)
        T *density_batch,          // Batch output (batch_size × nrxx)
        double omega,
        int nrxx,
        int batch_size);
};
}
#endif //CAL_DENSITY_REAL_OP_H
