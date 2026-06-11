#ifndef MUL_POTENTIAL_OP_H
#define MUL_POTENTIAL_OP_H
#include "source_base/macros.h"

namespace hamilt {

template <typename T, typename Device>
struct mul_potential_op
{
//     int npw = rhopw_dev->npw;
//     int nks = wfcpw->nks;
//     int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;
//     int nk = nks / nk_fac;
//
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static)
// #endif
//     for (int ig = 0; ig < npw; ig++)
//     {
//         int ig_kq = ik * nks * npw + iq * npw + ig;
//         density_recip[ig] *= pot[ig_kq];
//     }
    using FPTYPE = typename GetTypeReal<T>::type;

    // Single element-wise operation
    void operator()(const FPTYPE *pot, T *density_recip, int npw, int nks, int ik, int iq);

    // Batched element-wise operation (GPU only - CPU caller should loop)
    // Applies the same potential to multiple density arrays
    void operator_batch(
        const FPTYPE *pot,           // Potential (npw, same for all batch elements)
        T *density_recip_batch,      // Batch of density arrays (batch_size × npw)
        int npw,                     // Number of plane waves
        int batch_size);             // Number of batch elements
};

} // namespace hamilt

#endif // MUL_POTENTIAL_OP_H
