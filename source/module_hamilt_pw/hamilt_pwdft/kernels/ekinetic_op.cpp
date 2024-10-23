#include "module_hamilt_pw/hamilt_pwdft/kernels/ekinetic_op.h"

namespace hamilt {

template <typename FPTYPE>
struct ekinetic_pw_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* /*dev*/,
                    const int& nband,
                    const int& npw,
                    const int& max_npw,
                    const bool& is_first_node,
                    const FPTYPE& tpiba2,
                    const FPTYPE* gk2_ik,
                    std::complex<FPTYPE>* tmhpsi,
                    const std::complex<FPTYPE>* tmpsi_in)
    {
        if (is_first_node)
        {
            for (int ib = 0; ib < nband; ++ib)
            {
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (int ig = 0; ig < npw; ++ig)
                {
                    tmhpsi[ig] = gk2_ik[ig] * tpiba2 * tmpsi_in[ig];
                }
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (int ig = npw; ig < max_npw; ++ig)
                {
                    tmhpsi[ig] = 0.0;
                }
                tmpsi_in += max_npw;
                tmhpsi += max_npw;
            }
        }
        else
        {
            for (int ib = 0; ib < nband; ++ib)
            {
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (int ig = 0; ig < npw; ++ig)
                {
                    tmhpsi[ig] += gk2_ik[ig] * tpiba2 * tmpsi_in[ig];
                }
                tmpsi_in += max_npw;
                tmhpsi += max_npw;
            }
        }
    }
};

template struct ekinetic_pw_op<float, base_device::DEVICE_CPU>;
template struct ekinetic_pw_op<double, base_device::DEVICE_CPU>;

}  // namespace hamilt

