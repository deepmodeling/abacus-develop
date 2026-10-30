#include "source_pw/module_pwdft/kernels/nonlocal_op.h"

#include <cstdint>

namespace hamilt
{

template <typename FPTYPE>
struct nonlocal_pw_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* /*dev*/,
                    const int& l1,
                    const int& l2,
                    const int& l3,
                    int& sum,
                    int& iat,
                    const int& spin,
                    const int& nkb,
                    const int& deeq_x,
                    const int& deeq_y,
                    const int& deeq_z,
                    const FPTYPE* deeq,
                    std::complex<FPTYPE>* ps,
                    const std::complex<FPTYPE>* becp)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(3)
#endif
        for (int ii = 0; ii < l1; ii++)
        {
            // each atom has nproj, means this is with structure factor;
            // each projector (each atom) must multiply coefficient
            // with all the other projectors.
            for (int jj = 0; jj < l2; ++jj)
            {
                for (int kk = 0; kk < l3; kk++)
                {
                    for (int xx = 0; xx < l3; xx++)
                    {
                        const std::int64_t projector_offset = sum + static_cast<std::int64_t>(ii) * l3;
                        const std::int64_t output_index = (projector_offset + kk) * l2 + jj;
                        const std::int64_t input_index = static_cast<std::int64_t>(jj) * nkb + projector_offset + xx;
                        const std::int64_t deeq_index = ((static_cast<std::int64_t>(spin) * deeq_x + iat + ii) * deeq_y + xx) * deeq_z + kk;
                        ps[output_index] += deeq[deeq_index] * becp[input_index];
                    }
                }
            }
        }
        sum += l1 * l3;
        iat += l1;
    }

    void operator()(const base_device::DEVICE_CPU* dev,
                    const int& l1,
                    const int& l2,
                    const int& l3,
                    int& sum,
                    int& iat,
                    const int& nkb,
                    const int& deeq_x,
                    const int& deeq_y,
                    const int& deeq_z,
                    const std::complex<FPTYPE>* deeq_nc,
                    std::complex<FPTYPE>* ps,
                    const std::complex<FPTYPE>* becp)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(3)
#endif
        for (int ii = 0; ii < l1; ii++)
        {
            // each atom has nproj, means this is with structure factor;
            // each projector (each atom) must multiply coefficient
            // with all the other projectors.
            for (int jj = 0; jj < l2; jj += 2)
            {
                for (int kk = 0; kk < l3; kk++)
                {
                    for (int xx = 0; xx < l3; xx++)
                    {
                        const std::int64_t projector_offset = sum + static_cast<std::int64_t>(ii) * l3;
                        const std::int64_t psind = (projector_offset + kk) * l2 + jj;
                        const std::int64_t becpind = static_cast<std::int64_t>(jj) * nkb + projector_offset + xx;
                        const std::int64_t spin_stride = static_cast<std::int64_t>(deeq_x) * deeq_y * deeq_z;
                        const std::int64_t deeq_index = ((static_cast<std::int64_t>(iat) + ii) * deeq_y + kk) * deeq_z + xx;
                        auto& becp1 = becp[becpind];
                        auto& becp2 = becp[becpind + nkb];
                        ps[psind] += deeq_nc[deeq_index] * becp1 + deeq_nc[spin_stride + deeq_index] * becp2;
                        ps[psind + 1] += deeq_nc[2 * spin_stride + deeq_index] * becp1 + deeq_nc[3 * spin_stride + deeq_index] * becp2;
                    } // end jj
                }
            }
        }
        iat += l1;
        sum += l1 * l3;
    }
};

template struct nonlocal_pw_op<float, base_device::DEVICE_CPU>;
template struct nonlocal_pw_op<double, base_device::DEVICE_CPU>;

} // namespace hamilt
