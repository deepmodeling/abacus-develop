#ifndef EXX_Q_STATE_OP_H
#define EXX_Q_STATE_OP_H

#include "source_base/macros.h"
#include "source_base/vector3.h"
#include "source_base/module_device/types.h"
#include "source_cell/klist.h"

#include <complex>
#include <vector>

namespace ModulePW
{
class PW_Basis_K;
}

namespace hamilt
{
struct ExxSymmetryRemap
{
    std::vector<int> rep_igl;
    std::vector<int> fft_isz;
    std::vector<int> fft_ixyz;
    std::vector<std::complex<double>> phase;
    mutable int* rep_igl_device = nullptr;
    mutable int* fft_ixyz_device = nullptr;
    mutable std::complex<float>* phase_float_device = nullptr;
    mutable std::complex<double>* phase_double_device = nullptr;
};

ExxSymmetryRemap build_exx_symmetry_remap(const ModulePW::PW_Basis_K* wfcpw,
                                          const K_Vectors::ExxFullPoint& full_point,
                                          int rep_spin_index,
                                          bool need_gpu_fft_index);

template <typename T, typename Device>
struct exx_conjugate_real_op
{
    void operator()(const T* in, T* out, int nrxx);
};
} // namespace hamilt

#endif // EXX_Q_STATE_OP_H
