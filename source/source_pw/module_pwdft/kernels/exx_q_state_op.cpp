#include "source_pw/module_pwdft/kernels/exx_q_state_op.h"

#include "source_base/constants.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_basis_k.h"

#include <complex>
#include <cmath>
#include <unordered_map>

namespace hamilt
{
namespace
{
struct IntGKey
{
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const IntGKey& rhs) const
    {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }
};

struct IntGKeyHash
{
    std::size_t operator()(const IntGKey& key) const
    {
        std::size_t h = static_cast<std::size_t>(key.x * 73856093);
        h ^= static_cast<std::size_t>(key.y * 19349663);
        h ^= static_cast<std::size_t>(key.z * 83492791);
        return h;
    }
};

IntGKey make_int_g_key(const ModuleBase::Vector3<double>& g)
{
    IntGKey key;
    key.x = static_cast<int>(std::lround(g.x));
    key.y = static_cast<int>(std::lround(g.y));
    key.z = static_cast<int>(std::lround(g.z));
    return key;
}

ModuleBase::Vector3<double> generic_g_direct_from_ig(const ModulePW::PW_Basis_K* wfcpw, int ig)
{
    const int isz = wfcpw->ig2isz[ig];
    const int iz_raw = isz % wfcpw->nz;
    const int is = isz / wfcpw->nz;
    const int ixy = wfcpw->is2fftixy[is];
    int ix = ixy / wfcpw->fftny;
    int iy = ixy % wfcpw->fftny;
    int iz = iz_raw;
    if (ix >= int(wfcpw->nx / 2) + 1)
    {
        ix -= wfcpw->nx;
    }
    if (iy >= int(wfcpw->ny / 2) + 1)
    {
        iy -= wfcpw->ny;
    }
    if (iz >= int(wfcpw->nz / 2) + 1)
    {
        iz -= wfcpw->nz;
    }
    return ModuleBase::Vector3<double>(ix, iy, iz);
}

ModuleBase::Vector3<double> generic_g_cartesian_from_ig(const ModulePW::PW_Basis_K* wfcpw, int ig)
{
    return generic_g_direct_from_ig(wfcpw, ig) * wfcpw->G;
}
} // namespace

ExxSymmetryRemap build_exx_symmetry_remap(const ModulePW::PW_Basis_K* wfcpw,
                                          const K_Vectors::ExxFullPoint& full_point,
                                          int rep_spin_index,
                                          bool need_gpu_fft_index)
{
    const ModuleBase::Vector3<double> raw_rep = full_point.full_kvec_d * full_point.kgmatrix;
    const ModuleBase::Vector3<double> rep_shift = raw_rep - wfcpw->kvec_d[rep_spin_index];

    std::unordered_map<IntGKey, int, IntGKeyHash> rep_g_to_ig;
    rep_g_to_ig.reserve(wfcpw->npwk[rep_spin_index]);
    for (int ig_rep = 0; ig_rep < wfcpw->npwk[rep_spin_index]; ++ig_rep)
    {
        rep_g_to_ig.emplace(make_int_g_key(wfcpw->getgdirect(rep_spin_index, ig_rep)), ig_rep);
    }

    ExxSymmetryRemap remap;
    remap.rep_igl.reserve(wfcpw->npw);
    remap.fft_isz.reserve(wfcpw->npw);
    remap.phase.reserve(wfcpw->npw);
    if (need_gpu_fft_index)
    {
        remap.fft_ixyz.reserve(wfcpw->npw);
    }

    for (int ig = 0; ig < wfcpw->npw; ++ig)
    {
        const ModuleBase::Vector3<double> g_cart = generic_g_cartesian_from_ig(wfcpw, ig);
        const ModuleBase::Vector3<double> gplus_full = g_cart + full_point.full_kvec_c;
        if (gplus_full.norm2() > wfcpw->gk_ecut + 1e-10)
        {
            continue;
        }
        const ModuleBase::Vector3<double> g_full = generic_g_direct_from_ig(wfcpw, ig);
        const ModuleBase::Vector3<double> g_rep = g_full * full_point.kgmatrix + rep_shift;
        const auto it = rep_g_to_ig.find(make_int_g_key(g_rep));
        if (it == rep_g_to_ig.end())
        {
            ModuleBase::WARNING_QUIT("build_exx_symmetry_remap",
                                     "failed to map full-point G vector to representative G vector");
        }

        const int ig_rep = it->second;
        remap.rep_igl.push_back(ig_rep);
        remap.fft_isz.push_back(wfcpw->ig2isz[ig]);
        const ModuleBase::Vector3<double> gk_rep = wfcpw->getgdirect(rep_spin_index, ig_rep)
                                                   + wfcpw->kvec_d[rep_spin_index];
        const double phase_arg = ModuleBase::TWO_PI * (gk_rep * full_point.gtrans);
        remap.phase.push_back(std::complex<double>(std::cos(phase_arg), std::sin(phase_arg)));

        if (need_gpu_fft_index)
        {
            const int isz = wfcpw->ig2isz[ig];
            const int iz = isz % wfcpw->nz;
            const int is = isz / wfcpw->nz;
            const int ixy = wfcpw->is2fftixy[is];
            const int iy = ixy % wfcpw->ny;
            const int ix = ixy / wfcpw->ny;
            remap.fft_ixyz.push_back(iz + iy * wfcpw->nz + ix * wfcpw->ny * wfcpw->nz);
        }
    }

    if (remap.rep_igl.empty())
    {
        ModuleBase::WARNING_QUIT("build_exx_symmetry_remap", "empty full-point G-vector map");
    }
    if (need_gpu_fft_index && remap.fft_ixyz.size() != remap.rep_igl.size())
    {
        ModuleBase::WARNING_QUIT("build_exx_symmetry_remap", "incomplete full-point GPU FFT map");
    }
    return remap;
}

template <typename FPTYPE>
struct exx_conjugate_real_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>
{
    using T = std::complex<FPTYPE>;
    void operator()(const T* in, T* out, int nrxx)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < nrxx; ++ir)
        {
            out[ir] = std::conj(in[ir]);
        }
    }
};

template struct exx_conjugate_real_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct exx_conjugate_real_op<std::complex<double>, base_device::DEVICE_CPU>;
} // namespace hamilt
