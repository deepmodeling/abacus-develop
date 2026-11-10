#include "fft_cuda.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/device_check.h"
#if defined(__CUFFTDX)
#include "source_basis/module_pw/kernels/fused_veff.h"
#endif

#include <algorithm>
#include <cassert>

namespace ModuleBase
{
namespace
{
// Grids smaller than this are dominated by the extra kernel launches of the
// pruned transform, so the regular cuFFT path is cheaper.
constexpr int pruned_fft_min_grid_points = 1 << 20;
// Keep one quarter of the x planes inactive before using the pruned plan.
constexpr int pruned_fft_active_numerator = 3;
constexpr int pruned_fft_active_denominator = 4;
// Upper bound for the fused spinor workspace, in bytes.
constexpr size_t spinor_workspace_bytes = 256ULL * 1024 * 1024;
constexpr int max_spinor_batch = 16;
} // namespace

template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::initfft(int nx_in, int ny_in, int nz_in)
{
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
}

template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::initfft(int nx_in,
                               int ny_in,
                               int nz_in,
                               int lixy_in,
                               int rixy_in,
                               int ns_in,
                               int /*nplane_in*/,
                               int nproc_in,
                               bool gamma_only_in,
                               bool xprime_in)
{
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
    this->lixy_index = lixy_in;
    this->rixy_index = rixy_in;
    this->nproc = nproc_in;
    this->ns_index = ns_in;
    this->gamma_only = gamma_only_in;
    this->xprime = xprime_in;
}

template <>
void FFT_CUDA<float>::setupFFT()
{
}

template <>
void FFT_CUDA<double>::setupFFT()
{
}

template <>
void FFT_CUDA<float>::cleanFFT()
{
    if (c_handle)
    {
        cufftDestroy(c_handle);
        c_handle = {};
    }
}

template <>
void FFT_CUDA<double>::cleanFFT()
{
    if (z_handle)
    {
        cufftDestroy(z_handle);
        z_handle = {};
    }
    for (auto* handle : {&z_x_handle, &z_yz_left_handle, &z_yz_right_handle,
                         &z_stick_handle[0], &z_stick_handle[1], &z_stick_handle[2]})
    {
        if (*handle)
        {
            cufftDestroy(*handle);
            *handle = {};
        }
    }
}

template <>
void FFT_CUDA<float>::clear()
{
    this->cleanFFT();
    if (c_auxr_3d != nullptr)
    {
        delmem_cd_op()(c_auxr_3d);
        c_auxr_3d = nullptr;
    }
}

template <>
void FFT_CUDA<double>::clear()
{
    this->cleanFFT();
    if (z_auxr_3d != nullptr)
    {
        delmem_zd_op()(z_auxr_3d);
        z_auxr_3d = nullptr;
        z_auxr_3d_capacity = 0;
    }
}

template <>
void FFT_CUDA<float>::fft3D_forward(std::complex<float>* in, std::complex<float>* out) const
{
    if (!this->c_handle)
    {
        CHECK_CUFFT(cufftPlan3d(&this->c_handle, this->nx, this->ny, this->nz, CUFFT_C2C));
    }
    CHECK_CUFFT(cufftExecC2C(this->c_handle,
                             reinterpret_cast<cufftComplex*>(in),
                             reinterpret_cast<cufftComplex*>(out),
                             CUFFT_FORWARD));
}

template <>
void FFT_CUDA<double>::fft3D_forward(std::complex<double>* in, std::complex<double>* out) const
{
    if (!this->z_handle)
    {
        CHECK_CUFFT(cufftPlan3d(&this->z_handle, this->nx, this->ny, this->nz, CUFFT_Z2Z));
    }
    CHECK_CUFFT(cufftExecZ2Z(this->z_handle,
                             reinterpret_cast<cufftDoubleComplex*>(in),
                             reinterpret_cast<cufftDoubleComplex*>(out),
                             CUFFT_FORWARD));
}

template <>
void FFT_CUDA<float>::fft3D_backward(std::complex<float>* in, std::complex<float>* out) const
{
    if (!this->c_handle)
    {
        CHECK_CUFFT(cufftPlan3d(&this->c_handle, this->nx, this->ny, this->nz, CUFFT_C2C));
    }
    CHECK_CUFFT(cufftExecC2C(this->c_handle,
                             reinterpret_cast<cufftComplex*>(in),
                             reinterpret_cast<cufftComplex*>(out),
                             CUFFT_INVERSE));
}

template <>
void FFT_CUDA<double>::fft3D_backward(std::complex<double>* in, std::complex<double>* out) const
{
    if (!this->z_handle)
    {
        CHECK_CUFFT(cufftPlan3d(&this->z_handle, this->nx, this->ny, this->nz, CUFFT_Z2Z));
    }
    CHECK_CUFFT(cufftExecZ2Z(this->z_handle,
                             reinterpret_cast<cufftDoubleComplex*>(in),
                             reinterpret_cast<cufftDoubleComplex*>(out),
                             CUFFT_INVERSE));
}

template <typename FPTYPE>
bool FFT_CUDA<FPTYPE>::use_pruned_gpu_fft() const
{
    return false;
}

template <>
bool FFT_CUDA<double>::use_pruned_gpu_fft() const
{
    // Small grids are dominated by the extra kernel launches. Sparse large
    // grids avoid the y-z transforms in the unoccupied middle x slabs.
    return this->xprime && !this->gamma_only && this->nproc == 1
           && this->nx * this->ny * this->nz >= pruned_fft_min_grid_points
           && pruned_fft_active_denominator * (this->lixy_index + 1 + this->nx - this->rixy_index)
                  <= pruned_fft_active_numerator * this->nx;
}

template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::init_gpu_pruned_plan() const
{
}

template <>
void FFT_CUDA<double>::init_gpu_pruned_plan() const
{
    if (this->z_x_handle)
    {
        return;
    }
    int x_size = this->nx;
    int yz_size[2] = {this->ny, this->nz};
    const int plane = this->ny * this->nz;
    CHECK_CUFFT(cufftPlanMany(&this->z_x_handle, 1, &x_size, &x_size, plane, 1,
                             &x_size, plane, 1, CUFFT_Z2Z, plane));
    CHECK_CUFFT(cufftPlanMany(&this->z_yz_left_handle, 2, yz_size, yz_size, 1, plane,
                             yz_size, 1, plane, CUFFT_Z2Z, this->lixy_index + 1));
    if (this->rixy_index < this->nx)
    {
        CHECK_CUFFT(cufftPlanMany(&this->z_yz_right_handle, 2, yz_size, yz_size, 1, plane,
                                 yz_size, 1, plane, CUFFT_Z2Z, this->nx - this->rixy_index));
    }
}

template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::fft3D_forward_pw(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const
{
    this->fft3D_forward(in, out);
}

template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::fft3D_backward_pw(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const
{
    this->fft3D_backward(in, out);
}

template <>
void FFT_CUDA<double>::fft3D_forward_pw(std::complex<double>* in, std::complex<double>* out) const
{
    if (!this->use_pruned_gpu_fft())
    {
        this->fft3D_forward(in, out);
        return;
    }
    this->init_gpu_pruned_plan();
    auto* input = reinterpret_cast<cufftDoubleComplex*>(in);
    auto* output = reinterpret_cast<cufftDoubleComplex*>(out);
    CHECK_CUFFT(cufftExecZ2Z(this->z_x_handle, input, output, CUFFT_FORWARD));
    CHECK_CUFFT(cufftExecZ2Z(this->z_yz_left_handle, output, output, CUFFT_FORWARD));
    if (this->z_yz_right_handle)
    {
        auto* right = output + this->rixy_index * this->ny * this->nz;
        CHECK_CUFFT(cufftExecZ2Z(this->z_yz_right_handle, right, right, CUFFT_FORWARD));
    }
}

template <>
void FFT_CUDA<double>::fft3D_backward_pw(std::complex<double>* in, std::complex<double>* out) const
{
    if (!this->use_pruned_gpu_fft())
    {
        this->fft3D_backward(in, out);
        return;
    }
    this->init_gpu_pruned_plan();
    auto* input = reinterpret_cast<cufftDoubleComplex*>(in);
    auto* output = reinterpret_cast<cufftDoubleComplex*>(out);
    CHECK_CUFFT(cufftExecZ2Z(this->z_yz_left_handle, input, input, CUFFT_INVERSE));
    if (this->z_yz_right_handle)
    {
        auto* right = input + this->rixy_index * this->ny * this->nz;
        CHECK_CUFFT(cufftExecZ2Z(this->z_yz_right_handle, right, right, CUFFT_INVERSE));
    }
    CHECK_CUFFT(cufftExecZ2Z(this->z_x_handle, input, output, CUFFT_INVERSE));
}

#if defined(__CUFFTDX)
template <typename FPTYPE>
int FFT_CUDA<FPTYPE>::fused_veff_grid_size() const
{
    return 0;
}
template <>
int FFT_CUDA<double>::fused_veff_grid_size() const
{
    return this->use_pruned_gpu_fft() && fused_fft(this->nx) != nullptr && fused_fft(this->ny) != nullptr
                   && fused_fft(this->nz) != nullptr
               ? (this->lixy_index + 1 + this->nx - this->rixy_index) * this->ny * this->nz
               : 0;
}
template <typename FPTYPE>
int FFT_CUDA<FPTYPE>::fused_veff_stick_size() const
{
    return 0;
}
template <>
int FFT_CUDA<double>::fused_veff_stick_size() const
{
    return this->fused_veff_grid_size() > 0
                   && 2 * this->ns_index <= (this->lixy_index + 1 + this->nx - this->rixy_index) * this->ny
               ? this->ns_index * this->nz
               : 0;
}
template <typename FPTYPE>
int FFT_CUDA<FPTYPE>::spinor_batch_size() const
{
    return 1;
}
template <>
int FFT_CUDA<double>::spinor_batch_size() const
{
    const int sticks = this->fused_veff_stick_size();
    if (sticks == 0) { return 1; }
    const size_t bytes = 2 * sizeof(std::complex<double>)
                         * (size_t(sticks) + this->fused_veff_grid_size());
    return std::max(1, std::min(max_spinor_batch, int(spinor_workspace_bytes / bytes)));
}
template <typename FPTYPE>
std::complex<double>* FFT_CUDA<FPTYPE>::get_spinor_stick_data(int) const
{
    return nullptr;
}
template <>
std::complex<double>* FFT_CUDA<double>::get_spinor_stick_data(int bands) const
{
    assert(bands == 1 || bands == 2 || bands == 4);
    const int size = std::max(this->nx * this->ny * this->nz, 2 * bands * this->fused_veff_stick_size());
    if (this->z_auxr_3d_capacity < size)
    {
        resmem_zd_op()(this->z_auxr_3d, size);
        this->z_auxr_3d_capacity = size;
    }
    return this->z_auxr_3d;
}
template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::transform_spinor_sticks(std::complex<double>*, std::complex<double>*, const int*, bool, int)
    const
{
}
template <>
void FFT_CUDA<double>::transform_spinor_sticks(std::complex<double>* up, std::complex<double>* down,
                                               const int* reverse, bool inverse, int bands) const
{
    const int stick_size = this->fused_veff_stick_size();
    assert(stick_size > 0 && reverse != nullptr);
    assert(bands == 1 || bands == 2 || bands == 4);
    auto& plan = this->z_stick_handle[bands / 2];
    if (!plan)
    {
        int length = this->nz;
        CHECK_CUFFT(cufftPlanMany(&plan, 1, &length, &length, 1, length,
                                 &length, 1, length, CUFFT_Z2Z, 2 * bands * this->ns_index));
    }
    std::complex<double>* slabs[2] = {up, down};
    auto* workspace = this->get_spinor_stick_data(bands);
    auto* data = reinterpret_cast<cufftDoubleComplex*>(workspace);
    if (inverse) { CHECK_CUFFT(cufftExecZ2Z(plan, data, data, CUFFT_INVERSE)); }
    for (int spin = 0; spin < 2; ++spin)
    {
        auto* sticks = workspace + spin * bands * stick_size;
        ModulePW::transform_stick_slabs(sticks,
                                        slabs[spin],
                                        reverse,
                                        this->nx,
                                        this->ny,
                                        this->nz,
                                        this->lixy_index + 1,
                                        this->rixy_index,
                                        inverse,
                                        bands,
                                        stick_size);
    }
    if (!inverse) { CHECK_CUFFT(cufftExecZ2Z(plan, data, data, CUFFT_FORWARD)); }
}
template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::transform_spinor_slabs(std::complex<double>*, std::complex<double>*, const int*, bool, int) const
{
}
template <>
void FFT_CUDA<double>::transform_spinor_slabs(std::complex<double>* up, std::complex<double>* down,
                                              const int* reverse, bool inverse, int bands) const
{
    const int stick_size = this->fused_veff_stick_size();
    assert(stick_size > 0 && reverse != nullptr);
    assert(bands >= 1 && bands <= this->spinor_batch_size());
    auto* workspace = this->get_spinor_stick_data(bands);
    std::complex<double>* slabs[2] = {up, down};
    for (int spin = 0; spin < 2; ++spin)
    {
        auto* sticks = workspace + spin * bands * stick_size;
        ModulePW::transform_stick_slabs(sticks, slabs[spin], reverse, this->nx, this->ny, this->nz,
                                        this->lixy_index + 1, this->rixy_index, inverse, bands, stick_size);
    }
}
template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::apply_spinor_potential(const double*, std::complex<double>*, std::complex<double>*,
                                              const int*, int) const
{
}
template <>
void FFT_CUDA<double>::apply_spinor_potential(const double* potential, std::complex<double>* up,
                                              std::complex<double>* down, const int* reverse, int bands) const
{
    assert(this->fused_veff_grid_size() > 0);
    if (reverse != nullptr)
    {
        this->transform_spinor_sticks(up, down, reverse, true, bands);
        ModulePW::fused_veff_x(this->nx,
                               this->ny * this->nz,
                               this->lixy_index + 1,
                               this->rixy_index,
                               potential,
                               up,
                               down,
                               bands);
        this->transform_spinor_sticks(up, down, reverse, false, bands);
        return;
    }
    assert(bands == 1);
    this->init_gpu_pruned_plan();
    const int plane = this->ny * this->nz;
    const int left = this->lixy_index + 1;
    auto transform_yz = [&](std::complex<double>* grid, int direction)
    {
        auto* data = reinterpret_cast<cufftDoubleComplex*>(grid);
        CHECK_CUFFT(cufftExecZ2Z(this->z_yz_left_handle, data, data, direction));
        if (this->z_yz_right_handle)
        {
            data += left * plane;
            CHECK_CUFFT(cufftExecZ2Z(this->z_yz_right_handle, data, data, direction));
        }
    };
    transform_yz(up, CUFFT_INVERSE);
    transform_yz(down, CUFFT_INVERSE);
    ModulePW::fused_veff_x(this->nx, plane, left, this->rixy_index, potential, up, down);
    transform_yz(up, CUFFT_FORWARD);
    transform_yz(down, CUFFT_FORWARD);
}
#endif

template <>
std::complex<float>* FFT_CUDA<float>::get_auxr_3d_data() const
{
    if (this->c_auxr_3d == nullptr)
    {
        resmem_cd_op()(this->c_auxr_3d, this->nx * this->ny * this->nz);
    }
    return this->c_auxr_3d;
}

template <>
std::complex<double>* FFT_CUDA<double>::get_auxr_3d_data() const
{
    if (this->z_auxr_3d == nullptr)
    {
        this->z_auxr_3d_capacity = this->nx * this->ny * this->nz;
        resmem_zd_op()(this->z_auxr_3d, this->z_auxr_3d_capacity);
    }
    return this->z_auxr_3d;
}

template FFT_CUDA<float>::FFT_CUDA();
template FFT_CUDA<float>::~FFT_CUDA();
template FFT_CUDA<double>::FFT_CUDA();
template FFT_CUDA<double>::~FFT_CUDA();
} // namespace ModuleBase
