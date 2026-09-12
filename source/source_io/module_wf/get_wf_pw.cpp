#include "source_io/module_wf/get_wf_pw.h"

#include "source_base/constants.h"
#include "source_base/module_container/ATen/core/tensor.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_base/tool_quit.h"
#include "source_io/module_output/cube_io.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <type_traits>

namespace ModuleIO
{
// This nested class owns scratch storage for one begin() call and can access the output object's private data.
template <typename T, typename Device>
class Get_wf_pw<T, Device>::Workspace
{
  public:
    // typename marks a type selected from the Device-dependent mapping.
    using ContainerDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device_type = ct::DeviceTypeToEnum<ContainerDevice>::value;
    const bool is_cpu = device_type == ct::DeviceType::CpuDevice;
    const bool needs_host_copy = !is_cpu || !std::is_same<T, std::complex<double>>::value;
    const bool is_spinor;
    const bool needs_interpolation;
    const int smooth_nrxx;
    const int dense_nrxx;
    const int npwx;
    ct::Tensor smooth[2];
    ct::Tensor smooth_host[2];
    ct::Tensor dense_host[2];
    ct::Tensor reciprocal_host;
    std::vector<std::complex<double>> wfcr[2];

    // Norm and Re output use separate traversals and share the same values buffer.
    std::vector<std::vector<double>> values;
    std::vector<std::vector<double>> imag;
    std::vector<std::complex<double>> phase;

    // Brace initialization constructs both spinor slots; unused device/grid buffers have zero size.
    explicit Workspace(const Get_wf_pw& output)
        : is_spinor(output.nspin_ == 4), needs_interpolation(&output.pw_rhod_ != &output.pw_rho_), smooth_nrxx(output.pw_wfc_.nrxx),
          dense_nrxx(output.pw_rhod_.nrxx), npwx(output.psi_.get_nbasis() / (is_spinor ? 2 : 1)),
          smooth{ct::Tensor(ct::DataTypeToEnum<T>::value, device_type, ct::TensorShape({smooth_nrxx})),
                 ct::Tensor(ct::DataTypeToEnum<T>::value, device_type, ct::TensorShape({is_spinor ? smooth_nrxx : 0}))},
          smooth_host{
              ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({needs_host_copy ? smooth_nrxx : 0})),
              ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                         ct::DeviceType::CpuDevice,
                         ct::TensorShape({needs_host_copy && is_spinor ? smooth_nrxx : 0}))},
          dense_host{ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                ct::DeviceType::CpuDevice,
                                ct::TensorShape({needs_interpolation ? dense_nrxx : 0})),
                     ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                ct::DeviceType::CpuDevice,
                                ct::TensorShape({needs_interpolation && is_spinor ? dense_nrxx : 0}))},
          reciprocal_host(ct::DataType::DT_COMPLEX_DOUBLE,
                          ct::DeviceType::CpuDevice,
                          ct::TensorShape({needs_interpolation ? output.pw_rhod_.npw : 0})),
          wfcr{std::vector<std::complex<double>>(dense_nrxx), std::vector<std::complex<double>>(is_spinor ? dense_nrxx : 0)},
          values(output.nspin_, std::vector<double>(dense_nrxx)), imag(output.nspin_, std::vector<double>(dense_nrxx)), phase(dense_nrxx)
    {
    }
};

template <typename T, typename Device>
Get_wf_pw<T, Device>::Get_wf_pw(const psi::Psi<T, Device>& psi,
                                const ModulePW::PW_Basis_K& pw_wfc,
                                const ModulePW::PW_Basis& pw_rho,
                                const ModulePW::PW_Basis& pw_rhod,
                                const int nspin,
                                const int global_nbands)
    : psi_(psi), pw_wfc_(pw_wfc), pw_rho_(pw_rho), pw_rhod_(pw_rhod), nspin_(nspin), global_nbands_(global_nbands)
{
}

template <typename T, typename Device>
void Get_wf_pw<T, Device>::begin(const UnitCell& ucell,
                                 const Parallel_Grid& pgrid,
                                 const K_Vectors& kv,
                                 const std::vector<int>& out_wfc_norm,
                                 const std::vector<int>& out_wfc_re_im,
                                 const std::string& global_out_dir) const
{
    // Resolve global band ownership collectively before validating the selection.
    const BandParallelLayout layout(psi_.get_nbands(), global_nbands_);
    if (static_cast<int>(out_wfc_norm.size()) > global_nbands_ || static_cast<int>(out_wfc_re_im.size()) > global_nbands_)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::get_wf_pw",
                                 "The number of bands specified by `out_wfc_norm` or `out_wfc_re_im` in the INPUT file exceeds `nbands`!");
    }
    const std::vector<int> norm_mask = select_bands(out_wfc_norm, "out_wfc_norm");
    const std::vector<int> re_im_mask = select_bands(out_wfc_re_im, "out_wfc_re_im");
    Workspace work(*this);
    // Keep norm and Re/Im traversals separate to preserve the FFT and collective order.
    for (int band = 0; band < global_nbands_; ++band)
    {
        if (!norm_mask[band])
        {
            continue;
        }
        for (int is = 0; is < nspin_; ++is)
        {
            std::fill(work.values[is].begin(), work.values[is].end(), 0.0);
        }
        write_norm(band, ucell, pgrid, kv, global_out_dir, layout, &work);
    }
    for (int band = 0; band < global_nbands_; ++band)
    {
        if (!re_im_mask[band])
        {
            continue;
        }
        for (int is = 0; is < nspin_; ++is)
        {
            std::fill(work.values[is].begin(), work.values[is].end(), 0.0);
            std::fill(work.imag[is].begin(), work.imag[is].end(), 0.0);
        }
        write_complex(band, ucell, pgrid, kv, global_out_dir, layout, &work);
    }
}

template <typename T, typename Device>
std::vector<int> Get_wf_pw<T, Device>::select_bands(const std::vector<int>& selection, const std::string& parameter_name) const
{
    // begin() checks all selection lengths first; omitted bands remain unselected.
    std::vector<int> band_mask(global_nbands_, 0);
    for (int value: selection)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::get_wf_pw",
                                     "The elements of `" + parameter_name + "` must be either 0 or 1. Invalid values found!");
        }
    }
    std::copy(selection.begin(), selection.end(), band_mask.begin());
    return band_mask;
}

template <typename T, typename Device>
void Get_wf_pw<T, Device>::transform_band(const int global_band, const int ik, const BandParallelLayout& layout, Workspace* work) const
{
    const int owner = layout.owner_group(global_band);
    // All band groups must visit the same band/k/component sequence. Only the
    // owner indexes Psi; the broadcast replicates a slab, not the entire grid.
    for (int component = 0; component < (work->is_spinor ? 2 : 1); ++component)
    {
        if (layout.band_group() == owner)
        {
            const int local_band = layout.local_index(global_band);
            psi_.fix_k(ik);
            // Spinor coefficients occupy two consecutive blocks with stride npwx.
            const std::complex<double>* owner_wfcr = transform_wfc(&psi_(local_band, component * work->npwx), ik, component, work);
            std::copy(owner_wfcr, owner_wfcr + work->dense_nrxx, work->wfcr[component].begin());
        }
#ifdef __MPI
        Parallel_Common::bcast_data(work->wfcr[component].data(), work->dense_nrxx, BP_WORLD, owner);
#endif
    }
}

template <typename T, typename Device>
const std::complex<double>* Get_wf_pw<T, Device>::transform_wfc(const T* coefficients,
                                                                const int ik,
                                                                const int component,
                                                                Workspace* work) const
{
    // Reconstruct the periodic part u_nk(r) using the solver's precision and device.
    // .template identifies a member template when the object type depends on T or Device.
    pw_wfc_.template recip_to_real<T, Device>(coefficients, work->smooth[component].template data<T>(), ik);
    const std::complex<double>* smooth_data = nullptr;
    if (work->needs_host_copy)
    {
        // Convert only this state's FFT result to host double for grid processing.
        base_device::memory::cast_memory_op<std::complex<double>, T, base_device::DEVICE_CPU, Device>()(
            work->smooth_host[component].template data<std::complex<double>>(),
            work->smooth[component].template data<T>(),
            work->smooth_nrxx);
        smooth_data = work->smooth_host[component].template data<std::complex<double>>();
    }
    else
    {
        // CPU double output can use the FFT buffer directly.
        smooth_data = work->smooth[component].template data<std::complex<double>>();
    }
    if (!work->needs_interpolation)
    {
        return smooth_data;
    }

    // Keep the smooth-grid Fourier coefficients and pad the extra dense-grid coefficients with zero.
    // The inverse FFT samples the same periodic wavefunction on the dense rho grid.
    work->reciprocal_host.zero();
    pw_rho_.real2recip(smooth_data, work->reciprocal_host.template data<std::complex<double>>());
    pw_rhod_.recip2real(work->reciprocal_host.template data<std::complex<double>>(),
                        work->dense_host[component].template data<std::complex<double>>());
    return work->dense_host[component].template data<std::complex<double>>();
}

template <typename T, typename Device>
void Get_wf_pw<T, Device>::write_norm(const int band,
                                      const UnitCell& ucell,
                                      const Parallel_Grid& pgrid,
                                      const K_Vectors& kv,
                                      const std::string& out_dir,
                                      const BandParallelLayout& layout,
                                      Workspace* work) const
{
    // Collinear spin channels share the same physical k-point numbering in file names.
    const int nks_without_spin = nspin_ == 2 ? kv.get_nkstot() / 2 : kv.get_nkstot();
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        const int spin_index = kv.isk[ik];
        const int k_number = kv.ik2iktot[ik] % nks_without_spin + 1;
        transform_band(band, ik, layout, work);
        // The solver supplies L2-normalized NC states or S-normalized USPP states.
        // Wavefunction amplitudes carry the inverse square root of the cell volume.
        const double scale = std::sqrt(1.0 / ucell.omega);
        calc_norm(spin_index, scale, work);
        write_cube(band, spin_index, k_number, "", ucell, pgrid, out_dir, work->values[spin_index]);
    }
}

template <typename T, typename Device>
void Get_wf_pw<T, Device>::write_complex(const int band,
                                         const UnitCell& ucell,
                                         const Parallel_Grid& pgrid,
                                         const K_Vectors& kv,
                                         const std::string& out_dir,
                                         const BandParallelLayout& layout,
                                         Workspace* work) const
{
    // Collinear spin channels share the same physical k-point numbering in file names.
    const int nks_without_spin = nspin_ == 2 ? kv.get_nkstot() / 2 : kv.get_nkstot();
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        const int spin_index = kv.isk[ik];
        const int k_number = kv.ik2iktot[ik] % nks_without_spin + 1;
        transform_band(band, ik, layout, work);
        // The solver supplies L2-normalized NC states or S-normalized USPP states.
        // Wavefunction amplitudes carry the inverse square root of the cell volume.
        const double scale = std::sqrt(1.0 / ucell.omega);
        calc_phase(ik, kv, &work->phase);
        // Scalar/collinear output selects isk; spinors emit both components.
        const int component_begin = work->is_spinor ? 0 : spin_index;
        const int component_end = work->is_spinor ? 2 : spin_index + 1;
        for (int component = component_begin; component < component_end; ++component)
        {
            const int wfc_component = work->is_spinor && component == 1 ? 1 : 0;
            calc_component(work->wfcr[wfc_component], work->phase, scale, &work->values[component], &work->imag[component]);
            write_cube(band, component, k_number, "re", ucell, pgrid, out_dir, work->values[component]);
            write_cube(band, component, k_number, "im", ucell, pgrid, out_dir, work->imag[component]);
        }
    }
}

template <typename T, typename Device>
void Get_wf_pw<T, Device>::calc_norm(const int spin_index, const double scale, Workspace* work) const
{
    // Output the modulus, not its square; a spinor combines the squared moduli of both components.
    // The Bloch phase has unit modulus and therefore does not enter this field.
    for (int ir = 0; ir < work->dense_nrxx; ++ir)
    {
        const double norm
            = work->is_spinor ? std::sqrt(std::norm(work->wfcr[0][ir]) + std::norm(work->wfcr[1][ir])) : std::abs(work->wfcr[0][ir]);
        work->values[spin_index][ir] = norm * scale;
    }
}

template <typename T, typename Device>
void Get_wf_pw<T, Device>::calc_phase(const int ik, const K_Vectors& kv, std::vector<std::complex<double>>* phase) const
{
    // Build exp(i k.r) from fractional k-point and grid coordinates to recover the Bloch state.
    // The local slab is stored as [x][y][local_z]; startz_current restores the global z index.
    for (int ir = 0; ir < pw_rhod_.nrxx; ++ir)
    {
        const int ix = ir / (pw_rhod_.ny * pw_rhod_.nplane);
        const int iy = ir / pw_rhod_.nplane % pw_rhod_.ny;
        const int iz = ir % pw_rhod_.nplane + pw_rhod_.startz_current;
        const double phase_argument
            = ModuleBase::TWO_PI
              * (kv.kvec_d[ik].x * static_cast<double>(ix) / pw_rhod_.nx + kv.kvec_d[ik].y * static_cast<double>(iy) / pw_rhod_.ny
                 + kv.kvec_d[ik].z * static_cast<double>(iz) / pw_rhod_.nz);
        (*phase)[ir] = std::exp(std::complex<double>(0.0, phase_argument));
    }
}

template <typename T, typename Device>
void Get_wf_pw<T, Device>::calc_component(const std::vector<std::complex<double>>& component,
                                          const std::vector<std::complex<double>>& phase,
                                          const double scale,
                                          std::vector<double>* real,
                                          std::vector<double>* imag) const
{
    // Restore psi_nk(r) = exp(i k.r) u_nk(r) before separating and normalizing Re/Im.
    for (int ir = 0; ir < pw_rhod_.nrxx; ++ir)
    {
        const std::complex<double> bloch_wfc = component[ir] * phase[ir];
        (*real)[ir] = std::real(bloch_wfc) * scale;
        (*imag)[ir] = std::imag(bloch_wfc) * scale;
    }
}

template <typename T, typename Device>
void Get_wf_pw<T, Device>::write_cube(const int band,
                                      const int component,
                                      const int k_number,
                                      const std::string& part,
                                      const UnitCell& ucell,
                                      const Parallel_Grid& pgrid,
                                      const std::string& out_dir,
                                      const std::vector<double>& values) const
{
    std::stringstream filename;
    filename << out_dir << "wfi" << band + 1 << "s" << component + 1 << "k" << k_number << part << ".cube";
    ModuleIO::write_vdata_palgrid(pgrid, values.data(), component, nspin_, 0, filename.str(), 0.0, &ucell, 11, 0, false, true);
}

// Explicit instantiation emits both precisions for each supported device from this .cpp file.
template class Get_wf_pw<std::complex<float>, base_device::DEVICE_CPU>;
template class Get_wf_pw<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class Get_wf_pw<std::complex<float>, base_device::DEVICE_GPU>;
template class Get_wf_pw<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleIO
