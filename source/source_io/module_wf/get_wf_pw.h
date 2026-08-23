#ifndef GET_WF_PW_H
#define GET_WF_PW_H

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_base/parallel_comm.h"
#include "source_io/module_output/band_parallel_output.h"

namespace ModuleIO
{
/**
 * @brief Write real-space norms and complex components of selected PW states.
 *
 * Every selected band and k-point is written independently. Scalar or collinear
 * states produce one component, while an nspin=4 spinor produces a combined norm
 * cube and separate real/imaginary cubes for its up and down components.
 */
template <typename Device>
void get_wf_pw(const std::vector<int>& out_wfc_norm,
               const std::vector<int>& out_wfc_re_im,
               const int nspin,
               const int global_nbands,
               UnitCell* ucell,
               const psi::Psi<std::complex<double>, Device>* kspw_psi,
               const ModulePW::PW_Basis_K* pw_wfc,
               const ModulePW::PW_Basis* pw_rho,
               const ModulePW::PW_Basis* pw_rhod,
               const Parallel_Grid& pgrid,
               const std::string& global_out_dir,
               const K_Vectors& kv)
{
    const int nks = kv.get_nks();       // current process pool k-point count
    const int nkstot = kv.get_nkstot(); // total k-point count
    // The output masks address global bands, but BPCG stores only the current band group's
    // contiguous shard in Psi. The collective layout supplies the required mapping.
    const BandParallelLayout band_layout(kspw_psi->get_nbands(), global_nbands);

    const int nks_without_spin = nspin == 2 ? nkstot / 2 : nkstot;
    const int smooth_nrxx = pw_wfc->nrxx;
    const int dense_nrxx = pw_rhod->nrxx;
    // Avoid an extra forward and backward FFT unless a distinct dense grid is in use.
    const bool needs_interpolation = pw_rhod != pw_rho;

    // The two INPUT vectors are independent band masks for norm and Re/Im output.
    std::vector<int> bands_picked_norm(global_nbands, 0);
    std::vector<int> bands_picked_re_im(global_nbands, 0);

    if (static_cast<int>(out_wfc_norm.size()) > global_nbands || static_cast<int>(out_wfc_re_im.size()) > global_nbands)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::get_wf_pw",
                                 "The number of bands specified by `out_wfc_norm` or `out_wfc_re_im` in the "
                                 "INPUT file exceeds `nbands`!");
    }

    for (int value: out_wfc_norm)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::get_wf_pw",
                                     "The elements of `out_wfc_norm` must be either 0 or 1. "
                                     "Invalid values found!");
        }
    }
    for (int value: out_wfc_re_im)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::get_wf_pw",
                                     "The elements of `out_wfc_re_im` must be either 0 or 1. "
                                     "Invalid values found!");
        }
    }

    int length = std::min(static_cast<int>(out_wfc_norm.size()), global_nbands);
    for (int i = 0; i < length; ++i)
    {
        bands_picked_norm[i] = static_cast<int>(out_wfc_norm[i]);
    }
    length = std::min(static_cast<int>(out_wfc_re_im.size()), global_nbands);
    for (int i = 0; i < length; ++i)
    {
        bands_picked_re_im[i] = static_cast<int>(out_wfc_re_im[i]);
    }

    // Map the wavefunction backend type to the tensor library's device type.
    using ContainerDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device_type = ct::DeviceTypeToEnum<ContainerDevice>::value;
    const bool is_cpu = device_type == ct::DeviceType::CpuDevice;
    const bool is_spinor = nspin == 4;
    // Spinor coefficients store the up and down blocks consecutively.
    const int npwx = kspw_psi->get_nbasis() / (is_spinor ? 2 : 1);

    // Zero-length tensors avoid allocating buffers that a given backend or spin mode never uses.
    ct::Tensor wfcr_up_smooth(ct::DataType::DT_COMPLEX_DOUBLE, device_type, ct::TensorShape({smooth_nrxx}));
    ct::Tensor wfcr_down_smooth(ct::DataType::DT_COMPLEX_DOUBLE, device_type, ct::TensorShape({is_spinor ? smooth_nrxx : 0}));
    ct::Tensor wfcr_up_smooth_host(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, ct::TensorShape({is_cpu ? 0 : smooth_nrxx}));
    ct::Tensor wfcr_down_smooth_host(ct::DataType::DT_COMPLEX_DOUBLE,
                                     ct::DeviceType::CpuDevice,
                                     ct::TensorShape({!is_cpu && is_spinor ? smooth_nrxx : 0}));
    ct::Tensor wfcr_up_dense_host(ct::DataType::DT_COMPLEX_DOUBLE,
                                  ct::DeviceType::CpuDevice,
                                  ct::TensorShape({needs_interpolation ? dense_nrxx : 0}));
    ct::Tensor wfcr_down_dense_host(ct::DataType::DT_COMPLEX_DOUBLE,
                                    ct::DeviceType::CpuDevice,
                                    ct::TensorShape({needs_interpolation && is_spinor ? dense_nrxx : 0}));
    ct::Tensor reciprocal_buffer_host(ct::DataType::DT_COMPLEX_DOUBLE,
                                      ct::DeviceType::CpuDevice,
                                      ct::TensorShape({needs_interpolation ? pw_rhod->npw : 0}));

    // Capture the shared bases and scratch buffers by reference. The returned pointer is owned by
    // one of the tensor arguments and remains valid until that tensor is reused or destroyed.
    auto transform_wfc = [&](const std::complex<double>* coefficients,
                             const int ik,
                             ct::Tensor& smooth,
                             ct::Tensor& smooth_host,
                             ct::Tensor& dense_host) -> const std::complex<double>* {
        // Perform the wavefunction FFT on its native device, then expose host data for cube output.
        pw_wfc->template recip_to_real<std::complex<double>, Device>(coefficients, smooth.data<std::complex<double>>(), ik);
        const std::complex<double>* smooth_data = smooth.data<std::complex<double>>();
        if (!is_cpu)
        {
            ct::kernels::synchronize_memory<std::complex<double>, ct::DEVICE_CPU, ContainerDevice>()(
                smooth_host.data<std::complex<double>>(),
                smooth_data,
                smooth_nrxx);
            smooth_data = smooth_host.data<std::complex<double>>();
        }
        if (!needs_interpolation)
        {
            return smooth_data;
        }

        // Zero padding in reciprocal space transfers the smooth-grid field to the dense rho grid.
        reciprocal_buffer_host.zero();
        pw_rho->real2recip(smooth_data, reciprocal_buffer_host.data<std::complex<double>>());
        pw_rhod->recip2real(reciprocal_buffer_host.data<std::complex<double>>(), dense_host.data<std::complex<double>>());
        return dense_host.data<std::complex<double>>();
    };

    // These buffers hold one rank-local dense-grid slab after resolving band ownership;
    // the real-space grid itself remains distributed over POOL_WORLD.
    std::vector<std::complex<double>> wfcr_up_global(dense_nrxx);
    std::vector<std::complex<double>> wfcr_down_global(is_spinor ? dense_nrxx : 0);
    // The owner alone indexes local Psi and performs the FFT. BP_WORLD then gives every
    // band group the same slab without gathering the complete plane-wave wavefunction.
    // Callers on all ranks must preserve an identical collective order.
    auto transform_global_band = [&](const int global_band,
                                     const int basis_offset,
                                     const int ik,
                                     ct::Tensor& smooth,
                                     ct::Tensor& smooth_host,
                                     ct::Tensor& dense_host,
                                     std::vector<std::complex<double>>& global_wfcr) -> const std::complex<double>* {
        const int owner = band_layout.owner_group(global_band);
        if (band_layout.band_group() == owner)
        {
            const int local_band = band_layout.local_index(global_band);
            kspw_psi->fix_k(ik);
            const std::complex<double>* owner_wfcr
                = transform_wfc(&kspw_psi[0](local_band, basis_offset), ik, smooth, smooth_host, dense_host);
            std::copy(owner_wfcr, owner_wfcr + dense_nrxx, global_wfcr.begin());
        }
#ifdef __MPI
        MPI_Bcast(global_wfcr.data(), dense_nrxx, MPI_DOUBLE_COMPLEX, owner, BP_WORLD);
#endif
        return global_wfcr.data();
    };

    // Norm cubes are phase invariant; for spinors they contain sqrt(|up|^2 + |down|^2).
    // Iterating global bands on all groups keeps owner broadcasts synchronized.
    std::vector<std::vector<double>> rho_band_norm(nspin, std::vector<double>(dense_nrxx));
    for (int ib = 0; ib < global_nbands; ++ib)
    {
        if (!bands_picked_norm[ib])
        {
            continue;
        }

        for (int is = 0; is < nspin; ++is)
        {
            std::fill(rho_band_norm[is].begin(), rho_band_norm[is].end(), 0.0);
        }
        for (int ik = 0; ik < nks; ++ik)
        {
            const int ikstot = kv.ik2iktot[ik];
            const int spin_index = kv.isk[ik];
            // In collinear calculations the two spin channels share the same k-point numbering.
            const int k_number = ikstot % nks_without_spin + 1;

            const std::complex<double>* wfcr_up_host_data
                = transform_global_band(ib, 0, ik, wfcr_up_smooth, wfcr_up_smooth_host, wfcr_up_dense_host, wfcr_up_global);
            const std::complex<double>* wfcr_down_host_data = nullptr;
            if (is_spinor)
            {
                wfcr_down_host_data
                    = transform_global_band(ib, npwx, ik, wfcr_down_smooth, wfcr_down_smooth_host, wfcr_down_dense_host, wfcr_down_global);
            }

            const double spin_degeneracy = nspin == 1 ? 2.0 : 1.0;
            const double scale = std::sqrt(spin_degeneracy / ucell->omega);
            for (int ir = 0; ir < dense_nrxx; ++ir)
            {
                const double norm = is_spinor ? std::sqrt(std::norm(wfcr_up_host_data[ir]) + std::norm(wfcr_down_host_data[ir]))
                                              : std::abs(wfcr_up_host_data[ir]);
                rho_band_norm[spin_index][ir] = norm * scale;
            }

            std::stringstream ss_file;
            ss_file << global_out_dir << "wfi" << ib + 1 << "s" << spin_index + 1 << "k" << k_number << ".cube";
            ModuleIO::write_vdata_palgrid(pgrid,
                                          rho_band_norm[spin_index].data(),
                                          spin_index,
                                          nspin,
                                          0,
                                          ss_file.str(),
                                          0.0,
                                          ucell,
                                          11,
                                          0,
                                          false,
                                          true);
        }
    }

    // Re/Im cubes retain phase information and therefore write both spinor components separately.
    // This path deliberately reuses the same ownership rule as the norm path.
    std::vector<std::vector<double>> rho_band_re(nspin, std::vector<double>(dense_nrxx));
    std::vector<std::vector<double>> rho_band_im(nspin, std::vector<double>(dense_nrxx));
    for (int ib = 0; ib < global_nbands; ++ib)
    {
        if (!bands_picked_re_im[ib])
        {
            continue;
        }

        for (int is = 0; is < nspin; ++is)
        {
            std::fill(rho_band_re[is].begin(), rho_band_re[is].end(), 0.0);
            std::fill(rho_band_im[is].begin(), rho_band_im[is].end(), 0.0);
        }
        for (int ik = 0; ik < nks; ++ik)
        {
            const int ikstot = kv.ik2iktot[ik];
            const int spin_index = kv.isk[ik];
            const int k_number = ikstot % nks_without_spin + 1;

            const std::complex<double>* wfcr_up_host_data
                = transform_global_band(ib, 0, ik, wfcr_up_smooth, wfcr_up_smooth_host, wfcr_up_dense_host, wfcr_up_global);
            const std::complex<double>* wfcr_down_host_data = nullptr;
            if (is_spinor)
            {
                wfcr_down_host_data
                    = transform_global_band(ib, npwx, ik, wfcr_down_smooth, wfcr_down_smooth_host, wfcr_down_dense_host, wfcr_down_global);
            }

            const double spin_degeneracy = nspin == 1 ? 2.0 : 1.0;
            const double scale = std::sqrt(spin_degeneracy / ucell->omega);
            // For scalar/collinear states, select kv.isk[ik]; for spinors, emit both up and down.
            const int component_begin = is_spinor ? 0 : spin_index;
            const int component_end = is_spinor ? 2 : spin_index + 1;
            for (int component = component_begin; component < component_end; ++component)
            {
                const std::complex<double>* component_data = is_spinor && component == 1 ? wfcr_down_host_data : wfcr_up_host_data;
                for (int ir = 0; ir < dense_nrxx; ++ir)
                {
                    rho_band_re[component][ir] = std::real(component_data[ir]) * scale;
                    rho_band_im[component][ir] = std::imag(component_data[ir]) * scale;
                }

                std::stringstream ss_real;
                ss_real << global_out_dir << "wfi" << ib + 1 << "s" << component + 1 << "k" << k_number << "re.cube";
                ModuleIO::write_vdata_palgrid(pgrid,
                                              rho_band_re[component].data(),
                                              component,
                                              nspin,
                                              0,
                                              ss_real.str(),
                                              0.0,
                                              ucell,
                                              11,
                                              0,
                                              false,
                                              true);

                std::stringstream ss_imag;
                ss_imag << global_out_dir << "wfi" << ib + 1 << "s" << component + 1 << "k" << k_number << "im.cube";
                ModuleIO::write_vdata_palgrid(pgrid,
                                              rho_band_im[component].data(),
                                              component,
                                              nspin,
                                              0,
                                              ss_imag.str(),
                                              0.0,
                                              ucell,
                                              11,
                                              0,
                                              false,
                                              true);
            }
        }
    }
}
} // namespace ModuleIO

#endif // GET_WF_PW_H
