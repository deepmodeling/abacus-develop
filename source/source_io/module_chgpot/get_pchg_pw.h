#ifndef GET_PCHG_PW_H
#define GET_PCHG_PW_H

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_estate/module_charge/symm_rho.h"
#include "source_io/module_output/cube_io.h"

namespace ModuleIO
{
/**
 * @brief Write band-resolved partial charges from plane-wave coefficients.
 *
 * The selected bands are transformed to the dense real-space grid. Depending on
 * @p if_separate_k, the function either writes one cube per k-point or sums the
 * k-point contributions, restores symmetry, and writes one cube per spin/charge
 * component. For spinors, the four components are charge, m_x, m_y, and m_z.
 */
template <typename Device>
void get_pchg_pw(const std::vector<int>& out_pchg,
                 const int nspin,
                 UnitCell* ucell,
                 const psi::Psi<std::complex<double>, Device>* kspw_psi,
                 const ModulePW::PW_Basis* pw_rho,
                 const ModulePW::PW_Basis* pw_rhod,
                 const ModulePW::PW_Basis_K* pw_wfc,
                 const Parallel_Grid& pgrid,
                 const std::string& global_out_dir,
                 const bool if_separate_k,
                 const bool noncolin,
                 const K_Vectors& kv)
{
    const int nks = kv.get_nks();       // current process pool k-point count
    const int nkstot = kv.get_nkstot(); // total k-point count
    const int nbands = kspw_psi->get_nbands();

    const int nks_without_spin = nspin == 2 ? nkstot / 2 : nkstot;
    const int smooth_nrxx = pw_wfc->nrxx;
    const int dense_nrxx = pw_rhod->nrxx;
    // Avoid an extra forward and backward FFT unless a distinct dense grid is in use.
    const bool needs_interpolation = pw_rhod != pw_rho;

    // Expand the INPUT selection into a fixed-size mask indexed directly by band.
    std::vector<int> bands_picked(nbands, 0);
    if (static_cast<int>(out_pchg.size()) > nbands)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::get_pchg_pw",
                                 "The number of bands specified by `out_pchg` in the "
                                 "INPUT file exceeds `nbands`!");
    }
    for (int value: out_pchg)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::get_pchg_pw",
                                     "The elements of `out_pchg` must be either 0 or 1. "
                                     "Invalid values found!");
        }
    }
    const int length = std::min(static_cast<int>(out_pchg.size()), nbands);
    for (int i = 0; i < length; ++i)
    {
        bands_picked[i] = static_cast<int>(out_pchg[i]);
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

    std::vector<std::vector<double>> rho_band(nspin, std::vector<double>(dense_nrxx));
    // Convert a two-component spinor into the Pauli-basis fields (rho, m_x, m_y, m_z).
    // Per-k output overwrites the fields, whereas k-summed output accumulates weighted fields.
    auto accumulate_spinor_density
        = [&](const std::complex<double>* up, const std::complex<double>* down, const double weight, const bool accumulate) {
              for (int ir = 0; ir < dense_nrxx; ++ir)
              {
                  const double up_norm = std::norm(up[ir]);
                  const double down_norm = std::norm(down[ir]);
                  const double rho0 = (up_norm + down_norm) * weight;
                  const double mx = 2.0 * (up[ir].real() * down[ir].real() + up[ir].imag() * down[ir].imag()) * weight;
                  const double my = 2.0 * (up[ir].real() * down[ir].imag() - down[ir].real() * up[ir].imag()) * weight;
                  const double mz = (up_norm - down_norm) * weight;
                  if (accumulate)
                  {
                      rho_band[0][ir] += rho0;
                      rho_band[1][ir] += noncolin ? mx : 0.0;
                      rho_band[2][ir] += noncolin ? my : 0.0;
                      rho_band[3][ir] += mz;
                  }
                  else
                  {
                      rho_band[0][ir] = rho0;
                      rho_band[1][ir] = noncolin ? mx : 0.0;
                      rho_band[2][ir] = noncolin ? my : 0.0;
                      rho_band[3][ir] = mz;
                  }
              }
          };

    for (int ib = 0; ib < nbands; ++ib)
    {
        if (!bands_picked[ib])
        {
            continue;
        }

        for (int is = 0; is < nspin; ++is)
        {
            std::fill(rho_band[is].begin(), rho_band[is].end(), 0.0);
        }

        if (if_separate_k)
        {
            // Preserve each Bloch state's contribution; no Brillouin-zone weight is applied here.
            for (int ik = 0; ik < nks; ++ik)
            {
                const int ikstot = kv.ik2iktot[ik];
                const int spin_index = kv.isk[ik];
                // In collinear calculations the two spin channels share the same k-point numbering.
                const int k_number = ikstot % nks_without_spin + 1;

                kspw_psi->fix_k(ik);
                const std::complex<double>* wfcr_up
                    = transform_wfc(&kspw_psi[0](ib, 0), ik, wfcr_up_smooth, wfcr_up_smooth_host, wfcr_up_dense_host);
                const std::complex<double>* wfcr_up_host_data = wfcr_up;
                const std::complex<double>* wfcr_down_host_data = nullptr;
                if (is_spinor)
                {
                    const std::complex<double>* wfcr_down
                        = transform_wfc(&kspw_psi[0](ib, npwx), ik, wfcr_down_smooth, wfcr_down_smooth_host, wfcr_down_dense_host);
                    wfcr_down_host_data = wfcr_down;
                }

                const double spin_degeneracy = nspin == 1 ? 2.0 : 1.0;
                const double weight = spin_degeneracy / ucell->omega;
                if (is_spinor)
                {
                    accumulate_spinor_density(wfcr_up_host_data, wfcr_down_host_data, weight, false);
                }
                else
                {
                    for (int ir = 0; ir < dense_nrxx; ++ir)
                    {
                        rho_band[spin_index][ir] = std::norm(wfcr_up_host_data[ir]) * weight;
                    }
                }

                const int component_begin = is_spinor ? 0 : spin_index;
                const int component_end = is_spinor ? 4 : spin_index + 1;
                for (int component = component_begin; component < component_end; ++component)
                {
                    std::stringstream ssc;
                    ssc << global_out_dir << "pchgi" << ib + 1 << "s" << component + 1 << "k" << k_number << ".cube";
                    ModuleIO::write_vdata_palgrid(pgrid,
                                                  rho_band[component].data(),
                                                  component,
                                                  nspin,
                                                  0,
                                                  ssc.str(),
                                                  0.0,
                                                  ucell,
                                                  11,
                                                  0,
                                                  false,
                                                  true);
                }
            }
        }
        else
        {
            // Form a Brillouin-zone weighted partial density for the selected band.
            for (int ik = 0; ik < nks; ++ik)
            {
                const int spin_index = kv.isk[ik];

                kspw_psi->fix_k(ik);
                const std::complex<double>* wfcr_up
                    = transform_wfc(&kspw_psi[0](ib, 0), ik, wfcr_up_smooth, wfcr_up_smooth_host, wfcr_up_dense_host);
                const std::complex<double>* wfcr_up_host_data = wfcr_up;
                const std::complex<double>* wfcr_down_host_data = nullptr;
                if (is_spinor)
                {
                    const std::complex<double>* wfcr_down
                        = transform_wfc(&kspw_psi[0](ib, npwx), ik, wfcr_down_smooth, wfcr_down_smooth_host, wfcr_down_dense_host);
                    wfcr_down_host_data = wfcr_down;
                }

                const double weight = static_cast<double>(kv.wk[ik] / ucell->omega);
                if (is_spinor)
                {
                    accumulate_spinor_density(wfcr_up_host_data, wfcr_down_host_data, weight, true);
                }
                else
                {
                    for (int ir = 0; ir < dense_nrxx; ++ir)
                    {
                        rho_band[spin_index][ir] += std::norm(wfcr_up_host_data[ir]) * weight;
                    }
                }
            }

#ifdef __MPI
            if (kv.para_k.kpar > 1)
            {
                // Each pool owns only part of the k-point sum; assemble it before symmetrization.
                for (int is = 0; is < nspin; ++is)
                {
                    pgrid.reduce_across_pools(rho_band[is].data());
                }
            }
#endif

            // Symmetry_rho operates on arrays of component pointers and uses reciprocal workspaces.
            Symmetry_rho srho;
            std::vector<double*> rho_save_pointers(nspin);
            std::vector<std::vector<std::complex<double>>> rhog(nspin, std::vector<std::complex<double>>(pw_rhod->npw));
            std::vector<std::complex<double>*> rhog_pointers(nspin);
            for (int is = 0; is < nspin; ++is)
            {
                rho_save_pointers[is] = rho_band[is].data();
                rhog_pointers[is] = rhog[is].data();
            }
            if (is_spinor)
            {
                // Charge and magnetization components obey different spinor symmetry transformations.
                srho.begin(0, rho_save_pointers.data(), rhog_pointers.data(), pw_rhod->npw, nullptr, pw_rhod, ucell->symm);
                srho.begin_soc(rho_save_pointers.data(), rhog_pointers.data(), pw_rhod, ucell->symm);
            }
            else
            {
                for (int is = 0; is < nspin; ++is)
                {
                    srho.begin(is, rho_save_pointers.data(), rhog_pointers.data(), pw_rhod->npw, nullptr, pw_rhod, ucell->symm);
                }
            }

            for (int is = 0; is < nspin; ++is)
            {
                std::stringstream ssc;
                ssc << global_out_dir << "pchgi" << ib + 1 << "s" << is + 1 << ".cube";
                ModuleIO::write_vdata_palgrid(pgrid, rho_band[is].data(), is, nspin, 0, ssc.str(), 0.0, ucell, 11, 0, false, false);
            }
        }
    }
}
} // namespace ModuleIO

#endif // GET_PCHG_PW_H
