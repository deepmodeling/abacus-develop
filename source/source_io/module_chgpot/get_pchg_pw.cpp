#include "source_io/module_chgpot/get_pchg_pw.h"

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_base/tool_quit.h"
#include "source_estate/module_charge/symm_rho.h"
#include "source_estate/uspp_density.h"
#include "source_io/module_output/cube_io.h"

#include <algorithm>
#include <sstream>
#include <type_traits>

namespace ModuleIO
{
// This nested class owns scratch storage for one begin() call and can access the output object's private data.
template <typename T, typename Device>
class Get_pchg_pw<T, Device>::Workspace
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

    std::vector<std::vector<double>> density;

    // USPP states are normalized with the overlap operator S. Their valence density
    // combines the soft |psi|^2 term with atom-centered augmentation functions Q_ij.
    std::unique_ptr<elecstate::UsppProjector<T, Device>> projector;
    std::vector<double> state_weight;
    // Packed projector-pair weights for one (band,k,spin) state and its accumulated k sum.
    // Off-diagonal entries contain twice the real part to combine conjugate contributions.
    std::vector<double> state_becsum;
    std::vector<double> becsum;
    std::vector<std::vector<std::complex<double>>> augmentation_g;
    std::vector<std::complex<double>*> augmentation_pointers;
    std::vector<double> augmentation_r;

    // Brace initialization constructs both spinor slots; unused device/grid buffers have zero size.
    explicit Workspace(const Get_pchg_pw& output, const UnitCell& ucell)
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
          density(output.nspin_, std::vector<double>(dense_nrxx))
    {
        for (int it = 0; it < ucell.ntype; ++it)
        {
            if (ucell.atoms[it].ncpp.tvanp)
            {
                // tvanp identifies a USPP species. One accumulator covers all USPP atoms,
                // so allocate it once when the first such species is found, also in mixed NC/USPP cells.
                projector.reset(new elecstate::UsppProjector<T, Device>(ucell, output.ppcell_, 1));
                state_weight.resize(1);
                // nspin-sized storage alone does not provide spinor USPP support; the shared
                // projector and augmentation contracts must also be extended when upstream enables nspin=4.
                becsum.resize(elecstate::uspp_becsum_size(ucell, output.ppcell_, output.nspin_));
                state_becsum.resize(becsum.size());
                augmentation_g.assign(output.nspin_, std::vector<std::complex<double>>(output.pw_rhod_.npw));
                augmentation_pointers.resize(output.nspin_);
                augmentation_r.resize(dense_nrxx);
                for (int is = 0; is < output.nspin_; ++is)
                {
                    augmentation_pointers[is] = augmentation_g[is].data();
                }
                break;
            }
        }
    }
};

template <typename T, typename Device>
Get_pchg_pw<T, Device>::Get_pchg_pw(const psi::Psi<T, Device>& psi,
                                    const ModulePW::PW_Basis_K& pw_wfc,
                                    const ModulePW::PW_Basis& pw_rho,
                                    const ModulePW::PW_Basis& pw_rhod,
                                    const pseudopot_cell_vnl& ppcell,
                                    const int nspin,
                                    const int global_nbands)
    : psi_(psi), pw_wfc_(pw_wfc), pw_rho_(pw_rho), pw_rhod_(pw_rhod), ppcell_(ppcell), nspin_(nspin), global_nbands_(global_nbands)
{
}

template <typename T, typename Device>
void Get_pchg_pw<T, Device>::begin(UnitCell* ucell,
                                   const Parallel_Grid& pgrid,
                                   const K_Vectors& kv,
                                   const std::vector<int>& out_pchg,
                                   const std::string& global_out_dir,
                                   const bool if_separate_k,
                                   const bool noncolin) const
{
    // Resolve global band ownership collectively before validating the selection.
    const BandParallelLayout layout(psi_.get_nbands(), global_nbands_);
    if (static_cast<int>(out_pchg.size()) > global_nbands_)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::get_pchg_pw",
                                 "The number of bands specified by `out_pchg` in the INPUT file exceeds `nbands`!");
    }
    const std::vector<int> band_mask = select_bands(out_pchg, "out_pchg");
    Workspace work(*this, *ucell);
    for (int band = 0; band < global_nbands_; ++band)
    {
        if (!band_mask[band])
        {
            continue;
        }
        for (int is = 0; is < nspin_; ++is)
        {
            std::fill(work.density[is].begin(), work.density[is].end(), 0.0);
        }
        // Each output band starts a fresh sum of augmentation contributions over k points.
        std::fill(work.becsum.begin(), work.becsum.end(), 0.0);
        if (if_separate_k)
        {
            write_separate(band, *ucell, pgrid, kv, global_out_dir, noncolin, layout, &work);
        }
        else
        {
            write_summed(band, ucell, pgrid, kv, global_out_dir, noncolin, layout, &work);
        }
    }
}

template <typename T, typename Device>
std::vector<int> Get_pchg_pw<T, Device>::select_bands(const std::vector<int>& selection, const std::string& parameter_name) const
{
    // begin() checks all selection lengths first; omitted bands remain unselected.
    std::vector<int> band_mask(global_nbands_, 0);
    for (int value: selection)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::get_pchg_pw",
                                     "The elements of `" + parameter_name + "` must be either 0 or 1. Invalid values found!");
        }
    }
    std::copy(selection.begin(), selection.end(), band_mask.begin());
    return band_mask;
}

template <typename T, typename Device>
void Get_pchg_pw<T, Device>::transform_band(const int global_band, const int ik, const BandParallelLayout& layout, Workspace* work) const
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
const std::complex<double>* Get_pchg_pw<T, Device>::transform_wfc(const T* coefficients,
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
void Get_pchg_pw<T, Device>::write_separate(const int band,
                                            const UnitCell& ucell,
                                            const Parallel_Grid& pgrid,
                                            const K_Vectors& kv,
                                            const std::string& out_dir,
                                            const bool noncolin,
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
        // Per-k states carry spin degeneracy, without a Brillouin-zone weight.
        const double spin_degeneracy = nspin_ == 1 ? 2.0 : 1.0;
        // Divide by the cell volume to convert the squared FFT amplitudes to a density.
        calc_density(spin_index, spin_degeneracy / ucell.omega, noncolin, false, work);
        // Each separate-k file needs its own augmentation, weighted by the same spin degeneracy as the soft term.
        std::fill(work->becsum.begin(), work->becsum.end(), 0.0);
        accumulate_uspp(band, ik, spin_index, spin_degeneracy, layout, work);
        add_augmentation(ucell, work);
        // Scalar/collinear output selects isk; spinors emit charge and all magnetization components.
        const int component_begin = work->is_spinor ? 0 : spin_index;
        const int component_end = work->is_spinor ? 4 : spin_index + 1;
        for (int component = component_begin; component < component_end; ++component)
        {
            write_cube(band, component, k_number, ucell, pgrid, out_dir, true, work->density[component]);
        }
    }
}

template <typename T, typename Device>
void Get_pchg_pw<T, Device>::write_summed(const int band,
                                          UnitCell* ucell,
                                          const Parallel_Grid& pgrid,
                                          const K_Vectors& kv,
                                          const std::string& out_dir,
                                          const bool noncolin,
                                          const BandParallelLayout& layout,
                                          Workspace* work) const
{
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        transform_band(band, ik, layout, work);
        // wk supplies the k-point weight (including spin degeneracy); omega normalizes the density.
        calc_density(kv.isk[ik], kv.wk[ik] / ucell->omega, noncolin, true, work);
        // Use the same k-point weight for the soft and augmentation contributions.
        accumulate_uspp(band, ik, kv.isk[ik], kv.wk[ik], layout, work);
    }
    // Form rho_soft + rho_aug in each pool, then sum over k pools.
    // Symmetry must act on this complete density so both contributions receive the same transformation.
    add_augmentation(*ucell, work);
    sum_pools(pgrid, kv, work);
    symmetrize(ucell, work);
    for (int is = 0; is < nspin_; ++is)
    {
        write_cube(band, is, 0, *ucell, pgrid, out_dir, false, work->density[is]);
    }
}

template <typename T, typename Device>
void Get_pchg_pw<T, Device>::calc_density(const int spin_index,
                                          const double weight,
                                          const bool noncolin,
                                          const bool accumulate,
                                          Workspace* work) const
{
    // The Bloch phase cancels in the squared modulus; weight includes the inverse cell volume.
    std::vector<std::vector<double>>& rho = work->density;
    const std::complex<double>* up = work->wfcr[0].data();
    if (!work->is_spinor)
    {
        for (int ir = 0; ir < work->dense_nrxx; ++ir)
        {
            const double value = std::norm(up[ir]) * weight;
            if (accumulate)
            {
                rho[spin_index][ir] += value;
            }
            else
            {
                rho[spin_index][ir] = value;
            }
        }
        return;
    }

    // Convert the spinor into (charge, m_x, m_y, m_z). Per-k output overwrites
    // each field; k-summed output accumulates the explicitly weighted fields.
    // Magnetization is u^dagger sigma u; noncolin controls whether m_x and m_y are retained.
    const std::complex<double>* down = work->wfcr[1].data();
    for (int ir = 0; ir < work->dense_nrxx; ++ir)
    {
        const double up_norm = std::norm(up[ir]);
        const double down_norm = std::norm(down[ir]);
        const double rho0 = (up_norm + down_norm) * weight;
        const double mx = 2.0 * (up[ir].real() * down[ir].real() + up[ir].imag() * down[ir].imag()) * weight;
        const double my = 2.0 * (up[ir].real() * down[ir].imag() - down[ir].real() * up[ir].imag()) * weight;
        const double mz = (up_norm - down_norm) * weight;
        if (accumulate)
        {
            rho[0][ir] += rho0;
            rho[1][ir] += noncolin ? mx : 0.0;
            rho[2][ir] += noncolin ? my : 0.0;
            rho[3][ir] += mz;
        }
        else
        {
            rho[0][ir] = rho0;
            rho[1][ir] = noncolin ? mx : 0.0;
            rho[2][ir] = noncolin ? my : 0.0;
            rho[3][ir] = mz;
        }
    }
}

template <typename T, typename Device>
void Get_pchg_pw<T, Device>::accumulate_uspp(const int band,
                                             const int ik,
                                             const int spin,
                                             const double weight,
                                             const BandParallelLayout& layout,
                                             Workspace* work) const
{
    if (!work->projector)
    {
        return;
    }
    std::fill(work->state_becsum.begin(), work->state_becsum.end(), 0.0);
    const int owner = layout.owner_group(band);
    if (layout.band_group() == owner)
    {
        psi_.fix_k(ik);
        work->state_weight[0] = weight;
        // This call currently describes one scalar block and a collinear channel index.
        // Spinor support requires component layout information and cross-spin projector products.
        // Overlaps <beta_i|psi> measure this state's amplitudes in the atomic augmentation channels.
        // The helper sums them over plane-wave ranks, then forms weighted pairs <psi|beta_i><beta_j|psi>.
        work->projector->accumulate(ik,
                                    &psi_(layout.local_index(band), 0),
                                    psi_.get_nbasis(),
                                    psi_.get_current_ngk(),
                                    spin,
                                    work->state_weight,
                                    &work->state_becsum);
    }
#ifdef __MPI
    // Only the owner holds this band's coefficients; replicate its projector products to all band groups.
    // A broadcast gives every group the same augmentation without multiplying its charge by the group count.
    Parallel_Common::bcast_data(work->state_becsum.data(), static_cast<int>(work->state_becsum.size()), BP_WORLD, owner);
#endif
    for (std::size_t i = 0; i < work->becsum.size(); ++i)
    {
        work->becsum[i] += work->state_becsum[i];
    }
}

template <typename T, typename Device>
void Get_pchg_pw<T, Device>::add_augmentation(const UnitCell& ucell, Workspace* work) const
{
    if (!work->projector)
    {
        return;
    }
    for (int is = 0; is < nspin_; ++is)
    {
        std::fill(work->augmentation_g[is].begin(), work->augmentation_g[is].end(), std::complex<double>(0, 0));
    }
    // Build rho_aug(G) = sum_{I,i<=j} Q^I_ij(G) B^I_ij from the packed projector weights in becsum.
    // Q(G) already includes 1/omega, so B_ij carries only wk or spin degeneracy as its state weight.
    elecstate::add_uspp_density(ucell, ppcell_, pw_rhod_, nspin_, work->becsum, work->augmentation_pointers.data());
    // This component-wise addition can accommodate spinors once the builder supplies
    // the physical (rho, m_x, m_y, m_z) augmentation fields; looping over nspin alone does not construct them.
    for (int is = 0; is < nspin_; ++is)
    {
        // Add the augmentation on the dense grid to complete the S-normalized state's valence density.
        pw_rhod_.recip2real(work->augmentation_g[is].data(), work->augmentation_r.data());
        for (int ir = 0; ir < work->dense_nrxx; ++ir)
        {
            work->density[is][ir] += work->augmentation_r[ir];
        }
    }
}

template <typename T, typename Device>
void Get_pchg_pw<T, Device>::sum_pools(const Parallel_Grid& pgrid, const K_Vectors& kv, Workspace* work) const
{
#ifdef __MPI
    if (kv.para_k.kpar > 1)
    {
        // Complete the k sum before applying symmetry, without summing band replicas.
        for (int is = 0; is < nspin_; ++is)
        {
            pgrid.reduce_across_pools(work->density[is].data());
        }
    }
#endif
}

template <typename T, typename Device>
void Get_pchg_pw<T, Device>::symmetrize(UnitCell* ucell, Workspace* work) const
{
    Symmetry_rho srho;
    std::vector<double*> rho_pointers(nspin_);
    std::vector<std::vector<std::complex<double>>> rhog(nspin_, std::vector<std::complex<double>>(pw_rhod_.npw));
    std::vector<std::complex<double>*> rhog_pointers(nspin_);
    // These non-owning pointer arrays adapt vector storage to the symmetry interface.
    // The vectors retain ownership and remain alive throughout the symmetry calls.
    for (int is = 0; is < nspin_; ++is)
    {
        rho_pointers[is] = work->density[is].data();
        rhog_pointers[is] = rhog[is].data();
    }
    if (work->is_spinor)
    {
        // Charge and magnetization obey different spinor symmetry transformations.
        srho.begin(0, rho_pointers.data(), rhog_pointers.data(), pw_rhod_.npw, nullptr, &pw_rhod_, ucell->symm);
        srho.begin_soc(rho_pointers.data(), rhog_pointers.data(), &pw_rhod_, ucell->symm);
    }
    else
    {
        for (int is = 0; is < nspin_; ++is)
        {
            srho.begin(is, rho_pointers.data(), rhog_pointers.data(), pw_rhod_.npw, nullptr, &pw_rhod_, ucell->symm);
        }
    }
}

template <typename T, typename Device>
void Get_pchg_pw<T, Device>::write_cube(const int band,
                                        const int component,
                                        const int k_number,
                                        const UnitCell& ucell,
                                        const Parallel_Grid& pgrid,
                                        const std::string& out_dir,
                                        const bool separate_k,
                                        const std::vector<double>& values) const
{
    std::stringstream filename;
    filename << out_dir << "pchgi" << band + 1 << "s" << component + 1;
    if (separate_k)
    {
        filename << "k" << k_number;
    }
    filename << ".cube";
    ModuleIO::write_vdata_palgrid(pgrid, values.data(), component, nspin_, 0, filename.str(), 0.0, &ucell, 11, 0, false, separate_k);
}

// Explicit instantiation emits both precisions for each supported device from this .cpp file.
template class Get_pchg_pw<std::complex<float>, base_device::DEVICE_CPU>;
template class Get_pchg_pw<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class Get_pchg_pw<std::complex<float>, base_device::DEVICE_GPU>;
template class Get_pchg_pw<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleIO
