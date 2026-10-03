#include "td_current_pw.h"

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_hsolver/kernels/linear_op.h"

// Keep operator dependencies local to this implementation file.
#include "source_pw/module_pwdft/op_pw_vel.h"

#include <algorithm>
#include <fstream>
#include <iomanip>

namespace ModuleIO
{

template <typename FPTYPE, typename Device>
PWCurrentResult CurrentPW<FPTYPE, Device>::calculate(const UnitCell& ucell,
                                                     const ModulePW::PW_Basis_K* wfcpw,
                                                     psi::Psi<std::complex<FPTYPE>, Device>* psi,
                                                     const elecstate::ElecState* pelec,
                                                     const K_Vectors& kv,
                                                     pseudopot_cell_vnl* ppcell,
                                                     const int gauge,
                                                     const ModuleBase::Vector3<double>& vector_potential)
{
    ModuleBase::timer::start("ModuleIO", "calculate");

    using Complex = std::complex<FPTYPE>;
    using syncmem_complex_d2h_op = base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, Device>;
    using setmem_complex_op = base_device::memory::set_memory_op<Complex, Device>;

    double current_total[3] = {0.0, 0.0, 0.0};
    const int nks = wfcpw->nks;
    const int max_npw = wfcpw->npwk_max;
    const int nbands = psi->get_nbands();
    const int nkstot = kv.get_nkstot();
    const int* isk = kv.isk.data();

    const int n_npwx = nbands;

    // ==============================================================
    // Store Cartesian current components for each k point.
    // ==============================================================
    std::vector<double> current_k(3 * nkstot, 0.0);

    // Refresh the potential view used by the shared velocity operator.
    const FPTYPE* vtau = nullptr;
    int vtau_col = 0;
    int vtau_row = 0;
    if (gauge == 0 && XC_Functional::get_ked_flag())
    {
        if (pelec->pot == nullptr)
        {
            ModuleBase::WARNING_QUIT("calculate", "Missing potential for meta-GGA current.");
        }
        vtau = pelec->pot->template get_vofk_smooth_data<FPTYPE>();
        vtau_col = pelec->pot->get_vofk_smooth().nc;
        vtau_row = pelec->pot->get_vofk_smooth().nr;
        if (vtau_col != wfcpw->nrxx || (wfcpw->nrxx > 0 && (vtau == nullptr || vtau_row != kv.get_spin_mult())))
        {
            ModuleBase::WARNING_QUIT("calculate", "Invalid smooth-grid potential for meta-GGA current.");
        }
        // Preserve the enabled correction on ranks without real-space planes, which still join FFTs.
        vtau_row = kv.get_spin_mult();
    }
    if (!velocity_)
    {
        velocity_.reset(new hamilt::Velocity<FPTYPE, Device>(wfcpw, isk, ppcell, &ucell, true, vtau, vtau_col, vtau_row));
    }
    velocity_->set_state(isk, vtau, vtau_col, vtau_row);

    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int64_t velocity_size = std::max<int64_t>(1, 3LL * n_npwx * max_npw);
    const int64_t dot_size = std::max(1, n_npwx);
    if (vpsi_.NumElements() < velocity_size || vpsi_.data_type() != ct::DataTypeToEnum<Complex>::value || vpsi_.device_type() != device)
    {
        vpsi_ = ct::Tensor(ct::DataTypeToEnum<Complex>::value, device, {velocity_size});
    }
    if (dots_.NumElements() < dot_size || dots_.data_type() != ct::DataTypeToEnum<Complex>::value || dots_.device_type() != device)
    {
        dots_ = ct::Tensor(ct::DataTypeToEnum<Complex>::value, device, {dot_size});
    }
    Complex* d_vpsi = vpsi_.template data<Complex>();
    ct::Tensor& dot_buffer = dots_;
    band_current_.resize(n_npwx);
    std::vector<Complex>& band_current = band_current_;

    for (int ik = 0; ik < nks; ++ik)
    {
        psi->fix_k(ik);
        Complex* current_psi_ptr = psi->get_pointer();
        const int npw = wfcpw->npwk[ik];

        if (n_npwx * max_npw > 0)
        {
            setmem_complex_op()(d_vpsi, 0, 3 * n_npwx * max_npw);
        }

        velocity_->init(ik, vector_potential);
        velocity_->act(psi, n_npwx, current_psi_ptr, d_vpsi, false);

        for (int id = 0; id < 3; ++id)
        {
            hsolver::linear_op<Complex, Device>()
                .dot(max_npw, npw, n_npwx, current_psi_ptr, d_vpsi + id * n_npwx * max_npw, dot_buffer.data<Complex>());
            syncmem_complex_d2h_op()(band_current.data(), dot_buffer.data<Complex>(), n_npwx);
            for (int ib = 0; ib < nbands; ++ib)
            {
                const double contribution = -pelec->wg(ik, ib) * std::real(band_current[ib]);
                current_total[id] += contribution;
                current_k[kv.ik2iktot[ik] * 3 + id] += contribution;
            }
        }
    }

    // ==============================================================
    // Reduce all current components together across MPI ranks.
    // ==============================================================
    Parallel_Reduce::reduce_all(current_total, 3);
    Parallel_Reduce::reduce_all(current_k.data(), 3 * nkstot);

    PWCurrentResult result;
    for (int d = 0; d < 3; ++d)
    {
        result.total[d] = current_total[d] / ucell.omega;
    }
    for (double& value: current_k)
    {
        value /= ucell.omega;
    }
    result.per_k = std::move(current_k);
    ModuleBase::timer::end("ModuleIO", "calculate");
    return result;
}

void write_pw_current(const PWCurrentResult& current,
                      const int istep,
                      const int nk_per_spin,
                      const bool out_current_k,
                      const std::string& out_dir)
{
    // Always write the total current when this routine is enabled.
    std::string filename_tot = out_dir + "current_tot.txt";
    std::ofstream fout_tot;
    fout_tot.open(filename_tot, std::ios::app);
    fout_tot << std::setprecision(16) << std::scientific;
    fout_tot << istep + 1 << " " << current.total[0] << " " << current.total[1] << " " << current.total[2] << std::endl;
    fout_tot.close();

    // Optionally write the contribution from each k point.
    if (out_current_k)
    {
        for (int ik = 0; ik < static_cast<int>(current.per_k.size() / 3); ++ik)
        {
            // Use a one-based spin index in the output filename.
            int is = ik / nk_per_spin + 1;
            int k_idx = ik % nk_per_spin + 1;

            std::string filename_k = out_dir + "current_s" + std::to_string(is) + "k" + std::to_string(k_idx) + ".txt";
            std::ofstream fout_k;
            fout_k.open(filename_k, std::ios::app);
            fout_k << std::setprecision(16) << std::scientific;
            fout_k << istep + 1 << " " << current.per_k[ik * 3 + 0] << " " << current.per_k[ik * 3 + 1] << " " << current.per_k[ik * 3 + 2]
                   << std::endl;
            fout_k.close();
        }
    }
}

template <typename FPTYPE, typename Device>
void CurrentPW<FPTYPE, Device>::write(const int istep,
                                      const UnitCell& ucell,
                                      const ModulePW::PW_Basis_K* wfcpw,
                                      psi::Psi<std::complex<FPTYPE>, Device>* psi,
                                      const elecstate::ElecState* pelec,
                                      const K_Vectors& kv,
                                      pseudopot_cell_vnl* ppcell,
                                      const int gauge,
                                      const ModuleBase::Vector3<double>& A_right_ha,
                                      const bool out_current_k,
                                      const std::string& out_dir,
                                      const int world_rank)
{
    const PWCurrentResult result = calculate(ucell, wfcpw, psi, pelec, kv, ppcell, gauge, A_right_ha);
    if (world_rank == 0)
    {
        const int nk_per_spin = kv.get_nkstot() / kv.get_spin_mult();
        write_pw_current(result, istep, nk_per_spin, out_current_k, out_dir);
    }
}

template class CurrentPW<float, base_device::DEVICE_CPU>;
template class CurrentPW<double, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class CurrentPW<float, base_device::DEVICE_GPU>;
template class CurrentPW<double, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleIO
