#include "source_pw/module_rt/td_vel_pw.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_hsolver/module_linalg/kernels/linear_op.h"

#include <algorithm>
#include <cstdint>
#include <limits>

namespace pw
{
template <typename Real, typename Device>
void TDVelPW<Real, Device>::initialize(const UnitCell& cell, const ModulePW::PW_Basis_K& basis, const K_Vectors& kv, pseudopot_cell_vnl& pp)
{
    if (cell_ != nullptr)
    {
        ModuleBase::WARNING_QUIT("TDVelPW", "Initialize once before sampling.");
    }
    cell_ = &cell;
    basis_ = &basis;
    kv_ = &kv;
    pp_ = &pp;
}

template <typename Real, typename Device>
TDVelResultPW TDVelPW<Real, Device>::calculate(const psi::Psi<Complex, Device>& psi,
                                               const elecstate::ElecState& state,
                                               const ModuleBase::Vector3<double>& a,
                                               const CurrentRequest& request)
{
    ModuleBase::timer::start("TDVelPW", "calculate");
    if (cell_ == nullptr)
    {
        ModuleBase::WARNING_QUIT("TDVelPW", "Initialize before sampling observables.");
    }
    TDVelResultPW result;
    if (!request.charge_total && !request.charge_per_k && !request.spin_current)
    {
        ModuleBase::timer::end("TDVelPW", "calculate");
        return result;
    }
    const int npol = psi.get_npol();
    const int npwx = basis_->npwk_max;
    const int nbands = psi.get_nbands();
    const int ncols = nbands * npol;
    const int old_k = psi.get_current_k();
    if (psi.get_nbasis() != npwx * npol || npol != cell_->get_npol() || (request.spin_current && npol != 2))
    {
        ModuleBase::WARNING_QUIT("TDVelPW", "Incompatible wavefunction layout or spin-current request.");
    }

    if (!vel_)
    {
        vel_.reset(new hamilt::Velocity<Real, Device>(basis_, kv_->isk.data(), pp_, cell_, true, nullptr, 0, 0));
    }

    using CtDevice = typename ct::PsiToContainer<Device>::type;
    using CopyToHost = base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, Device>;
    using Zero = base_device::memory::set_memory_op<Complex, Device>;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const std::int64_t comp_size = static_cast<std::int64_t>(ncols) * npwx;
    const std::int64_t vel_size = 3 * comp_size;
    const std::int64_t capacity = std::max<std::int64_t>(1, vel_size);
    const std::int64_t dot_size = std::max(1, ncols);
    // Default tensors contain one CPU float element; capacity alone does not establish a usable buffer.
    if (vpsi_.NumElements() < capacity || vpsi_.data_type() != ct::DataTypeToEnum<Complex>::value || vpsi_.device_type() != device)
    {
        vpsi_ = ct::Tensor(ct::DataTypeToEnum<Complex>::value, device, {capacity});
    }
    if (dots_.NumElements() < dot_size || dots_.data_type() != ct::DataTypeToEnum<Complex>::value || dots_.device_type() != device)
    {
        dots_ = ct::Tensor(ct::DataTypeToEnum<Complex>::value, device, {dot_size});
    }
    Complex* vpsi = vpsi_.template data<Complex>();
    Complex* dots = dots_.template data<Complex>();
    diag_.resize(ncols);
    if (request.spin_current)
    {
        cross_.resize(nbands);
    }
    if (request.charge_per_k)
    {
        const std::int64_t size = 3 * static_cast<std::int64_t>(kv_->get_nkstot());
        result.per_k.resize(size, 0.0);
    }

    for (int ik = 0; ik < basis_->nks; ++ik)
    {
        psi.fix_k(ik);
        const Complex* data = psi.get_pointer();
        const int npw = basis_->npwk[ik];
        if (vel_size > 0)
        {
            Zero()(vpsi, 0, vel_size);
        }
        vel_->init(ik, a);
        vel_->act(&psi, ncols, data, vpsi, false);
        if (npw == 0)
        {
            // Projection reductions above are collective; only local inner products can be skipped.
            continue;
        }
        for (int id = 0; id < 3; ++id)
        {
            const Complex* component = vpsi + id * comp_size;
            // Each spin component is one npwx-strided column; reuse its diagonal product for both currents.
            hsolver::linear_op<Complex, Device>().dot(npwx, npw, ncols, data, component, dots);
            CopyToHost()(diag_.data(), dots, ncols);
            if (request.charge_total || request.charge_per_k)
            {
                for (int ib = 0; ib < nbands; ++ib)
                {
                    double value = 0.0;
                    for (int spin = 0; spin < npol; ++spin)
                    {
                        value += std::real(diag_[ib * npol + spin]);
                    }
                    const double charge = -state.wg(ik, ib) * value;
                    if (request.charge_total)
                    {
                        result.total[id] += charge;
                    }
                    if (request.charge_per_k)
                    {
                        const std::int64_t index = 3 * static_cast<std::int64_t>(kv_->ik2iktot[ik]) + id;
                        result.per_k[index] += charge;
                    }
                }
            }
            if (request.spin_current)
            {
                const int ld_psi = 2 * npwx;
                // Sigma_z uses the diagonal products already needed by charge current; S = sigma/2.
                for (int spin = 0; spin < 2; ++spin)
                {
                    const double sign = spin == 0 ? 1.0 : -1.0;
                    for (int ib = 0; ib < nbands; ++ib)
                    {
                        const double weight = 0.5 * state.wg(ik, ib);
                        result.spin[3 * id + 2] += weight * sign * std::real(diag_[2 * ib + spin]);
                    }
                }
                // Sigma_x adds the cross products, while sigma_y takes their imaginary parts with opposite signs.
                for (int left = 0; left < 2; ++left)
                {
                    const int right = 1 - left;
                    const Complex* bra = data + left * npwx;
                    const Complex* ket = component + right * npwx;
                    hsolver::linear_op<Complex, Device>().dot(ld_psi, npw, nbands, bra, ket, dots);
                    CopyToHost()(cross_.data(), dots, nbands);
                    const double sign = left == 0 ? 1.0 : -1.0;
                    for (int ib = 0; ib < nbands; ++ib)
                    {
                        const double weight = 0.5 * state.wg(ik, ib);
                        result.spin[3 * id] += weight * std::real(cross_[ib]);
                        result.spin[3 * id + 1] += weight * sign * std::imag(cross_[ib]);
                    }
                }
            }
        }
    }
    // Sampling must not leave a different k-point selected for the propagation driver.
    psi.fix_k(old_k);
    if (request.charge_total)
    {
        Parallel_Reduce::reduce_all(result.total.data(), 3);
    }
    const std::int64_t max_chunk = std::numeric_limits<int>::max();
    const std::int64_t k_size = result.per_k.size();
    for (std::int64_t offset = 0; offset < k_size;)
    {
        const std::int64_t remaining = k_size - offset;
        const int count = static_cast<int>(std::min(max_chunk, remaining));
        double* chunk = result.per_k.data() + offset;
        Parallel_Reduce::reduce_all(chunk, count);
        offset += count;
    }
    if (request.spin_current)
    {
        Parallel_Reduce::reduce_all(result.spin.data(), 9);
    }
    for (double& value: result.total)
    {
        value /= cell_->omega;
    }
    for (double& value: result.per_k)
    {
        value /= cell_->omega;
    }
    for (double& value: result.spin)
    {
        value /= cell_->omega;
    }
    ModuleBase::timer::end("TDVelPW", "calculate");
    return result;
}

template class TDVelPW<float, base_device::DEVICE_CPU>;
template class TDVelPW<double, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class TDVelPW<float, base_device::DEVICE_GPU>;
template class TDVelPW<double, base_device::DEVICE_GPU>;
#endif
} // namespace pw
