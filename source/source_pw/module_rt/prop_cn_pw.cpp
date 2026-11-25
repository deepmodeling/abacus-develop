#include "source_pw/module_rt/prop_cn_pw.h"

#include "source_base/kernels/math_kernel_op.h"

namespace pw
{
template <typename T, typename Device>
CNPropPW<T, Device>::CNPropPW(const ModulePW::PW_Basis_K& basis,
                              const hsolver::CNControls& cn_ctrl,
                              int npol,
                              const std::vector<int>& global_k,
                              const hsolver::diag_comm_info& comm,
                              std::ostream& log)
    : solver_(basis, cn_ctrl, npol, global_k, comm, log)
{
}

template <typename T, typename Device>
void CNPropPW<T, Device>::invalidate_basis()
{
    solver_.invalidate_basis();
}

template <typename T, typename Device>
void CNPropPW<T, Device>::evolve(TDHamiltPW<T, Device>& hamilt,
                                 const TDHistoryPW<T, Device>& history,
                                 elecstate::Potential& potential,
                                 const elecstate::TDFieldManager& field,
                                 psi::Psi<T, Device>* current,
                                 int elec_step,
                                 int scf_iter,
                                 bool log_detail)
{
    using Real = typename T::value_type;
    const Real* veff = potential.template get_veff_smooth_data<Real>();
    if (scf_iter > 1)
    {
        const ModuleBase::matrix& shape = potential.get_veff_smooth();
        const int size = shape.nr * shape.nc;
        hsolver::linear_buffer<Real, Device>(&veff_mid_, size);
        Real* buffer = veff_mid_.template data<Real>();
        const Real* previous = history.left_potential();
        ModuleBase::vector_add_vector_op<Real, Device, Real>()(size, buffer, previous, Real(0.5), veff, Real(0.5));
        veff = buffer;
    }

    ModuleBase::Vector3<double> vecpot_ha(0, 0, 0);
    if (field.gauge() == 1)
    {
        vecpot_ha = field.A_prop_ha();
    }

    hamilt.bind(veff, vecpot_ha);
    const psi::Psi<T, Device>& previous = history.previous();
    const double dt = field.interval_dt_ha();
    solver_.solve(hamilt, previous, current, dt, vecpot_ha, elec_step, scf_iter, log_detail);
}
template class CNPropPW<std::complex<float>, base_device::DEVICE_CPU>;
template class CNPropPW<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class CNPropPW<std::complex<float>, base_device::DEVICE_GPU>;
template class CNPropPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace pw
