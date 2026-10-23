#include "source_pw/module_rt/prop_cn_pw.h"

#include "source_base/kernels/math_kernel_op.h"

namespace pw
{
template <typename T, typename Device>
CNPropPW<T, Device>::CNPropPW(const ModulePW::PW_Basis_K& basis,
                              const hsolver::PWLinearOptions& cn_opts,
                              const hsolver::diag_comm_info& comm,
                              std::ostream& log)
    : solver_(basis, cn_opts, comm, log)
{
}
template <typename T, typename Device>
void CNPropPW<T, Device>::invalidate_basis()
{
    solver_.invalidate_basis();
}
template <typename T, typename Device>
void CNPropPW<T, Device>::evolve(PropArgs<T, Device>& args, psi::Psi<T, Device>* current)
{
    using Real = typename T::value_type;
    const Real* potential = args.potential.template get_veff_smooth_data<Real>();
    if (args.scf_iter > 1)
    {
        const ModuleBase::matrix& shape = args.potential.get_veff_smooth();
        const int size = shape.nr * shape.nc;
        hsolver::linear_buffer<Real, Device>(&veff_mid_, size);
        Real* buffer = veff_mid_.template data<Real>();
        const Real* previous = args.history.left_potential();
        ModuleBase::vector_add_vector_op<Real, Device, Real>()(size, buffer, previous, Real(0.5), potential, Real(0.5));
        potential = buffer;
    }
    ModuleBase::Vector3<double> vecpot_ha(0, 0, 0);
    if (args.field.gauge() == 1)
    {
        vecpot_ha = args.field.A_prop_ha();
    }
    args.hamilt.bind(potential, vecpot_ha);
    const auto& previous = args.history.previous();
    const double dt = args.field.interval_dt_ha();
    solver_.solve(args.hamilt, previous, current, dt, vecpot_ha, args.elec_step, args.scf_iter, args.log_detail);
}
template class CNPropPW<std::complex<float>, base_device::DEVICE_CPU>;
template class CNPropPW<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class CNPropPW<std::complex<float>, base_device::DEVICE_GPU>;
template class CNPropPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace pw
