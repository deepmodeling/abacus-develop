#include "source_lcao/module_lr/lr_zeq.h"

#include "source_lcao/module_lr/lr_dm.h"
#include "source_lcao/module_lr/window_response.h"
#include "source_lcao/module_lr/zeq_solver.h"

#include <complex>

namespace LR
{

template <typename T>
ct::Tensor solve_zeq(const GradientInputs<T>& in, const ct::Tensor& xz, const int ispin, const int nst,
                     const bool openshell, const std::string& grad_solver, const std::string& in_dir,
                     const std::vector<std::vector<ct::Tensor>>& window_density)
{
    ct::Tensor z = LR_Util::newTensor<T>({nst, in.nloc});
    std::weak_ptr<PotHxcLR> pot_weak = in.pot[ispin];
    Z_vector_equation(in, xz.template data<T>(), z.template data<T>(), nst, pot_weak, in.spin_types[ispin], in_dir,
                      openshell, grad_solver, window_density);
    return z;
}

namespace
{
/// Shared by both shells: the ground-state density matrix and the force terms.
template <typename T>
module_dm::DensityMatrix<T, double> gs_dm_for(const GradientInputs<T>& in, const GsForceCtx<T>& gs)
{
    return make_gs_dm<T>(gs.psi_all, gs.wg_all, gs.pv_all, in.pmat, in.ucell, in.gd, in.orb_cutoff, in.nspin,
                         in.kv.kvec_d, in.nk);
}
} // namespace

template <typename T>
std::vector<ModuleBase::matrix> force_xz(const GradientInputs<T>& in, const GsForceCtx<T>& gs, const int nactive,
                                         const int ispin, const ct::Tensor& xz, const std::vector<double>& omega,
                                         const int label_begin)
{
    const int nst = static_cast<int>(omega.size());
    const LR::WindowResponse<T> window = LR::build_window_response(in, xz, nst, nactive, ispin, gs.in_dir);
    const ct::Tensor z = solve_zeq(in, xz, ispin, nst, gs.openshell, gs.grad_solver, gs.in_dir, window.density);
    const auto dm_gs = gs_dm_for(in, gs);
    LR_Force<T> force_terms(in.ucell, in.kv.kvec_d, in.pmat, gs.rhod, gs.rho, gs.vloc, gs.sfac, in.gd, gs.tcb,
                            gs.settings, in.ofs
#ifdef __EXX
                            ,
                            in.exx_lri, in.hybrid_alpha
#endif
    );
    return evaluate_closed_shell_force(in, force_terms, dm_gs, xz, z, omega, label_begin, ispin, window);
}

template <typename T>
std::vector<ModuleBase::matrix> force_xz_os(const GradientInputs<T>& in, const GsForceCtx<T>& gs,
                                            const std::vector<int>& nactive, const ct::Tensor& xz,
                                            const std::vector<double>& omega, const int label_begin)
{
    const int nst = static_cast<int>(omega.size());
    const LR::WindowResponse<T> window = LR::build_open_shell_window_response(in, xz, nst, nactive);
    const ct::Tensor z = solve_zeq(in, xz, 0, nst, gs.openshell, gs.grad_solver, gs.in_dir, window.density);
    const auto dm_gs = gs_dm_for(in, gs);
    LR_Force<T> force_terms(in.ucell, in.kv.kvec_d, in.pmat, gs.rhod, gs.rho, gs.vloc, gs.sfac, in.gd, gs.tcb,
                            gs.settings, in.ofs
#ifdef __EXX
                            ,
                            in.exx_lri, in.hybrid_alpha
#endif
    );
    return evaluate_open_shell_force(in, force_terms, dm_gs, xz, z, omega, label_begin, window);
}

template ct::Tensor solve_zeq<double>(const GradientInputs<double>&, const ct::Tensor&, int, int, bool,
                                      const std::string&, const std::string&,
                                      const std::vector<std::vector<ct::Tensor>>&);
template ct::Tensor solve_zeq<std::complex<double>>(const GradientInputs<std::complex<double>>&, const ct::Tensor&,
                                                    int, int, bool, const std::string&, const std::string&,
                                                    const std::vector<std::vector<ct::Tensor>>&);
template std::vector<ModuleBase::matrix> force_xz<double>(const GradientInputs<double>&, const GsForceCtx<double>&,
                                                          int, int, const ct::Tensor&, const std::vector<double>&, int);
template std::vector<ModuleBase::matrix> force_xz<std::complex<double>>(const GradientInputs<std::complex<double>>&,
                                                                        const GsForceCtx<std::complex<double>>&, int,
                                                                        int, const ct::Tensor&,
                                                                        const std::vector<double>&, int);
template std::vector<ModuleBase::matrix> force_xz_os<double>(const GradientInputs<double>&, const GsForceCtx<double>&,
                                                             const std::vector<int>&, const ct::Tensor&,
                                                             const std::vector<double>&, int);
template std::vector<ModuleBase::matrix> force_xz_os<std::complex<double>>(
    const GradientInputs<std::complex<double>>&, const GsForceCtx<std::complex<double>>&, const std::vector<int>&,
    const ct::Tensor&, const std::vector<double>&, int);

} // namespace LR
