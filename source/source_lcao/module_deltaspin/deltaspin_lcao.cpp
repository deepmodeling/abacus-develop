#include "deltaspin_lcao.h"
#include "spin_constrain.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_estate/elecstate.h"

namespace ModuleESolver
{

template <typename TK>
void init_deltaspin_lcao(const UnitCell& ucell,
                          const Input_para& inp,
                          void* pv,
                          const K_Vectors& kv,
                          void* p_hamilt,
                          void* psi,
                          void* dm,
                          void* pelec)
{
    if (!inp.sc_mag_switch)
    {
        return;
    }
    
    spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
#ifdef __LCAO
    sc.init_sc(inp.sc_thr, inp.nsc, inp.nsc_min, inp.alpha_trial,
               inp.sccut, inp.sc_drop_thr, ucell,
               static_cast<Parallel_Orbitals*>(pv),
               inp.nspin, kv, p_hamilt, psi,
               static_cast<elecstate::DensityMatrix<TK, double>*>(dm),
               static_cast<elecstate::ElecState*>(pelec));
#else
    sc.init_sc(inp.sc_thr, inp.nsc, inp.nsc_min, inp.alpha_trial,
               inp.sccut, inp.sc_drop_thr, ucell,
               static_cast<Parallel_Orbitals*>(pv),
               inp.nspin, kv, p_hamilt, psi,
               static_cast<elecstate::ElecState*>(pelec));
#endif
}

template <typename TK>
void cal_mi_lcao_wrapper(const int iter, const Input_para& inp)
{
    if (!inp.sc_mag_switch)
    {
        return;
    }
    
#ifdef __LCAO
    spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
    sc.cal_mi_lcao(iter);
#endif
}

template <typename TK>
bool run_deltaspin_lambda_loop_lcao(const int iter,
                                     const double drho,
                                     const Input_para& inp)
{
    bool skip_solve = false;
    
    if (inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        
        if (!sc.mag_converged() && drho > 0 && drho < inp.sc_scf_thr)
        {
            /// optimize lambda to get target magnetic moments, but the lambda is not near target
            sc.run_lambda_loop(iter);
            sc.set_mag_converged(true);
            skip_solve = true;
        }
        else if (sc.mag_converged())
        {
            /// optimize lambda to get target magnetic moments, but the lambda is not near target
            sc.run_lambda_loop(iter);
            skip_solve = true;
        }
    }
    
    return skip_solve;
}

/// Template instantiation
template void init_deltaspin_lcao<double>(const UnitCell& ucell,
                                           const Input_para& inp,
                                           void* pv,
                                           const K_Vectors& kv,
                                           void* p_hamilt,
                                           void* psi,
                                           void* dm,
                                           void* pelec);
template void init_deltaspin_lcao<std::complex<double>>(const UnitCell& ucell,
                                                         const Input_para& inp,
                                                         void* pv,
                                                         const K_Vectors& kv,
                                                         void* p_hamilt,
                                                         void* psi,
                                                         void* dm,
                                                         void* pelec);

template void cal_mi_lcao_wrapper<double>(const int iter, const Input_para& inp);
template void cal_mi_lcao_wrapper<std::complex<double>>(const int iter, const Input_para& inp);

template bool run_deltaspin_lambda_loop_lcao<double>(const int iter,
                                                      const double drho,
                                                      const Input_para& inp);
template bool run_deltaspin_lambda_loop_lcao<std::complex<double>>(const int iter,
                                                                     const double drho,
                                                                     const Input_para& inp);

} // namespace ModuleESolver
