#include "deltaspin_lcao.h"
#include "spin_constrain.h"

namespace ModuleESolver
{

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
            sc.run_lambda_loop(iter - 1);
            sc.set_mag_converged(true);
            skip_solve = true;
        }
        else if (sc.mag_converged())
        {
            /// optimize lambda to get target magnetic moments, but the lambda is not near target
            sc.run_lambda_loop(iter - 1);
            skip_solve = true;
        }
    }
    
    return skip_solve;
}

/// Template instantiation
template bool run_deltaspin_lambda_loop_lcao<double>(const int iter,
                                                      const double drho,
                                                      const Input_para& inp);
template bool run_deltaspin_lambda_loop_lcao<std::complex<double>>(const int iter,
                                                                     const double drho,
                                                                     const Input_para& inp);

} // namespace ModuleESolver
