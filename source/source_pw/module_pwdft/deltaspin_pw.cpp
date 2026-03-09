#include "source_pw/module_pwdft/deltaspin_pw.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"

namespace pw
{

bool run_deltaspin_lambda_loop(const int iter,
                               const double drho,
                               const Input_para& inp)
{
    /// Return false if DeltaSpin is not enabled
    if (!inp.sc_mag_switch)
    {
        return false;
    }

    /// Get the singleton instance of SpinConstrain
    spinconstrain::SpinConstrain<std::complex<double>>& sc
        = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();

    /// Case 1: Magnetic moments not yet converged and SCF is close to convergence.
    /// This is the first time we enter the lambda loop after SCF is nearly converged.
    if (!sc.mag_converged() && drho > 0 && drho < inp.sc_scf_thr)
    {
        /// Optimize lambda to get target magnetic moments
        sc.run_lambda_loop(iter);
        sc.set_mag_converged(true);
        return true;
    }
    /// Case 2: Magnetic moments already converged in previous iteration.
    /// Continue to refine lambda in subsequent SCF iterations.
    else if (sc.mag_converged())
    {
        sc.run_lambda_loop(iter);
        return true;
    }

    /// Default: run the normal solver
    return false;
}

}
