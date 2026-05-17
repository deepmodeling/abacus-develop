#include "source_pw/module_pwdft/deltaspin_pw.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_estate/module_charge/charge_mixing.h"
#include "source_io/module_parameter/parameter.h"

namespace pw
{

bool run_deltaspin_lambda_loop(const int iter,
                               const double drho,
                               const Input_para& inp)
{
    if (!inp.sc_mag_switch)
    {
        return false;
    }

    spinconstrain::SpinConstrain<std::complex<double>>& sc
        = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();

    if (inp.sc_lambda_strategy == "linear_scan")
    {
        sc.run_lambda_linear_scan(iter);
        return true;
    }

    if (inp.sc_scf_thr_mode == "immediate")
    {
        // "immediate" mode: activate lambda loop from second SCF iteration.
        // First iteration (iter=0) is skipped because initial wavefunctions
        // are not available to compute initial magnetic moments.
        if (iter >= 1)
        {
            sc.run_lambda_loop(iter - 1);
            if (!sc.mag_converged()) { sc.set_mag_converged(true); }
            return true;
        }
        return false;
    }

    // "threshold" mode: activate when drho < sc_scf_thr
    // drho > 0 excludes iterations where drho has not been computed.
    if (!sc.mag_converged() && drho > 0 && drho < inp.sc_scf_thr)
    {
        sc.run_lambda_loop(iter - 1);
        sc.set_mag_converged(true);
        return true;
    }
    else if (sc.mag_converged())
    {
        sc.run_lambda_loop(iter - 1);
        return true;
    }

    return false;
}

void check_deltaspin_oscillation(const int iter,
                                 const double drho,
                                 Charge_Mixing* p_chgmix,
                                 const Input_para& inp)
{
    /// Return if DeltaSpin is not enabled
    if (!inp.sc_mag_switch)
    {
        return;
    }

    /// Get the singleton instance of SpinConstrain
    spinconstrain::SpinConstrain<std::complex<double>>& sc
        = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();

    /// Check if higher magnetization precision is needed
    if (!sc.higher_mag_prec)
    {
        /// Detect SCF oscillation
        sc.higher_mag_prec = p_chgmix->if_scf_oscillate(iter, drho, inp.sc_os_ndim, inp.scf_os_thr);

        /// If oscillation detected, set mixing restart step for next iteration
        if (sc.higher_mag_prec)
        {
            p_chgmix->mixing_restart_step = iter + 1;
        }
    }
}

}
