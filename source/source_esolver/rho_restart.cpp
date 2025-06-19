#include "rho_restart.h"

void ModuleESolver::rho_restart(const Input_para& inp,
		const UnitCell& ucell, 
        const elecstate::ElecState& elec,
		const int nrxx, // originally written as pw_rhod->nrxx  
        const int iter, // SCF iteration index
        const double& scf_thr,
        const double& scf_ene_thr,
        double& drho, // not sure how this is changed in this function
        Charge_Mixing& chr_mix,
        Charge &chr,
        bool &conv_esolver,
        bool &oscillate_esolver)
{
    ModuleBase::TITLE("ModuleESolver", "rho_restart");

    // mixing will restart once if drho is smaller than inp.mixing_restart
    const double restart_thr = inp.mixing_restart;

    // ks_run means this is KSDFT
    if (PARAM.globalv.ks_run)
    {
        //--------------------------------------------------------
        // step1: determine mixing_restart_step
        //--------------------------------------------------------
        // charge mixing restarts at chgmix.mixing_restart steps
        if (drho <= restart_thr 
            && restart_thr > 0.0
            && chgmix.mixing_restart_step > iter)
        {
            chgmix.mixing_restart_step = iter + 1;
        }


        //--------------------------------------------------------
        // step2: determine density oscillation 
        //--------------------------------------------------------
        // if density oscillation is detected, SCF will stop
        if (inp.scf_os_stop)
        {
			oscillate_esolver = chgmix.if_scf_oscillate(iter, 
					drho, 
					inp.scf_os_ndim, 
					inp.scf_os_thr);
        }

        //--------------------------------------------------------
        // step3: determine convergence of SCF: conv_esolver 
        //--------------------------------------------------------
        // drho will be 0 at the chgmix.mixing_restart step, 
        // which is not ground state
		bool is_mixing_restart_step = (iter == chgmix.mixing_restart_step);
		bool is_restart_thr_positive = (restart_thr > 0.0);
		bool is_restart_condition_met = is_mixing_restart_step && is_restart_thr_positive;
		bool not_restart_step =!is_restart_condition_met;

        // SCF will continue if U is not converged for uramping calculation
        bool is_U_converged = true;

        // to avoid unnecessary dependence on dft+u, refactor is needed
#ifdef __LCAO
        if (inp.dft_plus_u)
        {
            is_U_converged = GlobalC::dftu.u_converged();
        }
#endif

        conv_esolver = (drho < scf_thr && not_restart_step && is_U_converged);

        //--------------------------------------------------------
        // step4: determine conv_esolver if energy threshold is 
        // used in SCF
        //--------------------------------------------------------
        if (scf_ene_thr > 0.0)
        {
            // calculate energy of output charge density
            this->update_pot(ucell, istep, iter, conv_esolver);

            // '2' means Kohn-Sham functional
            elec.cal_energies(2);

            // now, etot_old is the energy of input density, while etot is the energy of output density
            elec.f_en.etot_delta = elec.f_en.etot - elec.f_en.etot_old;

            // output etot_delta
            GlobalV::ofs_running << " DeltaE_womix = " 
                                 << elec.f_en.etot_delta * ModuleBase::Ry_to_eV << " eV"
                                 << std::endl;

            // only check when density is converged
            if (iter > 1 && conv_esolver == 1)
            {
                // update the convergence flag
                conv_esolver
                    = (std::abs(elec.f_en.etot_delta * ModuleBase::Ry_to_eV) < scf_ene_thr);
            }
        }

        //--------------------------------------------------------
        // If drho < hsolver_error in the first iter or 
        // drho < scf_thr, we do nothing and do not change rho.
        //--------------------------------------------------------
        if (drho < hsolver_error  
            || conv_esolver   // SCF has been converged
            || inp.calculation == "nscf") // nscf calculations, do not change rho
        {
            if (drho < hsolver_error)
            {
                GlobalV::ofs_warning << " drho < hsolver_error, keep "
                                        "charge density unchanged."
                                     << std::endl;
            }
        }
        else
        {
            // mixing will restart after chgmix.mixing_restart steps
            if (restart_thr > 0.0 
                && iter == chgmix.mixing_restart_step - 1
                && drho <= restart_thr)
            {
                // do not mix charge density
            }
            else
            {
                // mix charge density (rho) 
                chgmix.mix_rho(&chr);
            }

            // renormalize rho in R-space would induce error in G space
            if (inp.scf_thr_type == 2)
            {
                chr.renormalize_rho();
            }
        }
    }

#ifdef __MPI
    // bcast drho in BP_WORLD (Band parallel world) 
    MPI_Bcast(&drho, 1, MPI_DOUBLE, 0, BP_WORLD);

    // be careful! conv_esolver is bool, not double !! Maybe a bug 20250302 by mohan 
    MPI_Bcast(&conv_esolver, 1, MPI_DOUBLE, 0, BP_WORLD);
    MPI_Bcast(chr.rho[0], nrxx, MPI_DOUBLE, 0, BP_WORLD);
#endif

}
