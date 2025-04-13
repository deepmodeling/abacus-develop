

void ModuleESolver::rho_mix(const Input_para& inp,
		const UnitCell& ucell, 
        Charge_Mixing& chr_mix)
{
    
    // ks_run means this is KSDFT, otherwise it is OFDFT
    if (PARAM.globalv.ks_run)
    {
        // mixing will restart at this->p_chgmix->mixing_restart steps
        if (drho <= inp.mixing_restart && inp.mixing_restart > 0.0
            && chgmix.mixing_restart_step > iter)
        {
            this->p_chgmix->mixing_restart_step = iter + 1;
        }

        if (inp.scf_os_stop) // if oscillation is detected, SCF will stop
        {
            this->oscillate_esolver
                = this->p_chgmix->if_scf_oscillate(iter, drho, inp.scf_os_ndim, inp.scf_os_thr);
        }

        // drho will be 0 at this->p_chgmix->mixing_restart step, which is
        // not ground state
        bool not_restart_step = !(iter == this->p_chgmix->mixing_restart_step && inp.mixing_restart > 0.0);
        // SCF will continue if U is not converged for uramping calculation
        bool is_U_converged = true;
        // to avoid unnecessary dependence on dft+u, refactor is needed
#ifdef __LCAO
        if (inp.dft_plus_u)
        {
            is_U_converged = GlobalC::dftu.u_converged();
        }
#endif

        conv_esolver = (drho < this->scf_thr && not_restart_step && is_U_converged);

        // add energy threshold for SCF convergence
        if (this->scf_ene_thr > 0.0)
        {
            // calculate energy of output charge density
            this->update_pot(ucell, istep, iter, conv_esolver);
            this->pelec->cal_energies(2); // 2 means Kohn-Sham functional
            // now, etot_old is the energy of input density, while etot is the energy of output density
            this->pelec->f_en.etot_delta = this->pelec->f_en.etot - this->pelec->f_en.etot_old;
            // output etot_delta
            GlobalV::ofs_running << " DeltaE_womix = " << this->pelec->f_en.etot_delta * ModuleBase::Ry_to_eV << " eV"
                                 << std::endl;
            if (iter > 1 && conv_esolver == 1) // only check when density is converged
            {
                // update the convergence flag
                conv_esolver
                    = (std::abs(this->pelec->f_en.etot_delta * ModuleBase::Ry_to_eV) < this->scf_ene_thr);
            }
        }

        // If drho < hsolver_error in the first iter or drho < scf_thr, we
        // do not change rho.
        if (drho < hsolver_error || conv_esolver || inp.calculation == "nscf")
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
            //----------charge mixing---------------
            // mixing will restart after this->p_chgmix->mixing_restart
            // steps
            if (inp.mixing_restart > 0 && iter == this->p_chgmix->mixing_restart_step - 1
                && drho <= inp.mixing_restart)
            {
                // do not mix charge density
            }
            else
            {
                p_chgmix->mix_rho(&this->chr); // update chr->rho by mixing
            }
            if (inp.scf_thr_type == 2)
            {
                this->chr.renormalize_rho(); // renormalize rho in R-space would
                                             // induce a error in K-space
            }
            //----------charge mixing done-----------
        }
    }

    // BCAST should not be here, mohan note 2025-04-13 
#ifdef __MPI
    MPI_Bcast(&drho, 1, MPI_DOUBLE, 0, BP_WORLD);

    // be careful! conv_esolver is bool, not double !! Maybe a bug 20250302 by mohan 
    MPI_Bcast(&conv_esolver, 1, MPI_DOUBLE, 0, BP_WORLD);
    MPI_Bcast(this->chr.rho[0], this->pw_rhod->nrxx, MPI_DOUBLE, 0, BP_WORLD);
#endif

}
