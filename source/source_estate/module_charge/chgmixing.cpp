#include "source_estate/module_charge/chgmixing.h"
#include "source_estate/update_pot.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"

void module_charge::chgmixing_ks(const int iter, // scf iteration number
		UnitCell& ucell,
		elecstate::ElecState* pelec, 
        Charge &chr, // charge density
        Charge_Mixing* p_chgmix, // charge mixing class
        const int nrxx, // charge density
        double &drho, // charge density deviation
        bool &oscillate_esolver, // whether the esolver has oscillation of charge density
        bool &conv_esolver,
        const double &hsolver_error,
        const double &scf_thr,
        const double &scf_ene_thr,
		const Input_para& inp) // input parameters
{

    if (PARAM.globalv.ks_run)
    {
        // mixing will restart at p_chgmix->mixing_restart steps
        if (drho <= inp.mixing_restart && inp.mixing_restart > 0.0
            && p_chgmix->mixing_restart_step > iter)
        {
            p_chgmix->mixing_restart_step = iter + 1;
        }

        if (inp.scf_os_stop) // if oscillation is detected, SCF will stop
        {
            oscillate_esolver = p_chgmix->if_scf_oscillate(iter, drho, 
               inp.scf_os_ndim, inp.scf_os_thr);
        }

        // drho will be 0 at p_chgmix->mixing_restart step, which is
        // not ground state
        bool not_restart_step = !(iter == p_chgmix->mixing_restart_step && inp.mixing_restart > 0.0);
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

        // add energy threshold for SCF convergence
        if (scf_ene_thr > 0.0)
        {
            // calculate energy of output charge density
            elecstate::update_pot(ucell, pelec, chr, conv_esolver);
            pelec->cal_energies(2); // 2 means Kohn-Sham functional
            // now, etot_old is the energy of input density, while etot is the energy of output density
            pelec->f_en.etot_delta = pelec->f_en.etot - pelec->f_en.etot_old;
            // output etot_delta
            GlobalV::ofs_running << " DeltaE_womix = " << pelec->f_en.etot_delta * ModuleBase::Ry_to_eV << " eV"
                                 << std::endl;
            if (iter > 1 && conv_esolver == 1) // only check when density is converged
            {
                // update the convergence flag
                conv_esolver
                    = (std::abs(pelec->f_en.etot_delta * ModuleBase::Ry_to_eV) < scf_ene_thr);
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
            // mixing will restart after p_chgmix->mixing_restart
            // steps
            if (inp.mixing_restart > 0 && iter == p_chgmix->mixing_restart_step - 1
                && drho <= inp.mixing_restart)
            {
                // do not mix charge density
            }
            else
            {
                p_chgmix->mix_rho(&chr); // update chr->rho by mixing
            }
            if (inp.scf_thr_type == 2)
            {
                chr.renormalize_rho(); // renormalize rho in R-space would
                                                  // induce a error in K-space
            }
            //----------charge mixing done-----------
        }
	}

#ifdef __MPI
    MPI_Bcast(&drho, 1, MPI_DOUBLE, 0, BP_WORLD);

    // change MPI_DOUBLE to MPI_C_BOOL, mohan 2025-04-13
    MPI_Bcast(&conv_esolver, 1, MPI_C_BOOL, 0, BP_WORLD);

    assert(nrxx>=0); // mohan add 2025-10-18
    MPI_Bcast(chr.rho[0], nrxx, MPI_DOUBLE, 0, BP_WORLD);
#endif

    // mohan move the following code here, 2025-10-18
    // SCF restart information
    if (PARAM.inp.mixing_restart > 0
        && iter == p_chgmix->mixing_restart_step - 1
        && iter != PARAM.inp.scf_nmax)
    {
        p_chgmix->mixing_restart_last = iter;
        std::cout << " SCF restart after this step!" << std::endl;
    }

    return;
}



void module_charge::chgmixing_ks_pw(const int iter, // scf iteration number
        Charge_Mixing* p_chgmix, // charge mixing class
		const Input_para& inp) // input parameters
{
    ModuleBase::TITLE("module_charge", "chgmixing_ks_pw");

    if (iter == 1)
    {
        p_chgmix->init_mixing();
        p_chgmix->mixing_restart_step = inp.scf_nmax + 1;
    }

    // For mixing restart
    if (iter == p_chgmix->mixing_restart_step && inp.mixing_restart > 0.0)
    {
        p_chgmix->init_mixing();
        p_chgmix->mixing_restart_count++;

        if (inp.dft_plus_u)
        {
            auto* dftu = ModuleDFTU::DFTU::get_instance();
            if (dftu->uramping > 0.01 && !dftu->u_converged())
            {
                p_chgmix->mixing_restart_step = inp.scf_nmax + 1;
            }
            if (dftu->uramping > 0.01)
            {
                bool do_uramping = true;
                if (inp.sc_mag_switch)
                {
                    spinconstrain::SpinConstrain<std::complex<double>>& sc
                        = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
                    if (!sc.mag_converged()) // skip uramping if mag not converged
                    {
						do_uramping = false;
					}
				}
				if (do_uramping)
				{
					dftu->uramping_update(); // update U by uramping if uramping > 0.01
					std::cout << " U-Ramping! Current U = ";
					for (int i = 0; i < dftu->U0.size(); i++)
					{
						std::cout << dftu->U[i] * ModuleBase::Ry_to_eV << " ";
					}
					std::cout << " eV " << std::endl;
				}
			}
		}
	}

    return;
}

void module_charge::chgmixing_ks_lcao(const int iter, // scf iteration number
        Charge_Mixing* p_chgmix, // charge mixing class
        const int nnr, // dimension of density matrix
		const Input_para& inp) // input parameters
{
    ModuleBase::TITLE("module_charge", "chgmixing_ks_lcao");

    if (iter == 1)
    {
        p_chgmix->mix_reset(); // init mixing
        p_chgmix->mixing_restart_step = inp.scf_nmax + 1;
        p_chgmix->mixing_restart_count = 0;
        // this output will be removed once the feeature is stable
        if (GlobalC::dftu.uramping > 0.01)
        {
            std::cout << " U-Ramping! Current U = ";
            for (int i = 0; i < GlobalC::dftu.U0.size(); i++)
            {
                std::cout << GlobalC::dftu.U[i] * ModuleBase::Ry_to_eV << " ";
            }
            std::cout << " eV " << std::endl;
        }
    }

    // for mixing restart
    if (iter == p_chgmix->mixing_restart_step && inp.mixing_restart > 0.0)
    {
        p_chgmix->init_mixing();
        p_chgmix->mixing_restart_count++;
        if (inp.dft_plus_u)
        {
            GlobalC::dftu.uramping_update(); // update U by uramping if uramping > 0.01
            if (GlobalC::dftu.uramping > 0.01)
            {
                std::cout << " U-Ramping! Current U = ";
                for (int i = 0; i < GlobalC::dftu.U0.size(); i++)
                {
                    std::cout << GlobalC::dftu.U[i] * ModuleBase::Ry_to_eV << " ";
                }
                std::cout << " eV " << std::endl;
            }
            if (GlobalC::dftu.uramping > 0.01 && !GlobalC::dftu.u_converged())
            {
                p_chgmix->mixing_restart_step = inp.scf_nmax + 1;
            }
        }
        if (inp.mixing_dmr) // for mixing_dmr
        {
            // allocate memory for dmr_mdata
            p_chgmix->allocate_mixing_dmr(nnr);
        }
    }
}
