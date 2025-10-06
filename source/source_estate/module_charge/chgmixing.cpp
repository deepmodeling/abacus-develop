#include "source_estate/module_charge/chgmixing.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"

void module_charge::chgmixing(const int iter, // scf iteration number
        Charge_Mixing* p_chgmix, // charge mixing class
		const Input_para& inp) // input parameters
{
    ModuleBase::TITLE("module_charge", "setup_chgmixing");

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
