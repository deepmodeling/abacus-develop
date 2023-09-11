#include "Exx_LRI_interface.h"
#include "module_ri/exx_abfs-jle.h"
#include "module_ri/exx_opt_orb.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"

template<typename Tdata>
void Exx_LRI_Interface<Tdata>::write_Hexxs(const std::string &file_name) const
{
	ModuleBase::TITLE("Exx_LRI","write_Hexxs");
	ModuleBase::timer::tick("Exx_LRI", "write_Hexxs");
	std::ofstream ofs(file_name, std::ofstream::binary);
	cereal::BinaryOutputArchive oar(ofs);
	oar(exx_lri->Hexxs);
	ModuleBase::timer::tick("Exx_LRI", "write_Hexxs");
}

template<typename Tdata>
void Exx_LRI_Interface<Tdata>::read_Hexxs(const std::string &file_name)
{
	ModuleBase::TITLE("Exx_LRI","read_Hexxs");
	ModuleBase::timer::tick("Exx_LRI", "read_Hexxs");
	std::ifstream ifs(file_name, std::ofstream::binary);
	cereal::BinaryInputArchive iar(ifs);
	iar(exx_lri->Hexxs);
	ModuleBase::timer::tick("Exx_LRI", "read_Hexxs");
}
template<typename Tdata>
void Exx_LRI_Interface<Tdata>::exx_beforescf(const K_Vectors& kv, const Charge_Mixing& chgmix)
{
#ifdef __MPI
		if ( GlobalC::exx_info.info_global.cal_exx )
		{
            if (GlobalC::ucell.atoms[0].ncpp.xc_func == "HF" || GlobalC::ucell.atoms[0].ncpp.xc_func == "PBE0" || GlobalC::ucell.atoms[0].ncpp.xc_func == "HSE")
            {
                XC_Functional::set_xc_type("pbe");
            }
            else if (GlobalC::ucell.atoms[0].ncpp.xc_func == "SCAN0")
            {
                XC_Functional::set_xc_type("scan");
            }

			exx_lri->cal_exx_ions();
		}

		if (Exx_Abfs::Jle::generate_matrix)
		{
			//program should be stopped after this judgement
			Exx_Opt_Orb exx_opt_orb;
			exx_opt_orb.generate_matrix(kv);
			ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");
			return;
		}
		
		// set initial parameter for mix_DMk_2D
		if(GlobalC::exx_info.info_global.cal_exx)
		{
			exx_lri->mix_DMk_2D.set_nks(kv.nks, GlobalV::GAMMA_ONLY_LOCAL);
			if(GlobalC::exx_info.info_global.separate_loop)
			{
				if(GlobalC::exx_info.info_global.mixing_beta_for_loop1==1.0)
					exx_lri->mix_DMk_2D.set_mixing_mode(Mixing_Mode::No);
				else
					exx_lri->mix_DMk_2D.set_mixing_mode(Mixing_Mode::Plain)
					                .set_mixing_beta(GlobalC::exx_info.info_global.mixing_beta_for_loop1);
			}
			else
			{
				if(chgmix.get_mixing_mode() == "plain")
					exx_lri->mix_DMk_2D.set_mixing_mode(Mixing_Mode::Plain);
				else if(chgmix.get_mixing_mode() == "pulay")
					exx_lri->mix_DMk_2D.set_mixing_mode(Mixing_Mode::Pulay);
				else
					throw std::invalid_argument(
						"mixing_mode = " + chgmix.get_mixing_mode() + ", mix_DMk_2D unsupported.\n"
						+ std::string(__FILE__) + " line " + std::to_string(__LINE__));
            }
        }
        // for exx two_level scf
        exx_lri->two_level_step = 0;
#endif // __MPI
}

template<typename Tdata>
void Exx_LRI_Interface<Tdata>::exx_eachiterinit(const Local_Orbital_Charge& loc, const Charge_Mixing& chgmix, const int& iter)
{
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        if (!GlobalC::exx_info.info_global.separate_loop && exx_lri->two_level_step)
        {
			exx_lri->mix_DMk_2D.set_mixing_beta(chgmix.get_mixing_beta());
			if(chgmix.get_mixing_mode() == "pulay")
				exx_lri->mix_DMk_2D.set_coef_pulay(iter, chgmix);
			const bool flag_restart = (iter==1) ? true : false;
			if(GlobalV::GAMMA_ONLY_LOCAL)
				exx_lri->mix_DMk_2D.mix(loc.dm_gamma, flag_restart);
			else
				exx_lri->mix_DMk_2D.mix(loc.dm_k, flag_restart);

            exx_lri->cal_exx_elec(*loc.LOWF->ParaV);
        }
    }
}

template<typename Tdata>
void Exx_LRI_Interface<Tdata>::exx_hamilt2density(elecstate::ElecState& elec, const Parallel_Orbitals& pv)
{
    // Peize Lin add 2020.04.04
    if (XC_Functional::get_func_type() == 4 || XC_Functional::get_func_type() == 5)
    {
        // add exx
        // Peize Lin add 2016-12-03
        elec.set_exx(this->get_Eexx());

        if (GlobalC::restart.info_load.load_H && GlobalC::restart.info_load.load_H_finish
            && !GlobalC::restart.info_load.restart_exx)
        {
            XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);

            exx_lri->cal_exx_elec(pv);
            GlobalC::restart.info_load.restart_exx = true;
        }
    }
    else
    {
        elec.f_en.exx = 0.;
    }
}

template<typename Tdata>
bool Exx_LRI_Interface<Tdata>::exx_after_converge(
    hamilt::Hamilt<double>& hamilt,
    LCAO_Matrix& lm,
    const Local_Orbital_Charge& loc,
    const K_Vectors& kv,
    int& iter)
{
    // Add EXX operator
    auto add_exx_operator = [&]() {
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            hamilt::Operator<double>* exx
                = new hamilt::OperatorEXX<hamilt::OperatorLCAO<double>>(&lm,
                                                                        nullptr, // no explicit call yet
                                                                        &(lm.Hloc),
                                                                        kv);
            hamilt.opsd->add(exx);
        }
        else
        {
            hamilt::Operator<std::complex<double>>* exx
                = new hamilt::OperatorEXX<hamilt::OperatorLCAO<std::complex<double>>>(&lm,
                                                                                      nullptr, // no explicit call yet
                                                                                      &(lm.Hloc2),
                                                                                      kv);
            hamilt.ops->add(exx);
        }
    };
    
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        // no separate_loop case
        if (!GlobalC::exx_info.info_global.separate_loop)
        {
            GlobalC::exx_info.info_global.hybrid_step = 1;

            // in no_separate_loop case, scf loop only did twice
            // in first scf loop, exx updated once in beginning,
            // in second scf loop, exx updated every iter

            if (exx_lri->two_level_step)
            {
                return true;
            }
            else
            {
                // update exx and redo scf
                XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);
                iter = 0;
                std::cout << " Entering 2nd SCF, where EXX is updated" << std::endl;
                exx_lri->two_level_step++;

                add_exx_operator();

                return false;
            }
        }
        // has separate_loop case
        // exx converged or get max exx steps
        else if (exx_lri->two_level_step == GlobalC::exx_info.info_global.hybrid_step
                 || (iter == 1 && exx_lri->two_level_step != 0))
        {
            return true;
        }
        else
        {
            // update exx and redo scf
            if (exx_lri->two_level_step == 0)
            {
                add_exx_operator();
                XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);
            }

			const bool flag_restart = (exx_lri->two_level_step==0) ? true : false;
			if (GlobalV::GAMMA_ONLY_LOCAL)
				exx_lri->mix_DMk_2D.mix(loc.dm_gamma, flag_restart);
			else
				exx_lri->mix_DMk_2D.mix(loc.dm_k, flag_restart);

            // GlobalC::exx_lcao.cal_exx_elec(p_esolver->LOC, p_esolver->LOWF.wfc_k_grid);
            exx_lri->cal_exx_elec(*loc.LOWF->ParaV);
            iter = 0;
            std::cout << " Updating EXX and rerun SCF" << std::endl;
            exx_lri->two_level_step++;
            return false;
        }
    }
    return true;
}