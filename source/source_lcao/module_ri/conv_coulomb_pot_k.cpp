#include "conv_coulomb_pot_k.h"
#include "../../source_base/constants.h"
#include "module_parameter/parameter.h"
#include "../../source_basis/module_ao/ORB_atomic_lm.h"
#include "../../source_pw/module_pwdft/global.h"

namespace Conv_Coulomb_Pot_K
{
	std::vector<double> cal_psi_fock_limits(
		const std::vector<double> & psif)
	{
		std::vector<double> psik2_ccp(psif.size());
		for( size_t ik=0; ik<psif.size(); ++ik )
			{ psik2_ccp[ik] = ModuleBase::FOUR_PI * psif[ik]; }
		return psik2_ccp;
	}

	// rongshi add 2022-07-27
	// Sphere truction -- Spencer
	std::vector<double> cal_psi_fock_spencer(
		const std::vector<double> &psif,
		const std::vector<double> &k_radial,
		const double rcut)
	{
		std::vector<double> psik2_ccp(psif.size());
		for (size_t ik = 0; ik < psif.size(); ++ik)
			{ psik2_ccp[ik] = ModuleBase::FOUR_PI * psif[ik] * (1 - std::cos(k_radial[ik] * rcut)); }
		return psik2_ccp;
	}


	std::vector<double> cal_psi_erfc_limits(
		const std::vector<double> & psif,
		const std::vector<double> & k_radial,
		const double erfc_omega)
	{
		std::vector<double> psik2_ccp(psif.size());
		for( size_t ik=0; ik<psif.size(); ++ik )
			{ psik2_ccp[ik] = ModuleBase::FOUR_PI * psif[ik] * (1-std::exp(-(k_radial[ik]*k_radial[ik])/(4*erfc_omega*erfc_omega))); }
		return psik2_ccp;
	}

	std::vector<double> cal_psi_erfc_spencer(
								const std::vector<double>& psif,
                                const std::vector<double>& k_radial,
                                const double erfc_omega,
                                const double rcut)
	{
		double eps = 1e-14;
		std::vector<double> psik2_ccp(psif.size());
		for (size_t ik = 0; ik < psif.size(); ++ik)
		{
	        double temp0 = std::cos(k_radial[ik] * rcut) * std::erfc(erfc_omega * rcut);
	        double temp1 = std::exp(-(k_radial[ik] * k_radial[ik]) / (4 * erfc_omega * erfc_omega));
			std::complex<double> temp2 = std::complex<double>(0, 0);
			std::complex<double> temp3 = std::complex<double>(0, 0);
			if (temp1 >= eps)
			{
	            temp2 = ModuleBase::ErrorFunc::erf(0.5 * (ModuleBase::IMAG_UNIT * k_radial[ik] + 2 * erfc_omega * erfc_omega * rcut)
													/ erfc_omega);
				temp3 = ModuleBase::NEG_IMAG_UNIT
	                    * ModuleBase::ErrorFunc::erfi(0.5 * k_radial[ik] / erfc_omega + ModuleBase::IMAG_UNIT * erfc_omega * rcut);
			}
	        std::complex<double> fock_part = -0.5 * (-2 + 2 * temp0 + temp1 * (temp2 + temp3));
        	psik2_ccp[ik] = ModuleBase::FOUR_PI * psif[ik] * fock_part.real();
		}
		return psik2_ccp;
	}


	template<>
	Numerical_Orbital_Lm cal_orbs_ccp<Numerical_Orbital_Lm>(
		const Numerical_Orbital_Lm &orbs,
		const std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>> &coulomb_param,
		const double rmesh_times)
	{
		std::vector<double> psik2_ccp(orbs.get_psif().size(), 0.0);

		for(const auto &param_list : coulomb_param)
		{
			switch(param_list.first)
			{
				case Conv_Coulomb_Pot_K::Coulomb_Type::Fock:
				{
					for(const auto &param : param_list.second)
					{
						if(param.at("singularity_correction") == "limits")
							{ psik2_ccp = psik2_ccp + std::stod(param.at("alpha")) * cal_psi_fock_limits( orbs.get_psif() ); }
						else if(param.at("singularity_correction") == "spencer" || param.at("singularity_correction") == "revised_spencer")
							{ psik2_ccp = psik2_ccp + std::stod(param.at("alpha")) * cal_psi_fock_spencer( orbs.get_psif(), orbs.get_k_radial(), std::stod(param.at("Rcut")) ); }
						else
							{ throw std::invalid_argument( "singularity_correction = " + param.at("singularity_correction") + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__) ); }
					}
					break;
				}
				case Conv_Coulomb_Pot_K::Coulomb_Type::Erfc:
				{
					for(const auto &param : param_list.second)
					{
						if(param.at("singularity_correction") == "limits")
							{ psik2_ccp = psik2_ccp + std::stod(param.at("alpha")) * cal_psi_erfc_limits( orbs.get_psif(), orbs.get_k_radial(), std::stod(param.at("omega")) ); }
						else if(param.at("singularity_correction") == "spencer" || param.at("singularity_correction") == "revised_spencer")
							{ psik2_ccp = psik2_ccp + std::stod(param.at("alpha")) * cal_psi_erfc_spencer( orbs.get_psif(), orbs.get_k_radial(), std::stod(param.at("omega")), std::stod(param.at("Rcut")) ); }
						else
							{ throw std::invalid_argument( "singularity_correction = " + param.at("singularity_correction") + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__) ); }
					}
					break;
				}
				default:
				{
					throw std::invalid_argument( std::string(__FILE__) + " line " + std::to_string(__LINE__) );
				}
			}
		}

		const double dr = orbs.get_rab().back();
		const int Nr = (static_cast<int>(orbs.getNr()*rmesh_times)) | 1;

		std::vector<double> rab(Nr);
		for( size_t ir=0; ir<std::min(orbs.getNr(),Nr); ++ir )
			{ rab[ir] = orbs.getRab(ir); }
		for( size_t ir=orbs.getNr(); ir<Nr; ++ir )
			{ rab[ir] = dr; }

		std::vector<double> r_radial(Nr);
		for( size_t ir=0; ir<std::min(orbs.getNr(),Nr); ++ir )
			{ r_radial[ir] = orbs.getRadial(ir); }
		for( size_t ir=orbs.getNr(); ir<Nr; ++ir )
			{ r_radial[ir] = orbs.get_r_radial().back() + (ir - orbs.getNr() + 1) * dr; }

		Numerical_Orbital_Lm orbs_ccp;
		orbs_ccp.set_orbital_info(
			orbs.getLabel(),
			orbs.getType(),
			orbs.getL(),
			orbs.getChi(),
			Nr,
			ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab),
			ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial),
			Numerical_Orbital_Lm::Psi_Type::Psik2,
			ModuleBase::GlobalFunc::VECTOR_TO_PTR(psik2_ccp),
			orbs.getNk(),
			orbs.getDk(),
			orbs.getDruniform(),
			false,
			true, PARAM.inp.cal_force);
		return orbs_ccp;
	}

	template<>
	double get_rmesh_proportion(
		const Numerical_Orbital_Lm &orbs,
		const double psi_threshold)
	{
		for(int ir=orbs.getNr()-1; ir>=0; --ir)
		{
			if(std::abs(orbs.getPsi(ir))>=psi_threshold)
				{ return static_cast<double>(ir)/orbs.getNr(); }
		}
		return 0.0;
	}

}
