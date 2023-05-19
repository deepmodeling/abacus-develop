#ifndef EXX_INFO_H
#define EXX_INFO_H

#include "xc_functional.h"
#include "module_ri/conv_coulomb_pot_k.h"

struct Exx_Info
{
	struct Exx_Info_Global
	{
		bool cal_exx = false;

		Conv_Coulomb_Pot_K::Ccp_Type ccp_type;
		double hybrid_alpha = 0.25;
		double hse_omega = 0.11;
		double mixing_beta_for_loop1 = 1.0;
		
		bool separate_loop = true;
		size_t hybrid_step = 1;
	};	
	Exx_Info_Global info_global;



	struct Exx_Info_Lip
	{
		const Conv_Coulomb_Pot_K::Ccp_Type &ccp_type;
		const double &hse_omega;
		double lambda;

		Exx_Info_Lip( const Exx_Info::Exx_Info_Global &info_global )
			:ccp_type(info_global.ccp_type),
			 hse_omega(info_global.hse_omega){}
	};
	Exx_Info_Lip info_lip;	



	struct Exx_Info_RI
	{
		const Conv_Coulomb_Pot_K::Ccp_Type &ccp_type;
		const double &hse_omega;
		
		bool real_number = false;
		
		double pca_threshold = 0;
		std::vector<std::string> files_abfs;
		double C_threshold  = 0;
		double V_threshold  = 0;
		double dm_threshold = 0;
		double cauchy_threshold = 0;
		double C_grad_threshold  = 0;
		double V_grad_threshold  = 0;
		double cauchy_force_threshold = 0;
		double cauchy_stress_threshold = 0;
		double ccp_threshold = 0;
		double ccp_rmesh_times = 10;
		double kmesh_times = 4;

		int abfs_Lmax = 0;			// tmp

		Exx_Info_RI( const Exx_Info::Exx_Info_Global &info_global )
			:ccp_type(info_global.ccp_type),
			 hse_omega(info_global.hse_omega){}
	};
	Exx_Info_RI info_ri;


	Exx_Info()
		:info_lip(this->info_global),
		 info_ri(this->info_global){}
};

#endif
