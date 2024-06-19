#ifndef OPEXXLCAO_HPP
#define OPEXXLCAO_HPP

#ifdef __EXX

#include "op_exx_lcao.h"
#include "module_ri/RI_2D_Comm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

namespace hamilt
{

template <typename TK, typename TR>
OperatorEXX<OperatorLCAO<TK, TR>>::OperatorEXX(
	LCAO_Matrix* LM_in,
	hamilt::HContainer<TR>* hR_in,
	std::vector<TK>* hK_in,
	const K_Vectors& kv_in,
	std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd_in,
	std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc_in,
	int* two_level_step_in,
	const bool restart_in)
	: OperatorLCAO<TK, TR>(LM_in, kv_in.kvec_d, hR_in, hK_in),
	  kv(kv_in),
	  Hexxd(Hexxd_in),
	  Hexxc(Hexxc_in), 
	  two_level_step(two_level_step_in),
	  restart(restart_in)
{
	this->cal_type = calculation_type::lcao_exx;
	if (this->restart)
	{///  Now only Hexx depends on DM, so we can directly read Hexx to reduce the computational cost.
	/// If other operators depends on DM, we can also read DM and then calculate the operators to save the memory to store operator terms.
		assert(this->two_level_step != nullptr);
		/// read in Hexx
		if (std::is_same<TK, double>::value)
		{
			this->LM->Hexxd_k_load.resize(this->kv.get_nks());
			for (int ik = 0; ik < this->kv.get_nks(); ik++)
			{
				this->LM->Hexxd_k_load[ik].resize(this->LM->ParaV->get_local_size(), 0.0);
				this->restart = GlobalC::restart.load_disk(
					"Hexx", ik,
					this->LM->ParaV->get_local_size(), this->LM->Hexxd_k_load[ik].data(), false);
				if (!this->restart) break;
			}
		}
		else
		{
			this->LM->Hexxc_k_load.resize(this->kv.get_nks());
			for (int ik = 0; ik < this->kv.get_nks(); ik++)
			{
				this->LM->Hexxc_k_load[ik].resize(this->LM->ParaV->get_local_size(), 0.0);
				this->restart = GlobalC::restart.load_disk(
					"Hexx", ik,
					this->LM->ParaV->get_local_size(), this->LM->Hexxc_k_load[ik].data(), false);
				if (!this->restart) break;
			}
		}
		if (!this->restart)
			std::cout << "WARNING: Hexx not found, restart from the non-exx loop." << std::endl
			          << "If the loaded charge density is EXX-solved, this may lead to poor convergence." << std::endl;
		GlobalC::restart.info_load.load_H_finish = this->restart;
	}
}

template<typename TK, typename TR>
void OperatorEXX<OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
	// Peize Lin add 2016-12-03
	if (GlobalV::CALCULATION != "nscf" && this->two_level_step != nullptr && *this->two_level_step == 0 && !this->restart) return;  //in the non-exx loop, do nothing 
	if (XC_Functional::get_func_type() == 4 || XC_Functional::get_func_type() == 5)
	{
		if (this->restart && this->two_level_step != nullptr)
		{
			if (*this->two_level_step == 0)
			{
				this->add_loaded_Hexx(ik);
				return;
			}
			else // clear loaded Hexx and release memory
			{
				if (this->LM->Hexxd_k_load.size() > 0)
				{
					this->LM->Hexxd_k_load.clear();
					this->LM->Hexxd_k_load.shrink_to_fit();
				}
				else if (this->LM->Hexxc_k_load.size() > 0)
				{
					this->LM->Hexxc_k_load.clear();
					this->LM->Hexxc_k_load.shrink_to_fit();
				}
			}
		}
		// cal H(k) from H(R) normally

		if (GlobalC::exx_info.info_ri.real_number)
			RI_2D_Comm::add_Hexx(
				this->kv,
				ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				this->Hexxd == nullptr ? *this->LM->Hexxd : *this->Hexxd,
				*this->LM->ParaV,
				*this->hK);
		else
			RI_2D_Comm::add_Hexx(
				this->kv,
				ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				this->Hexxc == nullptr ? *this->LM->Hexxc : *this->Hexxc,
				*this->LM->ParaV,
				*this->hK);
	}
}

} // namespace hamilt
#endif // __EXX
#endif // OPEXXLCAO_HPP