//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef RI_2D_COMM_HPP
#define RI_2D_COMM_HPP

#include "RI_2D_Comm.h"
#include "RI_Util.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/tool_title.h"
#include "module_base/timer.h"

#include <RI/global/Global_Func-2.h>

#include <cmath>
#include <string>
#include <stdexcept>

template<typename Tdata, typename Tmatrix>
auto RI_2D_Comm::split_m2D_ktoR(const K_Vectors &kv, const std::vector<Tmatrix> &mks_2D, const Parallel_Orbitals &pv)
-> std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>>
{
	ModuleBase::TITLE("RI_2D_Comm","split_m2D_ktoR");
	ModuleBase::timer::tick("RI_2D_Comm", "split_m2D_ktoR");

	const TC period = RI_Util::get_Born_vonKarmen_period(kv);
	const std::map<int,int> nspin_k = {{1,1}, {2,2}, {4,1}};
	const double SPIN_multiple = std::map<int,double>{{1,0.5}, {2,1}, {4,1}}.at(GlobalV::NSPIN);							// why?

	std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>> mRs_a2D(GlobalV::NSPIN);
	for(int is_k=0; is_k<nspin_k.at(GlobalV::NSPIN); ++is_k)
	{
		const std::vector<int> ik_list = RI_2D_Comm::get_ik_list(kv, is_k);
		for(const TC &cell : RI_Util::get_Born_von_Karmen_cells(period))
		{
			RI::Tensor<Tdata> mR_2D;
			for(const int ik : ik_list)
			{
				using Tdata_m = typename Tmatrix::type;
				RI::Tensor<Tdata_m> mk_2D = RI_Util::Matrix_to_Tensor<Tdata_m>(mks_2D[ik]);
				const Tdata_m frac = SPIN_multiple
					* RI::Global_Func::convert<Tdata_m>( std::exp(
						- ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (kv.kvec_c[ik] * (RI_Util::array3_to_Vector3(cell)*GlobalC::ucell.latvec))));
				if(mR_2D.empty())
					mR_2D = RI::Global_Func::convert<Tdata>(mk_2D * frac);
				else
					mR_2D = mR_2D + RI::Global_Func::convert<Tdata>(mk_2D * frac);
			}

			for(int iwt0_2D=0; iwt0_2D!=mR_2D.shape[0]; ++iwt0_2D)
			{
				const int iwt0 =
					ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER()
					? pv.MatrixInfo.col_set[iwt0_2D]
					: pv.MatrixInfo.row_set[iwt0_2D];
				int iat0, iw0_b, is0_b;
				std::tie(iat0,iw0_b,is0_b) = RI_2D_Comm::get_iat_iw_is_block(iwt0);
				const int it0 = GlobalC::ucell.iat2it[iat0];
				for(int iwt1_2D=0; iwt1_2D!=mR_2D.shape[1]; ++iwt1_2D)
				{
					const int iwt1 =
						ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER()
						? pv.MatrixInfo.row_set[iwt1_2D]
						: pv.MatrixInfo.col_set[iwt1_2D];
					int iat1, iw1_b, is1_b;
					std::tie(iat1,iw1_b,is1_b) = RI_2D_Comm::get_iat_iw_is_block(iwt1);
					const int it1 = GlobalC::ucell.iat2it[iat1];

					const int is_b = RI_2D_Comm::get_is_block(is_k, is0_b, is1_b);
					RI::Tensor<Tdata> &mR_a2D = mRs_a2D[is_b][iat0][{iat1,cell}];
					if(mR_a2D.empty())
						mR_a2D = RI::Tensor<Tdata>({static_cast<size_t>(GlobalC::ucell.atoms[it0].nw), static_cast<size_t>(GlobalC::ucell.atoms[it1].nw)});
					mR_a2D(iw0_b,iw1_b) = mR_2D(iwt0_2D, iwt1_2D);
				}
			}
		}
	}
	ModuleBase::timer::tick("RI_2D_Comm", "split_m2D_ktoR");
	return mRs_a2D;
}


template<typename Tdata>
void RI_2D_Comm::add_Hexx(
    const K_Vectors &kv,
    const int ik,
	const double alpha,
	const std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>> &Hs,
	const Parallel_Orbitals &pv,
	LCAO_Matrix &lm,
	std::vector<std::deque<std::vector<std::vector<Tdata>>>> &Hk_seq
	)
{
	ModuleBase::TITLE("RI_2D_Comm","add_Hexx");
	ModuleBase::timer::tick("RI_2D_Comm", "add_Hexx");

	// const std::map<int, std::vector<int>> is_list = {{1,{0}}, {2,{GlobalC::kv.isk[ik]}}, {4,{0,1,2,3}}};
	// for(const int is_b : is_list.at(GlobalV::NSPIN))
	// {
	// 	int is0_b, is1_b;
	// 	std::tie(is0_b,is1_b) = RI_2D_Comm::split_is_block(is_b);
	// 	for(const auto &Hs_tmpA : Hs[is_b])
	// 	{
	// 		const TA &iat0 = Hs_tmpA.first;
	// 		for(const auto &Hs_tmpB : Hs_tmpA.second)
	// 		{
	// 			const TA &iat1 = Hs_tmpB.first.first;
	// 			const TC &cell1 = Hs_tmpB.first.second;
	// 			const std::complex<double> frac = alpha
	// 				* std::exp( ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (GlobalC::kv.kvec_c[ik] * (RI_Util::array3_to_Vector3(cell1)*GlobalC::ucell.latvec)) );
	// 			const RI::Tensor<Tdata> &H = Hs_tmpB.second;
	// 			for(size_t iw0_b=0; iw0_b<H.shape[0]; ++iw0_b)
	// 			{
	// 				const int iwt0 = RI_2D_Comm::get_iwt(iat0, iw0_b, is0_b);
	// 				if(pv.trace_loc_row[iwt0]<0)	continue;
	// 				for(size_t iw1_b=0; iw1_b<H.shape[1]; ++iw1_b)
	// 				{
	// 					const int iwt1 = RI_2D_Comm::get_iwt(iat1, iw1_b, is1_b);
	// 					if(pv.trace_loc_col[iwt1]<0)	continue;

	// 					if(GlobalV::GAMMA_ONLY_LOCAL)
	// 						lm.set_HSgamma(iwt0, iwt1,
	// 							RI::Global_Func::convert<double>(H(iw0_b, iw1_b)) * RI::Global_Func::convert<double>(frac),
	// 							'L', lm.Hloc.data());
	// 					else
	// 						lm.set_HSk(iwt0, iwt1,
	// 							RI::Global_Func::convert<std::complex<double>>(H(iw0_b, iw1_b)) * frac,
	// 							'L', -1);
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	std::vector<std::vector<Tdata>> Hk;

	if(GlobalC::exx_info.info_global.separate_loop == 0)
	{
		double mixing_beta = GlobalC::CHR_MIX.get_mixing_beta();
		const std::string mixing_mode = GlobalC::CHR_MIX.get_mixing_mode();
		if(mixing_mode == "plain")
		{
			if(Hk_seq[ik].empty())
			{
				Hk = RI_2D_Comm::Hexxs_to_Hk(kv, pv, Hs, ik);
				Hk_seq[ik].emplace_back(Hk);
			}
			else
			{
				std::vector<std::vector<Tdata>> Hk_seq_tmp = Hk_seq[ik][0];
				std::vector<std::vector<Tdata>> Hk_tmp = RI_2D_Comm::Hexxs_to_Hk(kv, pv, Hs, ik);
				for(size_t iwt0=0; iwt0!=GlobalV::NLOCAL; ++iwt0)
					for(size_t iwt1=0; iwt1!=GlobalV::NLOCAL; ++iwt1)
						if(pv.in_this_processor(iwt0, iwt1))
							Hk[iwt0][iwt1] = (1-mixing_beta) * Hk_seq_tmp[iwt0][iwt1] + mixing_beta * Hk_tmp[iwt0][iwt1];
				Hk_seq[ik][0] = Hk;
			}
		}
		else if(mixing_mode == "pulay")
		{
			std::vector<std::vector<Tdata>> Hk_new = RI_2D_Comm::Hexxs_to_Hk(kv, pv, Hs, ik);
			Hk = RI_2D_Comm::pulay_mixing(pv, Hk_seq[ik], Hk_new, mixing_beta, mixing_mode);
		}
		else
			ModuleBase::WARNING_QUIT("RI_2D_Comm::add_Hexx",
    	                                 "`exx_separate_loop 0` only support plain and pulay mixing.");
	}
	else if(GlobalC::exx_info.info_global.separate_loop == 1)
	{
		double mixing_beta = GlobalC::exx_info.info_global.mixing_beta_for_loop1;
		if(Hk_seq[ik].empty())
		{
			Hk = RI_2D_Comm::Hexxs_to_Hk(kv, pv, Hs, ik);
			Hk_seq[ik].emplace_back(Hk);
		}
		else
		{
			std::vector<std::vector<Tdata>> Hk_seq_tmp = Hk_seq[ik][0];
			std::vector<std::vector<Tdata>> Hk_tmp = RI_2D_Comm::Hexxs_to_Hk(kv, pv, Hs, ik);
			for(size_t iwt0=0; iwt0!=GlobalV::NLOCAL; ++iwt0)
				for(size_t iwt1=0; iwt1!=GlobalV::NLOCAL; ++iwt1)
					if(pv.in_this_processor(iwt0, iwt1))
						Hk[iwt0][iwt1] = (1-mixing_beta) * Hk_seq_tmp[iwt0][iwt1] + mixing_beta * Hk_tmp[iwt0][iwt1];
			Hk_seq[ik][0] = Hk;
		}
	}

	for(size_t iwt0=0; iwt0!=GlobalV::NLOCAL; ++iwt0)
		for(size_t iwt1=0; iwt1!=GlobalV::NLOCAL; ++iwt1)
		{
			if(pv.in_this_processor(iwt0, iwt1))
			{
				const Tdata Hk_tmp = alpha * Hk[iwt0][iwt1];
				if(GlobalV::GAMMA_ONLY_LOCAL)
					lm.set_HSgamma(iwt0, iwt1, RI::Global_Func::convert<double>(Hk_tmp), 'L', lm.Hloc.data());
				else
					lm.set_HSk(iwt0, iwt1, RI::Global_Func::convert<std::complex<double>>(Hk_tmp), 'L', -1);
			}
		}

	ModuleBase::timer::tick("RI_2D_Comm", "add_Hexx");
}

template<typename Tdata>
std::vector<std::vector<Tdata>> RI_2D_Comm::Hexxs_to_Hk(
				const K_Vectors &kv,
				const Parallel_Orbitals &pv, 
				const std::vector< std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>> &Hexxs,
				const int ik
				)
{
	ModuleBase::TITLE("RI_2D_Comm", "Hexxs_to_Hk");
	ModuleBase::timer::tick("RI_2D_Comm", "Hexxs_to_Hk");

	std::vector<std::vector<Tdata>> Hk;
	Hk.resize(GlobalV::NLOCAL);
	for(size_t ir=0; ir!=GlobalV::NLOCAL; ++ir)
		Hk[ir].resize(GlobalV::NLOCAL);

	const std::map<int, std::vector<int>> is_list = {{1,{0}}, {2,{kv.isk[ik]}}, {4,{0,1,2,3}}};
	for(const int is_b : is_list.at(GlobalV::NSPIN))
	{
		int is0_b, is1_b;
		std::tie(is0_b, is1_b) = RI_2D_Comm::split_is_block(is_b);
		for(const auto &Hs_tmpA : Hexxs[is_b])
		{
			const TA &iat0 = Hs_tmpA.first;
			for(const auto &Hs_tmpB : Hs_tmpA.second)
			{
				const TA &iat1 = Hs_tmpB.first.first;
				const TC &cell1 = Hs_tmpB.first.second;
				const std::complex<double> frac = std::exp( ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (kv.kvec_c[ik] * (RI_Util::array3_to_Vector3(cell1)*GlobalC::ucell.latvec)));
				const RI::Tensor<Tdata> &H = Hs_tmpB.second;
				for(size_t iw0_b=0; iw0_b<H.shape[0]; ++iw0_b)
				{
					const int iwt0 = RI_2D_Comm::get_iwt(iat0, iw0_b, is0_b);
					if(pv.trace_loc_row[iwt0]<0) continue;
					for(size_t iw1_b=0; iw1_b<H.shape[1]; ++iw1_b)
					{
						const int iwt1 = RI_2D_Comm::get_iwt(iat1, iw1_b, is1_b);
						if(pv.trace_loc_col[iwt1]<0)	continue;
						Hk[iwt0][iwt1] += RI::Global_Func::convert<Tdata>(H(iw0_b, iw1_b)) * RI::Global_Func::convert<Tdata>(frac);
					}
				}
			}
		}
	}

	ModuleBase::timer::tick("RI_2D_Comm", "Hexxs_to_Hk");
	return Hk;
}

template<typename Tdata>
std::vector<std::vector<Tdata>> RI_2D_Comm::pulay_mixing(
	const Parallel_Orbitals &pv,
	std::deque<std::vector<std::vector<Tdata>>> &Hk_seq,
	const std::vector<std::vector<Tdata>> &Hk_new,
	const double mixing_beta,
	const std::string mixing_mode)
{
	ModuleBase::TITLE("RI_2D_Comm", "pulay_mixing");
	ModuleBase::timer::tick("RI_2D_Comm", "pulay_mixing");

	std::vector<std::vector<Tdata>> Hk;
	Hk.resize(GlobalV::NLOCAL);
	for(size_t ir=0; ir!=GlobalV::NLOCAL; ++ir)
		Hk[ir].resize(GlobalV::NLOCAL);

	if(GlobalC::CHR_MIX.get_totstep() == 0)
	{
		Hk_seq.clear();
		Hk = Hk_new;
		Hk_seq.emplace_back(Hk);
	}
	else if(GlobalC::CHR_MIX.get_totstep() == 1)
	{
		for(size_t iwt0=0; iwt0!=GlobalV::NLOCAL; ++iwt0)
			for(size_t iwt1=0; iwt1!=GlobalV::NLOCAL; ++iwt1)
				if(pv.in_this_processor(iwt0, iwt1))
				{
					Hk[iwt0][iwt1] = (1-mixing_beta) * Hk_seq[0][iwt0][iwt1] + mixing_beta * Hk_new[iwt0][iwt1];
					Hk_seq[0][iwt0][iwt1] = Hk[iwt0][iwt1];
				}
	}
	else
	{
		std::vector<std::vector<Tdata>> Hk_seq_tmp;
		Hk_seq_tmp.resize(GlobalV::NLOCAL);
		for(size_t ir=0; ir!=GlobalV::NLOCAL; ++ir)
			Hk_seq_tmp[ir].resize(GlobalV::NLOCAL);

		for(size_t iwt0=0; iwt0!=GlobalV::NLOCAL; ++iwt0)
			for(size_t iwt1=0; iwt1!=GlobalV::NLOCAL; ++iwt1)
				if(pv.in_this_processor(iwt0, iwt1))
					Hk_seq_tmp[iwt0][iwt1] = (1-mixing_beta) * Hk_seq[Hk_seq.size()-1][iwt0][iwt1] + mixing_beta * Hk_new[iwt0][iwt1];
		Hk_seq.emplace_back(Hk_seq_tmp);
		if(Hk_seq.size() > GlobalC::CHR_MIX.get_dstep()+1)
			Hk_seq.pop_front();
		
		auto pos_mod = [](const int i, const int N){ return (i%N+N)%N; };
		auto index = [&](const int i) -> int
		{
			const int alpha_size = Hk_seq.size()-1;
			const int alpha_begin = GlobalC::CHR_MIX.get_idstep() - alpha_size;
			return pos_mod( alpha_begin+i, GlobalC::CHR_MIX.get_dstep() );
		};
		for(size_t iwt0=0; iwt0!=GlobalV::NLOCAL; ++iwt0)
			for(size_t iwt1=0; iwt1!=GlobalV::NLOCAL; ++iwt1)
				if(pv.in_this_processor(iwt0, iwt1))
				{
					Hk[iwt0][iwt1] = (1+GlobalC::CHR_MIX.get_alpha()[index(Hk_seq.size()-2)]) * Hk_seq[Hk_seq.size()-1][iwt0][iwt1];
					for( int i=1; i<=Hk_seq.size()-2; ++i )
						Hk[iwt0][iwt1] += ( GlobalC::CHR_MIX.get_alpha()[index(i-1)] - GlobalC::CHR_MIX.get_alpha()[index(i)] ) * Hk_seq[i][iwt0][iwt1];
					Hk[iwt0][iwt1] -= GlobalC::CHR_MIX.get_alpha()[index(0)] * Hk_seq[0][iwt0][iwt1];
				}
	}

	ModuleBase::timer::tick("RI_2D_Comm", "pulay_mixing");
	return Hk;
}

std::tuple<int,int,int>
RI_2D_Comm::get_iat_iw_is_block(const int iwt)
{
	const int iat = GlobalC::ucell.iwt2iat[iwt];
	const int iw = GlobalC::ucell.iwt2iw[iwt];
	switch(GlobalV::NSPIN)
	{
		case 1: case 2:
			return std::make_tuple(iat, iw, 0);
		case 4:
			return std::make_tuple(iat, iw/2, iw%2);
		default:
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
}

int RI_2D_Comm::get_is_block(const int is_k, const int is_row_b, const int is_col_b)
{
	switch(GlobalV::NSPIN)
	{
		case 1:		return 0;
		case 2:		return is_k;
		case 4:		return is_row_b*2+is_col_b;
		default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
}

std::tuple<int,int>
RI_2D_Comm::split_is_block(const int is_b)
{
	switch(GlobalV::NSPIN)
	{
		case 1:	case 2:
			return std::make_tuple(0, 0);
		case 4:
			return std::make_tuple(is_b/2, is_b%2);
		default:
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
}



int RI_2D_Comm::get_iwt(const int iat, const int iw_b, const int is_b)
{
	const int it = GlobalC::ucell.iat2it[iat];
	const int ia = GlobalC::ucell.iat2ia[iat];
	int iw=-1;
	switch(GlobalV::NSPIN)
	{
		case 1: case 2:
			iw = iw_b;			break;
		case 4:
			iw = iw_b*2+is_b;	break;
		default:
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
	const int iwt = GlobalC::ucell.itiaiw2iwt(it,ia,iw);
	return iwt;
}

#endif
