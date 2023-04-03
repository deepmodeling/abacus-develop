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
auto RI_2D_Comm::split_m2D_ktoR(const std::vector<Tmatrix> &mks_2D, const Parallel_Orbitals &pv)
-> std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>>
{
	ModuleBase::TITLE("RI_2D_Comm","split_m2D_ktoR");
	ModuleBase::timer::tick("RI_2D_Comm", "split_m2D_ktoR");

	const TC period = RI_Util::get_Born_vonKarmen_period();
	const std::map<int,int> nspin_k = {{1,1}, {2,2}, {4,1}};
	const double SPIN_multiple = std::map<int,double>{{1,0.5}, {2,1}, {4,1}}.at(GlobalV::NSPIN);							// why?

	std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>> mRs_a2D(GlobalV::NSPIN);
	for(int is_k=0; is_k<nspin_k.at(GlobalV::NSPIN); ++is_k)
	{
		const std::vector<int> ik_list = RI_2D_Comm::get_ik_list(is_k);
		for(const TC &cell : RI_Util::get_Born_von_Karmen_cells(period))
		{
			RI::Tensor<Tdata> mR_2D;
			for(const int ik : ik_list)
			{
				using Tdata_m = typename Tmatrix::type;
				RI::Tensor<Tdata_m> mk_2D = RI_Util::Matrix_to_Tensor<Tdata_m>(mks_2D[ik]);
				const Tdata_m frac = SPIN_multiple
					* RI::Global_Func::convert<Tdata_m>( std::exp(
						- ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (GlobalC::kv.kvec_c[ik] * (RI_Util::array3_to_Vector3(cell)*GlobalC::ucell.latvec))));
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
	double mixing_beta = GlobalC::CHR_MIX.get_mixing_beta();
	if(GlobalC::CHR_MIX.get_mixing_mode() == "plain")
	{
		if(Hk_seq[ik].empty)
			Hk = RI_2D_Comm::Hexxs_to_Hk(pv, Hs, ik);
		else
		{
			std::vector<std::vector<Tdata>> Hk_seq_tmp = Hk_seq[ik][0];
			std::vector<std::vector<Tdata>> Hk_tmp = RI_2D_Comm::Hexxs_to_Hk(pv, Hs, ik);
			for(size_t iwt0=0; iwt0!=GlobalV::NLOCAL; ++iwt0)
				for(size_t iwt1=0; iwt1!=GlobalV::NLOCAL; ++iwt1)
					if(pv.in_this_processor(iwt0, iwt1))
						Hk[iwt0][iwt1] = (1-mixing_beta) * Hk_seq_tmp[iwt0][iwt1] + mixing_beta * Hk_tmp[iwt0][iwt1];
			Hk_seq[ik].pop_front();
		}
		Hk_seq[ik].emplace_back(Hk);
	}
	else
		throw std::invalid_argument("EXX only support plain mixing.");

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
std::vector<std::vector<Tdata>> RI_2D_Comm::Hexxs_to_Hk(const Parallel_Orbitals &pv, 
				const std::vector< std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>> &Hexxs,
				const int ik
				)
{
	ModuleBase::TITLE("Exx_LRI", "Hexxs_to_Hk");
	ModuleBase::timer::tick("Exx_LRI", "Hexxs_to_Hk");

	std::vector<std::vector<Tdata>> Hk;
	Hk.resize(GlobalV::NLOCAL);
	for(size_t ir=0; ir!=GlobalV::NLOCAL; ++ir)
		Hk[ir].resize(GlobalV::NLOCAL);

	const std::map<int, std::vector<int>> is_list = {{1,{0}}, {2,{GlobalC::kv.isk[ik]}}, {4,{0,1,2,3}}};
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
				const std::complex<double> frac = std::exp( ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (GlobalC::kv.kvec_c[ik] * (RI_Util::array3_to_Vector3(cell1)*GlobalC::ucell.latvec)));
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

	return Hk;
	ModuleBase::timer::tick("Exx_LRI", "Hexxs_to_Hk");
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
			return {0,0};
		case 4:
			return {is_b/2, is_b%2};
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