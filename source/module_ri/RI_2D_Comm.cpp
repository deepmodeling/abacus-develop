//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#include "RI_2D_Comm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_cell/klist.h"

#include <string>
#include <stdexcept>

// judge[is] = {s0, s1}
auto RI_2D_Comm::get_2D_judge(const Parallel_Orbitals &pv)
-> std::vector<std::tuple<std::set<TA>, std::set<TA>>>
{
	ModuleBase::TITLE("RI_2D_Comm","get_2D_judge");

	const std::map<int,int> nspin_b = {{1,1}, {2,1}, {4,2}};

	std::vector<std::set<TA>> iat0_list(nspin_b.at(GlobalV::NSPIN));
    for (int iwt0_2D = 0; iwt0_2D < pv.row_set.size(); ++iwt0_2D)
	{
        const int iwt0 = pv.row_set[iwt0_2D];
		int iat0, iw0_b, is0_b;
		std::tie(iat0,iw0_b,is0_b) = RI_2D_Comm::get_iat_iw_is_block(iwt0);
		iat0_list[is0_b].insert(iat0);
	}

	std::vector<std::set<TA>> iat1_list(nspin_b.at(GlobalV::NSPIN));
    for (int iwt1_2D = 0; iwt1_2D < pv.col_set.size(); ++iwt1_2D)
	{
        const int iwt1 = pv.col_set[iwt1_2D];
		int iat1, iw1_b, is1_b;
		std::tie(iat1,iw1_b,is1_b) = RI_2D_Comm::get_iat_iw_is_block(iwt1);
		iat1_list[is1_b].insert(iat1);
	}

	std::vector<std::tuple<std::set<TA>, std::set<TA>>> judge(GlobalV::NSPIN);
	switch(GlobalV::NSPIN)
	{
		case 1:
			judge[0] = std::make_tuple( std::move(iat0_list[0]), std::move(iat1_list[0]) );
			break;
		case 2:
			judge[0] = judge[1] = std::make_tuple( std::move(iat0_list[0]), std::move(iat1_list[0]) );
			break;
		case 4:
			for(int is0_b=0; is0_b<2; ++is0_b)
				for(int is1_b=0; is1_b<2; ++is1_b)
				{
					const int is_b = RI_2D_Comm::get_is_block(-1, is0_b, is1_b);
					judge[is_b] = std::make_tuple( iat0_list[is0_b], iat1_list[is1_b] );
				}
			break;
		default:
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
	return judge;
}


std::vector<int>
RI_2D_Comm::get_ik_list(const K_Vectors &kv, const int is_k)
{
	std::vector<int> ik_list;
	for(int ik=0; ik<kv.nks; ++ik)
		if(kv.isk[ik]==is_k)
			ik_list.push_back(ik);
	return ik_list;
}