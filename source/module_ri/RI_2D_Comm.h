//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef RI_2D_COMM_H
#define RI_2D_COMM_H

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_cell/klist.h"

#include <RI/global/Tensor.h>

#include <array>
#include <vector>
#include <tuple>
#include <map>
#include <set>
#include <deque>

namespace RI_2D_Comm
{
	using TA = int;
	using Tcell = int;
	static const size_t Ndim = 3;
	using TC = std::array<Tcell,Ndim>;
	using TAC = std::pair<TA,TC>;

//public:
	template<typename Tdata, typename Tmatrix>
	extern std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>>
	split_m2D_ktoR(const K_Vectors &kv, const std::vector<const Tmatrix*> &mks_2D, const Parallel_Orbitals &pv);

	// judge[is] = {s0, s1}
	extern std::vector<std::tuple<std::set<TA>, std::set<TA>>>
	get_2D_judge(const Parallel_Orbitals &pv);

	template<typename Tdata>
	extern void add_Hexx(
		const K_Vectors &kv,
		const int ik,
		const double alpha,
		const std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>> &Hs,
		LCAO_Matrix &lm);

	template<typename Tdata>
	extern std::vector<std::vector<Tdata>> Hexxs_to_Hk(
			const K_Vectors &kv,
			const Parallel_Orbitals &pv, 
			const std::vector< std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>> &Hexxs,
			const int ik);
	template<typename Tdata>
	std::vector<std::vector<Tdata>> pulay_mixing(
		const Parallel_Orbitals &pv,
		std::deque<std::vector<std::vector<Tdata>>> &Hk_seq,
		const std::vector<std::vector<Tdata>> &Hk_new,
		const double mixing_beta,
		const std::string mixing_mode);

//private:
	extern std::vector<int> get_ik_list(const K_Vectors &kv, const int is_k);
	extern inline std::tuple<int,int,int> get_iat_iw_is_block(const int iwt);
	extern inline int get_is_block(const int is_k, const int is_row_b, const int is_col_b);
	extern inline std::tuple<int,int> split_is_block(const int is_b);
	extern inline int get_iwt(const int iat, const int iw_b, const int is_b);
}

#include "RI_2D_Comm.hpp"

#endif