//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef RI_2D_COMM_H
#define RI_2D_COMM_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_cell/klist.h"

#include <RI/global/Tensor.h>
#include <RI/ri/Cell_Nearest.h>
#include <array>
#include <deque>
#include <map>
#include <mpi.h>
#include <set>
#include <tuple>
#include <vector>

namespace RI_2D_Comm
{
using TA = int;
using Tcell = int;
static const size_t Ndim = 3;
using TC = std::array<Tcell, Ndim>;
using TAC = std::pair<TA, TC>;

// public:
template <typename Tdata, typename Tmatrix>
extern std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>> split_m2D_ktoR(
    const UnitCell& ucell,
	const K_Vectors& kv, 
    const std::vector<const Tmatrix*>& mks_2D, 
    const Parallel_2D& pv, 
    const int nspin, 
    const bool spgsym = false);

	// judge[is] = {s0, s1}
	extern std::vector<std::tuple<std::set<TA>, std::set<TA>>>
        get_2D_judge(const UnitCell& ucell, const Parallel_2D& pv);

    template<typename Tdata, typename TK>
    extern void add_Hexx(
        const UnitCell& ucell,
        const K_Vectors& kv,
        const int ik,
        const double alpha,
        const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>& Hs,
        const Parallel_Orbitals& pv,
        TK* hk);
	
		
	template <typename Tdata, typename TK>
	extern void  add_Hexx_td(
		const UnitCell& ucell,
		const K_Vectors& kv,
		const int ik,
		const double alpha,
		const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>& Hs,
		const Parallel_Orbitals& pv,
		const ModuleBase::Vector3<double>& At,
		TK* hk);

    template<typename Tdata, typename TR>
    extern void add_HexxR(
        const int current_spin,
        const double alpha,
                      const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>& Hs,
                      const Parallel_Orbitals& pv,
                      const int npol,
                      hamilt::HContainer<TR>& HlocR,
                      const RI::Cell_Nearest<int, int, 3, double, 3>* const cell_nearest = nullptr);

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
	extern inline std::tuple<int,int,int> get_iat_iw_is_block(const UnitCell& ucell,const int& iwt);
	extern inline int get_is_block(const int is_k, const int is_row_b, const int is_col_b);
	extern inline std::tuple<int,int> split_is_block(const int is_b);
	extern inline int get_iwt(const UnitCell& ucell, const int iat, const int iw_b, const int is_b);

    template <typename TA, typename TAC, typename T>
    extern std::map<TA, std::map<TAC, T>> comm_map2_first(const MPI_Comm& mpi_comm,
                                                                const std::map<TA, std::map<TAC, T>>& Ds_in,
                                                                const std::set<TA>& s0,
                                                                const std::set<TA>& s1);
    template <typename TA, typename TAC, typename T, typename Tjudge>
    extern std::map<TA, std::map<TAC, T>> comm_map2(const MPI_Comm& mpi_comm,
                                                    const std::map<TA, std::map<TAC, T>>& Ds_in,
                                                    const Tjudge& judge);
    template <typename Tkey, typename Tvalue>
    extern void set_value_add(Tkey&& key, Tvalue&& value, std::map<Tkey, Tvalue>& data);
    template <typename Tkey0, typename Tkey1, typename Tvalue>
    extern void set_value_add(std::tuple<Tkey0, Tkey1>&& key,
                            Tvalue&& value,
                            std::map<Tkey0, std::map<Tkey1, Tvalue>>& data);
    template <typename Tkey, typename Tvalue>
    extern void add_datas(std::map<Tkey, Tvalue>&& data_local, std::map<Tkey, Tvalue>& data_recv);
    template <typename Tkey0, typename Tkey1, typename Tvalue>
    extern void add_datas(std::map<Tkey0, std::map<Tkey1, Tvalue>>&& data_local,
                        std::map<Tkey0, std::map<Tkey1, Tvalue>>& data_recv);
} // namespace RI_2D_Comm

#include "RI_2D_Comm.hpp"

#endif