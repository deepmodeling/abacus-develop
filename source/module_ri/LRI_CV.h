//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef LRI_CV_H
#define LRI_CV_H

#include "Matrix_Orbs11.h"
#include "Matrix_Orbs21.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "module_base/abfs-vector3_order.h"
#include "module_base/element_basis_index.h"

#include <RI/global/Tensor.h>
#include <RI/global/Global_Func-2.h>

#include <vector>
#include <map>
#include <functional>
#include <pthread.h>

template<typename Tdata>
class LRI_CV
{
private:
	using TA = int;
	using TC = std::array<int,3>;
	using TAC = std::pair<TA,TC>;
	using Tdata_real = RI::Global_Func::To_Real_t<Tdata>;

public:
	LRI_CV();
	~LRI_CV();

	void set_orbitals(
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &lcaos_in,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs_in,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs_ccp_in,
		const double &kmesh_times,
		const double &ccp_rmesh_times_in);
	inline std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>
	cal_Vs(
		const std::vector<TA> &list_A0,
		const std::vector<TAC> &list_A1,
		const std::map<std::string,bool> &flags);						// "writable_Vws"
	inline std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>
	cal_dVs(
		const std::vector<TA> &list_A0,
		const std::vector<TAC> &list_A1,
		const std::map<std::string,bool> &flags);						// "writable_dVws"
	std::pair<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,
	          std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>>
	cal_Cs_dCs(
		const std::vector<TA> &list_A0,
		const std::vector<TAC> &list_A1,
		const std::map<std::string,bool> &flags);						// "cal_dC", "writable_Cws", "writable_dCws", "writable_Vws", "writable_dVws"
	
	size_t get_index_abfs_size(const size_t &iat){return this->index_abfs[iat].count_size; }

private:
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> lcaos;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_ccp;
	ModuleBase::Element_Basis_Index::IndexLNM index_lcaos;
	ModuleBase::Element_Basis_Index::IndexLNM index_abfs;
	double ccp_rmesh_times;

public:
	std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,RI::Tensor<Tdata>>>> Vws;
	std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,RI::Tensor<Tdata>>>> Cws;
	std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,std::array<RI::Tensor<Tdata>,3>>>> dVws;
	std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,std::array<RI::Tensor<Tdata>,3>>>> dCws;
private:
	pthread_rwlock_t rwlock_Vw;
	pthread_rwlock_t rwlock_Cw;
	pthread_rwlock_t rwlock_dVw;
	pthread_rwlock_t rwlock_dCw;

	Matrix_Orbs11 m_abfs_abfs;
	Matrix_Orbs21 m_abfslcaos_lcaos;

	template<typename Tresult>
	using T_func_DPcal_data = std::function<Tresult(
		const int it0,
		const int it1,
		const Abfs::Vector3_Order<double> &R,
		const std::map<std::string,bool> &flags)>;
	template<typename Tresult>
	std::map<TA,std::map<TAC,Tresult>>
	cal_datas(
		const std::vector<TA> &list_A0,
		const std::vector<TAC> &list_A1,
		const std::map<std::string,bool> &flags,
		const double &rmesh_times,
		const T_func_DPcal_data<Tresult> &func_DPcal_data);

	inline RI::Tensor<Tdata>
	DPcal_V(
		const int it0,
		const int it1,
		const Abfs::Vector3_Order<double> &R,
		const std::map<std::string,bool> &flags);						// "writable_Vws"
	inline std::array<RI::Tensor<Tdata>,3>
	DPcal_dV(
		const int it0,
		const int it1,
		const Abfs::Vector3_Order<double> &R,
		const std::map<std::string,bool> &flags);						// "writable_dVws"
	std::pair<RI::Tensor<Tdata>, std::array<RI::Tensor<Tdata>,3>>
	DPcal_C_dC(
		const int it0,
		const int it1,
		const Abfs::Vector3_Order<double> &R,
		const std::map<std::string,bool> &flags);						// "cal_dC", "writable_Cws", "writable_dCws", "writable_Vws", "writable_dVws"

	template<typename To11, typename Tfunc>
	To11 DPcal_o11(
		const int it0,
		const int it1,
		const Abfs::Vector3_Order<double> &R,
		const bool &flag_writable_o11ws,
		pthread_rwlock_t &rwlock_o11,
		std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,To11>>> &o11ws,
		const Tfunc &func_cal_o11);
};

#include "LRI_CV.hpp"

#endif
