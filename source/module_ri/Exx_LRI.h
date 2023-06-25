//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef EXX_LRI_H
#define EXX_LRI_H

#include "LRI_CV.h"
#include "module_hamilt_general/module_xc/exx_info.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "module_base/matrix.h"
#include "module_ri/Mix_DMk_2D.h"
#include <RI/physics/Exx.h>

#include <vector>
#include <array>
#include <map>
#include <deque>
#include <mpi.h>

	class Local_Orbital_Charge;
	class Parallel_Orbitals;
	
	template<typename Tdata>
	class RPA_LRI;

    template<typename Tdata>
    class Exx_LRI_Interface;

template<typename Tdata>
class Exx_LRI
{
private:
	using TA = int;
	using Tcell = int;
	static constexpr std::size_t Ndim = 3;
	using TC = std::array<Tcell,Ndim>;
	using TAC = std::pair<TA,TC>;
	using TatomR = std::array<double,Ndim>;		// tmp

public:
	Exx_LRI( const Exx_Info::Exx_Info_RI &info_in ) :info(info_in){}

	void init(const MPI_Comm &mpi_comm_in, const K_Vectors &kv_in);
	void cal_exx_force();
	void cal_exx_stress();

	std::vector< std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>> Hexxs;
	Tdata Eexx;
	ModuleBase::matrix force_exx;
	ModuleBase::matrix stress_exx;

private:
	const Exx_Info::Exx_Info_RI &info;
    MPI_Comm mpi_comm;
    const K_Vectors *p_kv;

	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> lcaos;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_ccp;

	LRI_CV<Tdata> cv;
	RI::Exx<TA,Tcell,Ndim,Tdata> exx_lri;

	void cal_exx_ions();
	void cal_exx_elec(const Parallel_Orbitals &pv);
	void post_process_Hexx( std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> &Hexxs_io ) const;
	Tdata post_process_Eexx( const Tdata &Eexx_in ) const;

    int two_level_step = 0;
    Mix_DMk_2D mix_DMk_2D;
    
	friend class RPA_LRI<Tdata>;
    friend class Exx_LRI_Interface<Tdata>;
};

#include "Exx_LRI.hpp"

#endif