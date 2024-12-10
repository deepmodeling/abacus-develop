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
#include <RI/physics/Exx.h>

#include <vector>
#include <array>
#include <map>
#include <deque>
#include <mpi.h>

#include "module_exx_symmetry/symmetry_rotation.h"

	class Parallel_Orbitals;
	
	template<typename T, typename Tdata>
	class RPA_LRI;

	template<typename T, typename Tdata>
	class Exx_LRI_Interface;

    namespace LR
    {
        template<typename T, typename TR>
        class ESolver_LR;

        template<typename T>
        class OperatorLREXX;
    }

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
    Exx_LRI(const Exx_Info::Exx_Info_RI& info_in) :info(info_in) {}
    Exx_LRI operator=(const Exx_LRI&) = delete;
    Exx_LRI operator=(Exx_LRI&&);

    void reset_Cs(const std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in) { this->exx_lri.set_Cs(Cs_in, this->info.C_threshold); }
    void reset_Vs(const std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in) { this->exx_lri.set_Vs(Vs_in, this->info.V_threshold); }

	void init(const MPI_Comm &mpi_comm_in, 
			  const UnitCell &ucell,
			  const K_Vectors &kv_in, 
			  const LCAO_Orbitals& orb);
	void cal_exx_force(const int& nat);
    void cal_exx_stress(const double& omega, const double& lat0);
	void cal_exx_ions(const UnitCell& ucell, const bool write_cv = false);
    void cal_exx_elec(const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>& Ds,
		const UnitCell& ucell,
        const Parallel_Orbitals& pv,
        const ModuleSymmetry::Symmetry_rotation* p_symrot = nullptr);
    std::vector<std::vector<int>> get_abfs_nchis() const;

	std::vector< std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>> Hexxs;
    double Eexx;
	ModuleBase::matrix force_exx;
	ModuleBase::matrix stress_exx;
    

private:
	const Exx_Info::Exx_Info_RI &info;
	MPI_Comm mpi_comm;
	const K_Vectors *p_kv = nullptr;
    std::vector<double> orb_cutoff_;

	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> lcaos;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_ccp;

	LRI_CV<Tdata> cv;
	RI::Exx<TA,Tcell,Ndim,Tdata> exx_lri;

	void post_process_Hexx( std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> &Hexxs_io ) const;
    double post_process_Eexx(const double& Eexx_in) const;

	friend class RPA_LRI<double, Tdata>;
	friend class RPA_LRI<std::complex<double>, Tdata>;
	friend class Exx_LRI_Interface<double, Tdata>;
	friend class Exx_LRI_Interface<std::complex<double>, Tdata>;
    friend class LR::ESolver_LR<double, double>;
    friend class LR::ESolver_LR<std::complex<double>, double>;
    friend class LR::OperatorLREXX<double>;
    friend class LR::OperatorLREXX<std::complex<double>>;
};

#include "Exx_LRI.hpp"

#endif
