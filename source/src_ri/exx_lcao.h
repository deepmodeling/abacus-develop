#ifndef EXX_LCAO_H
#define EXX_LCAO_H

#include "abfs.h"
#include "abfs-vector3_order.h"
#include "exx_abfs-matrix_orbs11.h" 
#include "exx_abfs-matrix_orbs21.h" 
#include "exx_abfs-dm.h"
#include "exx_abfs-parallel-communicate-hexx.h"
#include "exx_abfs-screen-schwarz.h"
#include "exx_abfs-screen-cauchy.h"
#include "../module_base/element_basis_index.h"
#include "../module_xc/exx_global.h"
#include "src_lcao/local_orbital_charge.h"

#if EXX_DM==1
#include "exx_abfs-parallel-communicate-dm.h"
#elif EXX_DM==2
#include "exx_abfs-parallel-communicate-dm2.h"
#elif EXX_DM==3
#include "exx_abfs-parallel-communicate-dm3.h"
#endif

#include <set>
#include <vector>
#include <map>
#include <deque>
#include<memory>
#include<limits>
#include<omp.h>

class Exx_Lcao
{
public:
	struct{ std::string process; std::string thread; std::string matrix; } test_dir;	// Peize Lin test
	Exx_Lcao( const Exx_Global::Exx_Info &info_global );				// Peize Lin test
public:
	void init();
	void cal_exx_ions(const Parallel_Orbitals &pv);
	void cal_exx_elec(Local_Orbital_Charge &loc, complex<double>*** wfc_k_grid);
	void cal_exx_elec_nscf(const Parallel_Orbitals &pv);
	void add_Hexx(const size_t ik, const double alpha, LCAO_Matrix &lm) const;
private:
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> cal_Hexx() const;
	double cal_energy(
		const std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &HexxR ) const;
	void init_radial_table_ions( const std::set<size_t> &atom_centres_core, const std::vector<std::pair<size_t,size_t>> &atom_pairs_core );
		
public:
	enum class Distribute_Type {Htime,Kmeans2,Kmeans1,Order};
	
public:
	double get_energy() const { return energy; }
	
private:
	
	int kmesh_times = 4;				// Peize Lin test
	
#if EXX_DM==1
	Exx_Abfs::Parallel::Communicate::DM DM_para;
#elif EXX_DM==2
	Exx_Abfs::Parallel::Communicate::DM2 DM_para;
#elif EXX_DM==3
	Exx_Abfs::Parallel::Communicate::DM3 DM_para;
#endif

#ifdef __MPI
	Exx_Abfs::Parallel::Communicate::Hexx Hexx_para;
#endif
	double energy = 0.0;
	
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> lcaos;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_ccp;
	ModuleBase::Element_Basis_Index::IndexLNM index_lcaos;
	ModuleBase::Element_Basis_Index::IndexLNM index_abfs;

	Exx_Abfs::Matrix_Orbs11 m_abfs_abfs;
	Exx_Abfs::Matrix_Orbs21 m_abfslcaos_lcaos;	
	
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> Cws;		// Cws[it1][it2][R]
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> Vws;		// Vws[it1][it2][R]
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> Cs;		// Cs[iat1][iat2][box2]
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> Cps;		// Cs[iat1][iat2][box2]
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> Vs;		// Vs[iat1][iat2][box2]
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> Vps;		// Vps[iat1][iat2][box2]

	std::vector<std::pair<size_t,size_t>> atom_pairs_core_origin;
	std::vector<std::pair<size_t,size_t>> atom_pairs_core;
	std::set<std::pair<size_t,size_t>> H_atom_pairs_core;
	Abfs::Vector3_Order<int> Born_von_Karman_period;

	Exx_Abfs::Screen::Schwarz schwarz;
	Exx_Abfs::Screen::Cauchy cauchy;
		
public:
	struct Exx_Info
	{
		const Exx_Global::Hybrid_Type &hybrid_type;
		
		const double &hse_omega;
		
		double pca_threshold = 0;
		std::vector<std::string> files_abfs;
		double c_threshold  = 0;
		double v_threshold  = 0;
		double dm_threshold = 0;
		double schwarz_threshold = 0;
		double cauchy_threshold = 0;
		double ccp_threshold = 0;
		double ccp_rmesh_times = 10;
		
		Exx_Lcao::Distribute_Type distribute_type;

		Exx_Info( const Exx_Global::Exx_Info &info_global );
	};
	Exx_Info info;
	
	friend class Local_Orbital_Charge;
	friend class LCAO_Hamilt;
};

#endif	// EXX_LCAO_H
