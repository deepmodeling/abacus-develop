//=========================================================
//AUTHOR : Peize Lin 
//DATE : 2016-02-03
//=========================================================

#ifndef EXX_ABFS_MATRIX_ORBS21_H
#define EXX_ABFS_MATRIX_ORBS21_H 

#include <vector>
#include <map>
#include <set>

#include "exx_abfs.h"
#include "module_orbital/ORB_table_phi.h"
#include "module_orbital/ORB_gaunt_table.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb21.h"

class LCAO_Orbitals;

class Exx_Abfs::Matrix_Orbs21
{	
public:
	// mode: 
	//    1: <jYs lcaos|lcaos>  <abfs lcaos|lcaos>
	void init(
		const int mode, 
		const double kmesh_times=1, 		// extend Kcut, keep dK
		const double rmesh_times=1);		// extend Rcut, keep dR
	
	void init_radial(	
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A1, 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A2, 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B );
	void init_radial(
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A1, 
		const LCAO_Orbitals &orb_A2, 
		const LCAO_Orbitals &orb_B );	
		
	void init_radial_table();
	void init_radial_table( const std::map<size_t,std::map<size_t,std::set<double>>> &Rs );		// unit: ucell.lat0
		
	enum class Matrix_Order{A2B_A1,BA2_A1};
	ModuleBase::matrix cal_overlap_matrix(  
		const size_t TA, 
		const size_t TB, 
		const ModuleBase::Vector3<double> &tauA,											// unit: ucell.lat0
		const ModuleBase::Vector3<double> &tauB,											// unit: ucell.lat0  
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A1, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B,
		const Matrix_Order &matrix_order) const;
	std::vector<ModuleBase::matrix> cal_overlap_matrix(  
		const size_t TA, 
		const size_t TB, 
		const ModuleBase::Vector3<double> &tauA,											// unit: ucell.lat0
		const ModuleBase::Vector3<double> &tauB,											// unit: ucell.lat0  
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A1, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B) const;
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> cal_overlap_matrix(  
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A1, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B) const;
	
protected:
	ORB_table_phi MOT;
	ORB_gaunt_table MGT;
	
	std::map<size_t,                                  // TA
		std::map<size_t,                              // TB
			std::map<int,                             // LA1
				std::map<size_t,                      // NA1
					std::map<int,                     // LA2
						std::map<size_t,              // NA2
							std::map<int,             // LB
								std::map<size_t,      // NB
									Center2_Orb::Orb21>>>>>>>> center2_orb21_s;
};

// this->center2_orb21_s[TA][TB][LA1][NA1][LA2][NA2][LB][NB]

#endif
