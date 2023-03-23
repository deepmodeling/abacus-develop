#ifndef EXX_ABFS_MATRIX_ORBS22_H
#define EXX_ABFS_MATRIX_ORBS22_H

#include <map>
using std::map;

#include "exx_abfs.h"
#include "module_basis/module_ao/ORB_table_phi.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb22.h"

class LCAO_Orbitals;

class Exx_Abfs::Matrix_Orbs22
{
public:
	// mode: 
	//    1: <lcaos lcaos|lcaos lcaos>
	void init(
		const int mode, 
		const double kmesh_times=1, 		// extend Kcut, keep dK
		const double rmesh_times=1);		// extend Rcut, keep dR

	void init_radial( 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A1, 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A2, 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B1, 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B2 );		
	void init_radial( 
		const LCAO_Orbitals &orb_A1, 
		const LCAO_Orbitals &orb_A2, 
		const LCAO_Orbitals &orb_B1, 
		const LCAO_Orbitals &orb_B2 );
		
	void init_radial_table();
	void init_radial_table( const std::map<size_t,std::map<size_t,std::set<double>>> &Rs );		// unit is ucell.lat0

	enum class Matrix_Order{A1B1_A2B2};
	ModuleBase::matrix cal_overlap_matrix(  
		const size_t TA, 
		const size_t TB, 
		const ModuleBase::Vector3<double> &tauA,
		const ModuleBase::Vector3<double> &tauB,  
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A1, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B1,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B2,
		const Matrix_Order &matrix_order) const;
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> cal_overlap_matrix(
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A1, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A2, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B1, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B2 ) const;

private:
	
	ORB_table_phi MOT;
	ORB_gaunt_table MGT;

	std::map<size_t,                                        // TA
		std::map<size_t,                                    // TB
			std::map<int,                                   // LA1
				std::map<size_t,                            // NA1
					std::map<int,                           // LA2
						std::map<size_t,                    // NA2
							std::map<int,                   // LB1
								std::map<size_t,            // NB1
									std::map<int,           // LB2
										std::map<size_t,    // NB2
											Center2_Orb::Orb22>>>>>>>>>> center2_orb22_s;
};
#endif
