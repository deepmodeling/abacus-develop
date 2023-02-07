#ifndef EXX_ABFS_MATRIX_LCAOSLCAOS_LCAOSLCAOS_H
#define EXX_ABFS_MATRIX_LCAOSLCAOS_LCAOSLCAOS_H

#include <vector>
#include <map>
#include <set>

#include "exx_abfs.h"
#include "module_orbital/ORB_table_phi.h"
#include "module_orbital/ORB_gaunt_table.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb22.h"

class LCAO_Orbitals;

class Exx_Abfs::Matrix_Lcaoslcaos_Lcaoslcaos
{
public:
	// mode: 
	//    1: <lcaos lcaos|lcaos lcaos>
	void init(
		const int mode, 
		const double kmesh_times=1, 		// extend Kcut, keep dK
		const double rmesh_times=1);		// extend Rcut, keep dR

	void init_radial( 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A, 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B );		
	void init_radial( 
		const LCAO_Orbitals &orb_A, 
		const LCAO_Orbitals &orb_B );
		
	void init_radial_table();
	void init_radial_table( std::map<size_t,std::map<size_t,set<double>>> &Rs );		// unit is GlobalC::ucell.lat0
	
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> cal_overlap_matrix(
		const ModuleBase::Element_Basis_Index::IndexLNM &index_r, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_c );

private:
	
	ORB_table_phi MOT;
	ORB_gaunt_table MGT;

	std::map<size_t,                                // TA
		std::map<size_t,                            // TB
			std::map<size_t,                        // LA
				std::map<size_t,                    // NA
					std::map<size_t,                // LB
						std::map<size_t,            // NB
							Center2_Orb::Orb22>>>>>> center2_orb22_s;
};
#endif
