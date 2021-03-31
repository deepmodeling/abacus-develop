#ifndef EXX_ABFS_MATRIX_LCAOSLCAOS_LCAOSLCAOS_H
#define EXX_ABFS_MATRIX_LCAOSLCAOS_LCAOSLCAOS_H

#include<vector>
using std::vector;
#include<map>
using std::map;
#include<set>
using std::set;

#include "exx_abfs.h"
#include "src_lcao/ORB_table_phi.h"
#include "src_lcao/ORB_gaunt_table.h"
#include "src_lcao/center2_orb-orb22.h"

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
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_A, 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_B );		
	void init_radial( 
		const LCAO_Orbitals &orb_A, 
		const LCAO_Orbitals &orb_B );
		
	void init_radial_table();
	void init_radial_table( map<size_t,map<size_t,set<double>>> &Rs );		// unit is ucell.lat0
	
	map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> cal_overlap_matrix(
		const Element_Basis_Index::IndexLNM &index_r, 
		const Element_Basis_Index::IndexLNM &index_c );

private:
	
	ORB_table_phi MOT;
	ORB_gaunt_table MGT;

	map<size_t,                                // TA
		map<size_t,                            // TB
			map<size_t,                        // LA
				map<size_t,                    // NA
					map<size_t,                // LB
						map<size_t,            // NB
							Center2_Orb::Orb22>>>>>> center2_orb22_s;
};
#endif
