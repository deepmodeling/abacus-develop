#ifndef EXX_ABFS_MATRIX_ORBS22_H
#define EXX_ABFS_MATRIX_ORBS22_H

#include<map>
using std::map;

#include "exx_abfs.h"
#include "module_ORB/ORB_table_phi.h"
#include "module_ORB/ORB_gaunt_table.h"
#include "src_lcao/center2_orb-orb22.h"

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
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_A1, 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_A2, 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_B1, 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_B2 );		
	void init_radial( 
		const LCAO_Orbitals &orb_A1, 
		const LCAO_Orbitals &orb_A2, 
		const LCAO_Orbitals &orb_B1, 
		const LCAO_Orbitals &orb_B2 );
		
	void init_radial_table();
	void init_radial_table( const map<size_t,map<size_t,set<double>>> &Rs );		// unit is ucell.lat0

	enum class Matrix_Order{A1B1_A2B2};
	matrix cal_overlap_matrix(  
		const size_t TA, 
		const size_t TB, 
		const Vector3<double> &tauA,
		const Vector3<double> &tauB,  
		const Element_Basis_Index::IndexLNM &index_A1, 
		const Element_Basis_Index::IndexLNM &index_A2,
		const Element_Basis_Index::IndexLNM &index_B1,
		const Element_Basis_Index::IndexLNM &index_B2,
		const Matrix_Order &matrix_order) const;
	map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> cal_overlap_matrix(
		const Element_Basis_Index::IndexLNM &index_A1, 
		const Element_Basis_Index::IndexLNM &index_A2, 
		const Element_Basis_Index::IndexLNM &index_B1, 
		const Element_Basis_Index::IndexLNM &index_B2 ) const;

private:
	
	ORB_table_phi MOT;
	ORB_gaunt_table MGT;

	map<size_t,                                        // TA
		map<size_t,                                    // TB
			map<int,                                   // LA1
				map<size_t,                            // NA1
					map<int,                           // LA2
						map<size_t,                    // NA2
							map<int,                   // LB1
								map<size_t,            // NB1
									map<int,           // LB2
										map<size_t,    // NB2
											Center2_Orb::Orb22>>>>>>>>>> center2_orb22_s;
};
#endif
