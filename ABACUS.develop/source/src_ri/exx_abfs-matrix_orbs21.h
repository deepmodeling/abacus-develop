//=========================================================
//AUTHOR : Peize Lin 
//DATE : 2016-02-03
//=========================================================

#ifndef EXX_ABFS_MATRIX_ORBS21_H
#define EXX_ABFS_MATRIX_ORBS21_H 

#include<vector>
using std::vector;
#include<map>
using std::map;
#include<set>
using std::set;

#include "exx_abfs.h"
#include "src_lcao/ORB_table_phi.h"
#include "src_lcao/ORB_gaunt_table.h"
#include "src_lcao/center2_orb-orb21.h"

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
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_A1, 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_A2, 
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_B );
	void init_radial(
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_A1, 
		const LCAO_Orbitals &orb_A2, 
		const LCAO_Orbitals &orb_B );	
		
	void init_radial_table();
	void init_radial_table( const map<size_t,map<size_t,set<double>>> &Rs );		// unit: ucell.lat0
		
	enum class Matrix_Order{A2B_A1,BA2_A1};
	matrix cal_overlap_matrix(  
		const size_t TA, 
		const size_t TB, 
		const Vector3<double> &tauA,											// unit: ucell.lat0
		const Vector3<double> &tauB,											// unit: ucell.lat0  
		const Element_Basis_Index::IndexLNM &index_A1, 
		const Element_Basis_Index::IndexLNM &index_A2,
		const Element_Basis_Index::IndexLNM &index_B,
		const Matrix_Order &matrix_order) const;
	vector<matrix> cal_overlap_matrix(  
		const size_t TA, 
		const size_t TB, 
		const Vector3<double> &tauA,											// unit: ucell.lat0
		const Vector3<double> &tauB,											// unit: ucell.lat0  
		const Element_Basis_Index::IndexLNM &index_A1, 
		const Element_Basis_Index::IndexLNM &index_A2,
		const Element_Basis_Index::IndexLNM &index_B) const;
	map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> cal_overlap_matrix(  
		const Element_Basis_Index::IndexLNM &index_A1, 
		const Element_Basis_Index::IndexLNM &index_A2,
		const Element_Basis_Index::IndexLNM &index_B) const;
	
protected:
	ORB_table_phi MOT;
	ORB_gaunt_table MGT;
	
	map<size_t,                                  // TA
		map<size_t,                              // TB
			map<int,                             // LA1
				map<size_t,                      // NA1
					map<int,                     // LA2
						map<size_t,              // NA2
							map<int,             // LB
								map<size_t,      // NB
									Center2_Orb::Orb21>>>>>>>> center2_orb21_s;
};

// this->center2_orb21_s[TA][TB][LA1][NA1][LA2][NA2][LB][NB]

#endif
