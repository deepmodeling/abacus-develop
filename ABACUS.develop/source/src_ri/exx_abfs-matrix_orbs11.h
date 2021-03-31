#ifndef EXX_ABFS_MATRIX_ORBS11_H
#define EXX_ABFS_MATRIX_ORBS11_H

#include<vector>
using std::vector;
#include<map>
using std::map;
#include<set>
using std::set;

#include "exx_abfs.h"
#include "src_lcao/ORB_table_phi.h"
#include "src_lcao/ORB_gaunt_table.h"
#include "src_lcao/center2_orb-orb11.h"

class LCAO_Orbitals;

class Exx_Abfs::Matrix_Orbs11
{
public:
	// mode:
	//    1: <lcaos|lcaos>
	//    2: <jYs|jYs>  <abfs|abfs>
	void init(
		const int mode,
		const double kmesh_times=1, 		// extend Kcut, keep dK
		const double rmesh_times=1);		// extend Rcut, keep dR

	void init_radial(
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_A,
		const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_B);
	void init_radial(
		const LCAO_Orbitals &orb_A,
		const LCAO_Orbitals &orb_B);

	void init_radial_table();
	void init_radial_table( const map<size_t,map<size_t,set<double>>> &Rs );		// unit: ucell.lat0

	matrix cal_overlap_matrix(
		const size_t TA,
		const size_t TB,
		const Vector3<double> &tauA,											// unit: ucell.lat0
		const Vector3<double> &tauB,											// unit: ucell.lat0
		const Element_Basis_Index::IndexLNM &index_r,
		const Element_Basis_Index::IndexLNM &index_c ) const;
	map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> cal_overlap_matrix(
		const Element_Basis_Index::IndexLNM &index_r,
		const Element_Basis_Index::IndexLNM &index_c ) const;

protected:
	ORB_table_phi MOT;
	ORB_gaunt_table MGT;
                                               
	map<size_t,                                // TA
		map<size_t,                            // TB
			map<int,                           // LA
				map<size_t,                    // NA
					map<int,                   // LB
						map<size_t,            // NB
							Center2_Orb::Orb11>>>>>> center2_orb11_s;
};

// this->center2_orb11_s[TA][TB][LA][NA][LB][NB]

#endif
