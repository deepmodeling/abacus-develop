#ifndef EXX_ABFS_MATRIX_ORBS11_H
#define EXX_ABFS_MATRIX_ORBS11_H

#include <vector>
#include <map>
#include <set>

#include "exx_abfs.h"
#include "module_orbital/ORB_table_phi.h"
#include "module_orbital/ORB_gaunt_table.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb11.h"

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
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B);
	void init_radial(
		const LCAO_Orbitals &orb_A,
		const LCAO_Orbitals &orb_B);

	void init_radial_table();
	void init_radial_table( const std::map<size_t,std::map<size_t,std::set<double>>> &Rs );		// unit: ucell.lat0

	ModuleBase::matrix cal_overlap_matrix(
		const size_t TA,
		const size_t TB,
		const ModuleBase::Vector3<double> &tauA,											// unit: ucell.lat0
		const ModuleBase::Vector3<double> &tauB,											// unit: ucell.lat0
		const ModuleBase::Element_Basis_Index::IndexLNM &index_r,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_c ) const;
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> cal_overlap_matrix(
		const ModuleBase::Element_Basis_Index::IndexLNM &index_r,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_c ) const;

protected:
	ORB_table_phi MOT;
	ORB_gaunt_table MGT;
                                               
	std::map<size_t,                                // TA
		std::map<size_t,                            // TB
			std::map<int,                           // LA
				std::map<size_t,                    // NA
					std::map<int,                   // LB
						std::map<size_t,            // NB
							Center2_Orb::Orb11>>>>>> center2_orb11_s;
};

// this->center2_orb11_s[TA][TB][LA][NA][LB][NB]

#endif
