//=======================
// AUTHOR : Peize Lin
// DATE :   2023-02-23
//=======================

#ifndef MATRIX_ORB22_H
#define MATRIX_ORB22_H

#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb22.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_basis/module_ao/ORB_table_phi.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"
#include "module_base/vector3.h"
#include "module_base/element_basis_index.h"

#include <RI/global/Tensor.h>

#include <vector>
#include <map>
#include <set>

class Matrix_Orbs22
{
public:
	// mode:
	//    1: <lcaos lcaos|lcaos lcaos>
	void init(
		const int mode,
		const double kmesh_times, 		// extend Kcut, keep dK
		const double rmesh_times);		// extend Rcut, keep dR

	void init_radial(
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A1,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A2,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B1,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B2);
	void init_radial(
		const LCAO_Orbitals &orb_A1,
		const LCAO_Orbitals &orb_A2,
		const LCAO_Orbitals &orb_B1,
		const LCAO_Orbitals &orb_B2);

	void init_radial_table();
	void init_radial_table( const std::map<size_t,std::map<size_t,std::set<double>>> &Rs );		// unit: ucell.lat0

	enum class Matrix_Order{ A1A2B1B2, A1A2B2B1, A1B1A2B2, A1B1B2A2, A1B2A2B1, A1B2B1A2, A2A1B1B2, A2A1B2B1, A2B1A1B2, A2B1B2A1, A2B2A1B1, A2B2B1A1, B1A1A2B2, B1A1B2A2, B1A2A1B2, B1A2B2A1, B1B2A1A2, B1B2A2A1, B2A1A2B1, B2A1B1A2, B2A2A1B1, B2A2B1A1, B2B1A1A2, B2B1A2A1 };

	template<typename Tdata>
	RI::Tensor<Tdata> cal_overlap_matrix(
		const size_t TA,
		const size_t TB,
		const ModuleBase::Vector3<double> &tauA,												// unit: ucell.lat0
		const ModuleBase::Vector3<double> &tauB,												// unit: ucell.lat0
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B1,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B2,
		const Matrix_Order &matrix_order) const;
	template<typename Tdata>
	std::array<RI::Tensor<Tdata>,3> cal_grad_overlap_matrix(
		const size_t TA,
		const size_t TB,
		const ModuleBase::Vector3<double> &tauA,												// unit: ucell.lat0
		const ModuleBase::Vector3<double> &tauB,												// unit: ucell.lat0
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B1,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B2,
		const Matrix_Order &matrix_order) const;

    template<typename Tdata>
    std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, RI::Tensor<Tdata>>>>> cal_overlap_matrix_all(
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A1, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A2, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B1, 
        const ModuleBase::Element_Basis_Index::IndexLNM& index_B2) const;
private:
	ORB_table_phi MOT;
	ORB_gaunt_table MGT;

	std::map<size_t,                                       // TA
		std::map<size_t,                                   // TB
			std::map<int,                                  // LA1
				std::map<size_t,                           // NA1
					std::map<int,                          // LA2
						std::map<size_t,                   // NA2
							std::map<int,                  // LB1
								std::map<size_t,           // NB1
									std::map<int,          // LB2
										std::map<size_t,   // NB2
											Center2_Orb::Orb22>>>>>>>>>> center2_orb22_s;
	// this->center2_orb22_s[TA][TB][LA1][NA1][LA2][NA2][LB1][NB1][LB2][NB2]
};

#include "Matrix_Orbs22.hpp"

#endif
