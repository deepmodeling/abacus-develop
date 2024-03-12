//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef MATRIX_ORB11_H
#define MATRIX_ORB11_H

#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb11.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_basis/module_ao/ORB_table_phi.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"
#include "module_base/vector3.h"
#include "module_base/element_basis_index.h"

#include <RI/global/Tensor.h>

#include <vector>
#include <map>
#include <set>

class Matrix_Orbs11
{
public:
	// mode:
	//    1: <lcaos|lcaos>
	//    2: <jYs|jYs>  <abfs|abfs>
	void init(
		const int mode,
		const double kmesh_times,  		// extend Kcut, keep dK
		const double rmesh_times);		// extend Rcut, keep dR

	void init_radial(
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B);
	void init_radial(
		const LCAO_Orbitals &orb_A,
		const LCAO_Orbitals &orb_B);

	void init_radial_table();
	void init_radial_table( const std::map<size_t,std::map<size_t,std::set<double>>> &Rs );		// unit: ucell.lat0

	enum class Matrix_Order{ AB, BA };

	template<typename Tdata>
	RI::Tensor<Tdata> cal_overlap_matrix(
		const size_t TA,
		const size_t TB,
		const ModuleBase::Vector3<double> &tauA,												// unit: ucell.lat0
		const ModuleBase::Vector3<double> &tauB,												// unit: ucell.lat0
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B,
		const Matrix_Order &matrix_order) const;
	template<typename Tdata>
	std::array<RI::Tensor<Tdata>,3> cal_grad_overlap_matrix(
		const size_t TA,
		const size_t TB,
		const ModuleBase::Vector3<double> &tauA,												// unit: ucell.lat0
		const ModuleBase::Vector3<double> &tauB,												// unit: ucell.lat0
		const ModuleBase::Element_Basis_Index::IndexLNM &index_A,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_B,
		const Matrix_Order &matrix_order) const;

	template<typename Tdata>
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<Tdata>>>>> cal_overlap_matrix_all(
		const ModuleBase::Element_Basis_Index::IndexLNM &index_r,
        const ModuleBase::Element_Basis_Index::IndexLNM& index_c) const;
    
private:
	ORB_table_phi MOT;
	ORB_gaunt_table MGT;

	std::map<size_t,                                // TA
		std::map<size_t,                            // TB
			std::map<int,                           // LA
				std::map<size_t,                    // NA
					std::map<int,                   // LB
						std::map<size_t,            // NB
							Center2_Orb::Orb11>>>>>> center2_orb11_s;
	// this->center2_orb11_s[TA][TB][LA][NA][LB][NB]
};

#include "Matrix_Orbs11.hpp"

#endif
