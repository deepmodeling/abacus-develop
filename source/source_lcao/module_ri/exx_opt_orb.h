#ifndef EXX_OPT_ORB_H
#define EXX_OPT_ORB_H

#include "../../source_hamilt/module_xc/exx_info.h"
#include "../../source_base/matrix.h"
#include "../../source_base/element_basis_index.h"
#include "source_cell/klist.h"
#include "source_basis/module_ao/ORB_read.h"
#include <RI/global/Tensor.h>
#include <vector>
#include <map>
#include <set>

class Exx_Opt_Orb
{
public:
	void generate_matrix(
		const Exx_Info::Exx_Info_Opt_ABFs &info,
		const K_Vectors &kv,
		const UnitCell &ucell,
		const LCAO_Orbitals &orb) const;
private:
	std::vector<std::vector<RI::Tensor<double>>> cal_I( 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>> &ms, 
		const size_t TA, const size_t IA, const size_t TB, const size_t IB ) const;
	RI::Tensor<double> cal_mul_22(
		const std::vector<RI::Tensor<double>> & m_left,
		const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
		const std::vector<RI::Tensor<double>> & m_right ) const;
	RI::Tensor<double> cal_mul_21(
		const std::vector<RI::Tensor<double>> & m_left,
		const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
		const std::vector<RI::Tensor<double>> & m_right ) const;
	RI::Tensor<double> cal_mul_12(
		const std::vector<RI::Tensor<double>> & m_left,
		const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
		const std::vector<RI::Tensor<double>> & m_right ) const;
	RI::Tensor<double> cal_mul_11(
		const std::vector<RI::Tensor<double>> & m_left,
		const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
		const std::vector<RI::Tensor<double>> & m_right ) const;
	void print_matrix(
		const Exx_Info::Exx_Info_Opt_ABFs &info,
		const UnitCell& ucell,
		const K_Vectors &kv,
		const int Lmax,
		const std::vector<std::size_t> &ecut_number,
		const std::string& file_name,
		const std::vector<RI::Tensor<double>> &matrix_Q, 
		const std::vector<std::vector<RI::Tensor<double>>> &matrix_S,
		const RI::Tensor<double> &matrix_V,
		const size_t TA, const size_t IA, const size_t TB, const size_t IB,
		const std::vector<double>& orb_cutoff,
		const ModuleBase::Element_Basis_Index::Range &range_jles, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_jles) const;
	std::map<size_t,std::map<size_t,std::set<double>>> get_radial_R(const UnitCell& ucell) const;
};
#endif
