#ifndef EXX_OPT_ORB_H
#define EXX_OPT_ORB_H

#include "../module_hamilt_general/module_xc/exx_info.h"
#include "../module_base/matrix.h"
#include "../module_base/element_basis_index.h"
#include "module_cell/klist.h"
#include "module_basis/module_ao/ORB_read.h"
#include <RI/global/Tensor.h>
#include <vector>
#include <map>
#include <set>

namespace Exx_Opt_Orb
{
	extern void generate_matrix(
		const Exx_Info::Exx_Info_Opt_ABFs &info,
		const K_Vectors &kv,
		const LCAO_Orbitals& orb);

	// private:
	extern std::vector<std::vector<RI::Tensor<double>>> cal_I(
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>> &ms,
		const size_t TA, const size_t IA, const size_t TB, const size_t IB );
	extern RI::Tensor<double> cal_proj_22(
		const RI::Tensor<double> & m_big,
		const std::vector<RI::Tensor<double>> & m_left,
		const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
		const std::vector<RI::Tensor<double>> & m_right );
	extern RI::Tensor<double> cal_proj_21(
		const RI::Tensor<double> & m_big,
		const std::vector<RI::Tensor<double>> & m_left,
		const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
		const std::vector<RI::Tensor<double>> & m_right );
	/*
	extern RI::Tensor<double> cal_proj_12(
		const RI::Tensor<double> & m_big,
		const std::vector<RI::Tensor<double>> & m_left,
		const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
		const std::vector<RI::Tensor<double>> & m_right );
	*/
	extern RI::Tensor<double> cal_proj_11(
		const RI::Tensor<double> & m_big,
		const std::vector<RI::Tensor<double>> & m_left,
		const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
		const std::vector<RI::Tensor<double>> & m_right );
    extern void print_matrix(
		const Exx_Info::Exx_Info_Opt_ABFs &info,
        const K_Vectors &kv,
        const std::string& file_name,
		const std::vector<RI::Tensor<double>> &matrix_Q,
		const std::vector<std::vector<RI::Tensor<double>>> &matrix_S,
		const RI::Tensor<double> &matrix_V,
		const size_t TA, const size_t IA, const size_t TB, const size_t IB,
        const std::vector<double>& orb_cutoff,
		const ModuleBase::Element_Basis_Index::Range &range_jles,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_jles);
	extern std::map<size_t,std::map<size_t,std::set<double>>> get_radial_R();
};
#endif
