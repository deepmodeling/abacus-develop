#ifndef EXX_OPT_ORB_H
#define EXX_OPT_ORB_H

#include "../module_base/matrix.h"
#include "../module_base/element_basis_index.h"
#include "module_cell/klist.h"
#include <vector>
#include <map>
#include <set>
#include <RI/global/Tensor.h>

class Exx_Opt_Orb
{
public:
	void generate_matrix(const K_Vectors &kv) const;
private:
	std::vector<std::vector<RI::Tensor<double>>> cal_I( 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>> &ms, 
		const size_t TA, const size_t IA, const size_t TB, const size_t IB ) const;
	RI::Tensor<double> cal_proj( 
		const RI::Tensor<double> & m_big, 
		const std::vector<RI::Tensor<double>> & m_left, 
		const std::vector<std::vector<RI::Tensor<double>>> & m_middle, 
		const std::vector<RI::Tensor<double>> & m_right ) const;
    void print_matrix(
        const K_Vectors &kv,
        const std::string& file_name,
		const std::vector<RI::Tensor<double>> &matrix_Q, 
		const std::vector<std::vector<RI::Tensor<double>>> &matrix_S,
		const RI::Tensor<double> &matrix_V,
		const size_t TA, const size_t IA, const size_t TB, const size_t IB,
		const ModuleBase::Element_Basis_Index::Range &range_jles, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_jles) const;
	std::map<size_t,std::map<size_t,std::set<double>>> get_radial_R() const;
		
	int kmesh_times = 4;
};
#endif