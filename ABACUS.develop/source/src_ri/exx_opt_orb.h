#ifndef EXX_OPT_ORB_H
#define EXX_OPT_ORB_H

#include "module_base/matrix.h"
#include "src_global/element_basis_index.h"
#include <vector>
#include <map>
#include <set>
using namespace std;

class Exx_Opt_Orb
{
public:
	void generate_matrix() const;
private:
	vector<vector<matrix>> cal_I( 
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &ms, 
		const size_t TA, const size_t IA, const size_t TB, const size_t IB ) const;
	matrix cal_proj( 
		const matrix & m_big, 
		const vector<matrix> & m_left, 
		const vector<vector<matrix>> & m_middle, 
		const vector<matrix> & m_right ) const;
	void print_matrix(
		const string &file_name,
		const vector<matrix> &matrix_Q, 
		const vector<vector<matrix>> &matrix_S,
		const matrix &matrix_V,
		const size_t TA, const size_t IA, const size_t TB, const size_t IB,
		const Element_Basis_Index::Range &range_jles, 
		const Element_Basis_Index::IndexLNM &index_jles) const;
	map<size_t,map<size_t,set<double>>> get_radial_R() const;
		
	int kmesh_times = 4;
};

#endif