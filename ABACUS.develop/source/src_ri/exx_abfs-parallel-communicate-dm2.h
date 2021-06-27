#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_DM2_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_DM2_H

#include "exx_abfs-parallel.h"
#include "module_base/matrix.h"
#include "src_ri/abfs-vector3_order.h"
#include <vector>
#include <map>

class Exx_Abfs::Parallel::Communicate::DM2
{
public:

	void init( 
		const set<pair<size_t,size_t>> &H_atom_pairs_core,
		const double threshold_in);
	void clear_DMr();
	void set_DM_gamma( const matrix &DM_2D, const int is, const pair<int,int> &index_begin );
	
	void cal_DM_k( 
		const Abfs::Vector3_Order<int> &Born_von_Karman_period,
		const set<pair<size_t,size_t>> &H_atom_pairs_core,
		const double threshold );

	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> DMr;
	
private:

	double threshold;
	
	struct
	{
		vector<bool> row;
		vector<bool> col;
		vector<vector<bool>> row_col;
	}atom_in_exx;
};

#endif