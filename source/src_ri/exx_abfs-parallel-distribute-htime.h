#ifndef EXX_ABFS_PARALLEL_DISTRIBUTE_HTIME_H
#define EXX_ABFS_PARALLEL_DISTRIBUTE_HTIME_H

#include "exx_abfs-parallel.h"

#include "abfs.h"

#include <vector>
#include <utility>
#ifdef __MPI
class Exx_Abfs::Parallel::Distribute::Htime
{
public:
	static std::vector<std::pair<size_t,size_t>> distribute( 
		const Abfs::Vector3_Order<int> & Born_von_Karman_period,
		const double rmesh_times );
private:
	static std::vector<size_t> cal_Nadj( 
		const Abfs::Vector3_Order<int> & Born_von_Karman_period );
	static std::vector<std::pair<size_t,std::pair<size_t,size_t>>> cal_pair_costs( 
		const std::vector<size_t> &Nadj,
		const double rmesh_times );
	static std::vector<std::vector<std::pair<size_t,size_t>>> cal_rank_work( 
		const std::vector<std::pair<size_t,std::pair<size_t,size_t>>> & pair_costs );
};
#endif
#endif