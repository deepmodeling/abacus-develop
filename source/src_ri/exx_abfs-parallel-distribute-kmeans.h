#ifndef EXX_ABFS_PARALLEL_DISTRIBUTE_KMEANS_H
#define EXX_ABFS_PARALLEL_DISTRIBUTE_KMEANS_H

#include "exx_abfs-parallel.h"
#include "../module_base/vector3.h"
#include <vector>
#include <utility>
#include "mpi.h"
#ifdef __MPI
class Exx_Abfs::Parallel::Distribute::Kmeans
{
public:

	static std::vector<std::pair<size_t,size_t>> distribute_kmeans2( const MPI_Comm & mpi_comm, const int multiple_core=1 );
	static std::vector<std::pair<size_t,size_t>> distribute_kmeans1( const MPI_Comm & mpi_comm, const double rmesh_times );

private:

	struct Cluster
	{
		ModuleBase::Vector3<double> tau;
		size_t size;
		ModuleBase::Vector3<double> tau_sum;
	};

	struct Atom
	{
		ModuleBase::Vector3<double> tau;
		int center;
		double distance;
	};
	
	static std::pair< std::vector<Exx_Abfs::Parallel::Distribute::Kmeans::Atom>, std::vector<Exx_Abfs::Parallel::Distribute::Kmeans::Cluster> > 
		cluster( const int Nc );
};
#endif
#endif
