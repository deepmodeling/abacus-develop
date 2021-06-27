#ifndef EXX_ABFS_PARALLEL_DISTRIBUTE_KMEANS_H
#define EXX_ABFS_PARALLEL_DISTRIBUTE_KMEANS_H

#include "exx_abfs-parallel.h"
#include "module_base/vector3.h"
#include <vector>
#include <utility>
#include <mpi.h>

class Exx_Abfs::Parallel::Distribute::Kmeans
{
public:

	static vector<pair<size_t,size_t>> distribute_kmeans2( const MPI_Comm & mpi_comm, const int multiple_core=1 );
	static vector<pair<size_t,size_t>> distribute_kmeans1( const MPI_Comm & mpi_comm, const double rmesh_times );

private:

	struct Cluster
	{
		Vector3<double> tau;
		size_t size;
		Vector3<double> tau_sum;
	};

	struct Atom
	{
		Vector3<double> tau;
		int center;
		double distance;
	};
	
	static pair< vector<Exx_Abfs::Parallel::Distribute::Kmeans::Atom>, vector<Exx_Abfs::Parallel::Distribute::Kmeans::Cluster> > 
		cluster( const int Nc );
};

#endif