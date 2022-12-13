#ifndef SERIALIZATION_CEREAL_H
#define SERIALIZATION_CEREAL_H

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/functional.hpp>

#include "../module_base/vector3.h"
#include "../module_base/abfs-vector3_order.h"
#include "../module_base/matrix.h"



template<class Archive, typename T> void serialize( Archive & ar, Abfs::Vector3_Order<T> & v ){ ar(v.x); ar(v.y); ar(v.z); }

namespace ModuleBase
{
	template<class Archive, typename T> void serialize( Archive & ar, Vector3<T> & v ){ ar(v.x); ar(v.y); ar(v.z); }

	template<class Archive> void save( Archive & ar, const matrix & m )
	{
		ar(m.nr);	ar(m.nc);
		ar(cereal::binary_data(m.c, m.nr*m.nc*sizeof(double)));
	}
	template<class Archive> void load( Archive & ar, matrix & m )
	{
		int nr, nc;
		ar(nr);	ar(nc);
		m.create(nr,nc);
		ar(cereal::binary_data(m.c, m.nr*m.nc*sizeof(double)));
	}
}

#include <mpi.h>
#include <sstream>

namespace ModuleBase
{
	template<typename T>
	void bcast_data_cereal(T &data, const MPI_Comm &mpi_comm, const int &rank_bcast)
	{
		int my_rank;	MPI_Comm_rank( mpi_comm, &my_rank );
		if(my_rank==rank_bcast)
		{
			std::stringstream ss;
			{
				cereal::BinaryOutputArchive ar(ss);
				ar(data);
			}
			const int size = ss.str().size();
			MPI_Bcast( const_cast<int*>(&size), 1, MPI_INT, rank_bcast, mpi_comm );
			MPI_Bcast( const_cast<char*>(ss.str().c_str()), size, MPI_CHAR, rank_bcast, mpi_comm ); 
		}
		else
		{
			int size;
			MPI_Bcast( &size, 1, MPI_INT, rank_bcast, mpi_comm );
			std::vector<char> c(size);
			MPI_Bcast( c.data(), size, MPI_CHAR, rank_bcast, mpi_comm );   
			std::stringstream ss;  
			ss.rdbuf()->pubsetbuf(c.data(),size);
			{
				cereal::BinaryInputArchive ar(ss);
				ar(data);
			}
		}	
	}
}

#endif
