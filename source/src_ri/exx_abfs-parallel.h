#ifndef EXX_ABFS_PARALLEL_H
#define EXX_ABFS_PARALLEL_H

#include "exx_abfs.h"

#ifdef __MPI
class Exx_Abfs::Parallel
{
	public:
	
	class Distribute
	{
		public:
		class Htime;
		class Kmeans;
		class Order;
	};
	
	class Communicate
	{
		public:
		class Allreduce;
		class Function;
		class Hexx;
		#if EXX_DM==1
		class DM;
		#elif EXX_DM==2
		class DM2;
		#elif EXX_DM==3
		class DM3;
		#endif
	};
};
#endif

#endif