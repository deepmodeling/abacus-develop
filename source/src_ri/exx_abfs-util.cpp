#include "exx_abfs-util.h"
#include <numeric>
#include <deque>

#ifdef __MPI
void Exx_Abfs::Util::bcast( std::vector<std::string> &v, const int rank_src, const MPI_Comm &mpi_comm )
{
	int my_rank;	MPI_Comm_rank( mpi_comm, &my_rank );
	if(my_rank==rank_src)
	{
		int size_v = v.size();
		MPI_Bcast( &size_v, 1, MPI_INT, rank_src, mpi_comm );
		
		std::vector<int> size_s(size_v);
		transform( v.begin(), v.end(), size_s.begin(), [](const std::string &s){return s.size();} );
		MPI_Bcast( ModuleBase::GlobalFunc::VECTOR_TO_PTR(size_s), size_s.size(), MPI_INT, rank_src, mpi_comm );
		
		std::string s_all = accumulate( v.begin(), v.end(), std::string() );
		MPI_Bcast( const_cast<char*>(s_all.c_str()), s_all.size(), MPI_CHAR, rank_src, mpi_comm );
	}
	else
	{
		int size_v;
		MPI_Bcast( &size_v, 1, MPI_INT, rank_src, mpi_comm );
		v.resize(size_v);
		
		std::vector<int> size_s(size_v);
		MPI_Bcast( ModuleBase::GlobalFunc::VECTOR_TO_PTR(size_s), size_s.size(), MPI_INT, rank_src, mpi_comm );
		
		const int size_all = accumulate( size_s.begin(), size_s.end(), 0 );
		char *s_all_tmp = new char[size_all];
		MPI_Bcast( s_all_tmp, size_all, MPI_CHAR, rank_src, mpi_comm );
		std::string s_all(s_all_tmp, size_all);
		delete[]s_all_tmp;
		
		std::deque<int> index_begin(size_v);
		partial_sum( size_s.begin(), size_s.end(), index_begin.begin() );
		index_begin.push_front(0);
		
		for(int i=0; i<size_v; ++i)
			v[i] = s_all.substr(index_begin[i], size_s[i]);
	}
}
#endif