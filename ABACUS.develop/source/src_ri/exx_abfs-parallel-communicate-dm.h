#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_DM_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_DM_H

#include "exx_abfs-parallel.h"
#include "module_base/matrix.h"
#include "src_global/matrix_wrapper.h"
#include "src_ri/abfs-vector3_order.h"
#include <vector>
#include <map>
#include <set>
#include <mpi.h>
#include <atomic>
#include <boost/dynamic_bitset.hpp>
#include <boost/mpi.hpp>

class Exx_Abfs::Parallel::Communicate::DM
{
public:

	void cal_DM( 
		const Abfs::Vector3_Order<int> &Born_von_Karman_period,
		const set<pair<size_t,size_t>> &H_atom_pairs_core,
		const double threshold );

	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> DMr;
	
private:

	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> LOC_to_grid( 
		const Abfs::Vector3_Order<int> &Born_von_Karman_period,
		const double threshold ) const;

	class Allreduce
	{
	public:
		Allreduce( 
			const MPI_Comm & mpi_comm_in, 
			vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local_in,
			const Abfs::Vector3_Order<int> &Born_von_Karman_period,
			const set<pair<size_t,size_t>> &H_atom_pairs_core);
		~Allreduce();
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> grid_to_exx();
		
	private:
		void ask( const int rank_delta_now );
		void recv_data_process( const int rank_data );
		void insert_data( vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_rank );
		void insert_data( vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> &data_rank );
		void send_data_process( const int rank_asked );
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> get_data_local_wrapper( const map<size_t,map<size_t,set<Abfs::Vector3_Order<int>>>> & atom_asked ) const;
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> get_data_local_wrapper() const;

	private:
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> data_all;
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local;
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> data_localw;
		
		map<size_t,map<size_t,set<Abfs::Vector3_Order<int>>>> atom_unset;
		
		const MPI_Comm & mpi_comm;
		int comm_sz;
		int my_rank;
		
		static constexpr int tag_ask  = 1;
		static constexpr int tag_data = 2;
		
		atomic_flag lock_insert;
		atomic<int> lock_atom_unset_read;
		atomic<int> rank_delta;
		
		vector<boost::mpi::packed_oarchive*> oarps_isend_data;
		vector<boost::mpi::packed_oarchive*> oarps_atom_unset;
		vector<boost::mpi::packed_iarchive*> iarps_recv_data;
		vector<boost::mpi::packed_iarchive*> iarps_atom_asked;
				
		vector<atomic<int>*> flags_isend_data;
		vector<atomic<int>*> flags_ask_atom;
		boost::dynamic_bitset<> flags_recv_data;
	};
		
		
/*
private:
	struct
	{
		vector<bool> row;
		vector<bool> col;
	}atom_in_exx;
	void set_atom_in_exx( const set<pair<size_t,size_t>> &H_atom_pairs_core );	
*/
};

#endif