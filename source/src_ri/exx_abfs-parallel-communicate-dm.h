#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_DM_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_DM_H

#include "exx_abfs-parallel.h"
#include "../module_base/matrix.h"
#include "../module_base/matrix_wrapper.h"
#include "abfs-vector3_order.h"
#include <vector>
#include <map>
#include <set>
#ifdef __MPI
#include "mpi.h"
#endif
#include <atomic>
#include <boost/dynamic_bitset.hpp>
#include <boost/mpi.hpp>

class Exx_Abfs::Parallel::Communicate::DM
{
public:

	void cal_DM( 
		const Abfs::Vector3_Order<int> &Born_von_Karman_period,
		const set<std::pair<size_t,size_t>> &H_atom_pairs_core,
        const double threshold,
        complex<double>*** wfc_k_grid,
        double*** DM,
        double** DM_R );

	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> DMr;
	
private:

	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> LOC_to_grid( 
		const Abfs::Vector3_Order<int> &Born_von_Karman_period,
		const double threshold,
        double*** DM,
        double** DM_R) const;

	class Allreduce
	{
	public:
		Allreduce( 
			const MPI_Comm & mpi_comm_in, 
			std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_local_in,
			const Abfs::Vector3_Order<int> &Born_von_Karman_period,
			const set<std::pair<size_t,size_t>> &H_atom_pairs_core);
		~Allreduce();
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> grid_to_exx();
		
	private:
		void ask( const int rank_delta_now );
		void recv_data_process( const int rank_data );
		void insert_data( std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_rank );
		void insert_data( std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> &data_rank );
		void send_data_process( const int rank_asked );
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> get_data_local_wrapper( const std::map<size_t,std::map<size_t,set<Abfs::Vector3_Order<int>>>> & atom_asked ) const;
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> get_data_local_wrapper() const;

	private:
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> data_all;
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_local;
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> data_localw;
		
		std::map<size_t,std::map<size_t,set<Abfs::Vector3_Order<int>>>> atom_unset;
		
		const MPI_Comm & mpi_comm;
		int comm_sz;
		int my_rank;
		
		static constexpr int tag_ask  = 1;
		static constexpr int tag_data = 2;
		
		atomic_flag lock_insert;
		atomic<int> lock_atom_unset_read;
		atomic<int> rank_delta;
		
		std::vector<boost::mpi::packed_oarchive*> oarps_isend_data;
		std::vector<boost::mpi::packed_oarchive*> oarps_atom_unset;
		std::vector<boost::mpi::packed_iarchive*> iarps_recv_data;
		std::vector<boost::mpi::packed_iarchive*> iarps_atom_asked;
				
		std::vector<atomic<int>*> flags_isend_data;
		std::vector<atomic<int>*> flags_ask_atom;
		boost::dynamic_bitset<> flags_recv_data;
	};
		
		
/*
private:
	struct
	{
		std::vector<bool> row;
		std::vector<bool> col;
	}atom_in_exx;
	void set_atom_in_exx( const set<std::pair<size_t,size_t>> &H_atom_pairs_core );	
*/
};

#endif