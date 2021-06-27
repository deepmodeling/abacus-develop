#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_H

#include "exx_abfs-parallel.h"
#include "module_base/matrix.h"
#include "src_global/matrix_wrapper.h"
#include "module_base/complexmatrix.h"
#include "src_ri/abfs-vector3_order.h"
#include <vector>
#include <valarray>
#include <map>
#include <deque>
#include <mpi.h>
#include <atomic>

// mohan comment out 2021-02-06
//#include <boost/dynamic_bitset.hpp>
//##include <boost/mpi.hpp>

class Exx_Abfs::Parallel::Communicate::Hexx
{
public:
	void Rexx_to_Km2D(
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &HR_exx,
		const pair<bool,bool> &io_HR_a2D );

private:
//	map<size_t,map<size_t,matrix>> R_to_K(
//		map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>> & HR) const;
//	map<size_t,map<size_t,ComplexMatrix>> R_to_K(
//		const map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>> & HR,
//		const size_t ik) const;
//	template< typename T >
//	T a2D_to_m2D( const map<size_t,map<size_t,T>> & H_a2D ) const;
	template<typename T>
	inline T H_phase(matrix &&HR, const int ik, const Abfs::Vector3_Order<int> &box2) const;
	template<typename Tmatrix>
	Tmatrix Ra2D_to_Km2D(
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &HR_a2D, const int ik) const;
	template<typename Tmatrix>
	void Ra2D_to_Km2D_mixing(
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &HR_a2D,
		vector<Tmatrix> &HK_m2D,
		vector<deque<Tmatrix>> &HK_m2D_pulay_seq) const;
	template<typename T>
	T pulay_mixing( const T &H_pulay_old, deque<T> &H_seq, const T &H_new ) const;

	#if EXX_H_COMM==1
	class Allreduce
	{
	public:
		Allreduce(
			const MPI_Comm &mpi_comm_in,
			vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local_in );
		~Allreduce();
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> exx_to_a2D();

	private:
		void set_atom_in_2D();
		void ask( const int rank_delta_now );
		void recv_data_process( const int rank_data );
		void insert_data( vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> &data_rank );
		void insert_data();
		void send_data_process( const int rank_asked );
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,Matrix_Wrapper>>>> get_data_local_wrapper( 
			const pair<vector<bool>,vector<bool>> &atom_asked ) const;

	private:
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> data_all;
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local;

		pair< vector<bool>, vector<bool> > atom_in_2D;

		const MPI_Comm & mpi_comm;
		int comm_sz;
		int my_rank;

		static constexpr int tag_ask  = 1;
		static constexpr int tag_data = 2;

		atomic_flag lock_insert;
		atomic<int> rank_delta;

		boost::mpi::packed_oarchive *oarp_atom_in_2D;
		vector<boost::mpi::packed_oarchive*> oarps_isend_data;
		vector<boost::mpi::packed_iarchive*> iarps_recv_data;
		vector<boost::mpi::packed_iarchive*> iarps_atom_asked;

		vector<atomic<int>*> flags_isend_data;
		vector<atomic<int>*> flags_ask_atom;
		boost::dynamic_bitset<> flags_recv_data;
	};
	#elif EXX_H_COMM==2
	class Allreduce2
	{
	public:
		void init(
			const MPI_Comm &mpi_comm_in,
			const set<pair<size_t,size_t>> &H_atom_pairs_core);
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> exx_to_a2D(
			vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local) const;

	private:
		enum class Flag_Send {undo, begin_oar, finish_oar, begin_isend, finish_isend};
		enum class Flag_Recv {undo, begin_irecv, begin_iar, finish_iar};
		
		//vector<pair<vector<bool>,vector<bool>>> get_atom_in_2D_list() const;
		vector<size_t> get_send_size_list(
			const set<pair<size_t,size_t>> &H_atom_pairs_core,
			const vector<pair<vector<bool>,vector<bool>>> &atom_in_2D_list) const;
		void init_flags(
			vector<atomic<Flag_Send>> &flags_send,
			vector<atomic<Flag_Recv>> &flags_recv) const;
		bool finish_judge(
			const vector<atomic<Flag_Send>> &flags_send,
			const vector<atomic<Flag_Recv>> &flags_recv) const;
		bool memory_enough(
			const int rank_send_next,
			const vector<atomic<Flag_Send>> &flags_send) const;
		void send_data_process(
			const int rank_send_now,
			const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local,
			vector<valarray<double>> &oars_isend,
			vector<atomic<Flag_Send>> &flags_send) const;
		void recv_data_process(
			const int rank_recv,
			vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_all,
			vector<valarray<double>> &iarps_irecv,
			vector<atomic<Flag_Recv>> &flags_recv,
			atomic_flag &lock_insert) const;			
		void insert_data(
			vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_rank,
			vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_all) const;

	private:
		MPI_Comm mpi_comm;
		int comm_sz;
		int my_rank;
		
		vector<pair<vector<bool>,vector<bool>>> atom_in_2D_list;
		vector<size_t> send_size_list;
		size_t recv_size;
	};
	#endif

public:
	vector<matrix>        HK_Gamma_m2D;
	vector<ComplexMatrix> HK_K_m2D;
	vector<deque<matrix>>        HK_Gamma_m2D_pulay_seq;
	vector<deque<ComplexMatrix>> HK_K_m2D_pulay_seq;

	enum class Mixing_Mode{ No, Plain, Pulay };
	Mixing_Mode mixing_mode;
	double mixing_beta;
	
	Allreduce2 allreduce2;
};

#endif
