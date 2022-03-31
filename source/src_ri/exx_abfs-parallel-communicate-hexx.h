#ifdef __MPI
#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_H

#include "exx_abfs-parallel.h"
#include "../module_base/matrix.h"
#include "../module_base/matrix_wrapper.h"
#include "../module_base/complexmatrix.h"
#include "abfs-vector3_order.h"
#include <vector>
#include <valarray>
#include <map>
#include <deque>
#include "mpi.h"
#include <atomic>
#include "module_orbital/parallel_orbitals.h"

// mohan comment out 2021-02-06
//#include <boost/dynamic_bitset.hpp>
//##include <boost/mpi.hpp>

class Exx_Abfs::Parallel::Communicate::Hexx
{
public:
	void Rexx_to_Km2D(const Parallel_Orbitals &pv, 
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &HR_exx,
		const std::pair<bool,bool> &io_HR_a2D );

private:
//	std::map<size_t,std::map<size_t,matrix>> R_to_K(
//		std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,matrix>>> & HR) const;
//	std::map<size_t,std::map<size_t,ModuleBase::ComplexMatrix>> R_to_K(
//		const std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,matrix>>> & HR,
//		const size_t ik) const;
//	template< typename T >
//	T a2D_to_m2D( const std::map<size_t,std::map<size_t,T>> & H_a2D ) const;
	template<typename T>
	inline T H_phase(ModuleBase::matrix &&HR, const int ik, const Abfs::Vector3_Order<int> &box2) const;
	template<typename Tmatrix>
	Tmatrix Ra2D_to_Km2D(const Parallel_Orbitals &pv, 
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &HR_a2D, const int ik) const;
	template<typename Tmatrix>
	void Ra2D_to_Km2D_mixing(const Parallel_Orbitals &pv, 
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &HR_a2D,
		std::vector<Tmatrix> &HK_m2D,
		std::vector<std::deque<Tmatrix>> &HK_m2D_pulay_seq) const;
	template<typename T>
	T pulay_mixing( const T &H_pulay_old, std::deque<T> &H_seq, const T &H_new ) const;

	#if EXX_H_COMM==1
	class Allreduce
	{
	public:
		Allreduce(
			const MPI_Comm &mpi_comm_in,
			std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_local_in );
		~Allreduce();
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> exx_to_a2D();

	private:
		void set_atom_in_2D(Parallel_Orbitals &pv);
		void ask( const int rank_delta_now );
		void recv_data_process( const int rank_data );
		void insert_data( std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> &data_rank );
		void insert_data();
		void send_data_process( const int rank_asked );
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> get_data_local_wrapper( 
			const std::pair<std::vector<bool>,std::vector<bool>> &atom_asked ) const;

	private:
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> data_all;
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_local;

		std::pair< std::vector<bool>, std::vector<bool> > atom_in_2D;

		const MPI_Comm & mpi_comm;
		int comm_sz;
		int my_rank;

		static constexpr int tag_ask  = 1;
		static constexpr int tag_data = 2;

		std::atomic_flag lock_insert;
		std::atomic<int> rank_delta;

		boost::mpi::packed_oarchive *oarp_atom_in_2D;
		std::vector<boost::mpi::packed_oarchive*> oarps_isend_data;
		std::vector<boost::mpi::packed_iarchive*> iarps_recv_data;
		std::vector<boost::mpi::packed_iarchive*> iarps_atom_asked;

		std::vector<std::atomic<int>*> flags_isend_data;
		std::vector<std::atomic<int>*> flags_ask_atom;
		boost::dynamic_bitset<> flags_recv_data;
	};
	#elif EXX_H_COMM==2
	class Allreduce2
	{
	public:
		void init(
			const MPI_Comm &mpi_comm_in,
            const std::set<std::pair<size_t, size_t>>& H_atom_pairs_core,
            const Parallel_Orbitals &pv);
        std::vector<std::map<size_t, std::map<size_t, std::map<Abfs::Vector3_Order<int>, ModuleBase::matrix>>>> exx_to_a2D(
			std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_local) const;

	private:
		enum class Flag_Send {undo, begin_oar, finish_oar, begin_isend, finish_isend};
		enum class Flag_Recv {undo, begin_irecv, begin_iar, finish_iar};
		
		//std::vector<std::pair<std::vector<bool>,std::vector<bool>>> get_atom_in_2D_list() const;
		std::vector<size_t> get_send_size_list(
			const std::set<std::pair<size_t,size_t>> &H_atom_pairs_core,
			const std::vector<std::pair<std::vector<bool>,std::vector<bool>>> &atom_in_2D_list) const;
		void init_flags(
			std::vector<std::atomic<Flag_Send>> &flags_send,
			std::vector<std::atomic<Flag_Recv>> &flags_recv) const;
		bool finish_judge(
			const std::vector<std::atomic<Flag_Send>> &flags_send,
			const std::vector<std::atomic<Flag_Recv>> &flags_recv) const;
		bool memory_enough(
			const int rank_send_next,
			const std::vector<std::atomic<Flag_Send>> &flags_send) const;
		void send_data_process(
			const int rank_send_now,
			const std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_local,
			std::vector<std::valarray<double>> &oars_isend,
			std::vector<std::atomic<Flag_Send>> &flags_send) const;
		void recv_data_process(
			const int rank_recv,
			std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_all,
			std::vector<std::valarray<double>> &iarps_irecv,
			std::vector<std::atomic<Flag_Recv>> &flags_recv,
			std::atomic_flag &lock_insert) const;			
		void insert_data(
			std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_rank,
			std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_all) const;

	private:
		MPI_Comm mpi_comm;
		int comm_sz;
		int my_rank;
		
		std::vector<std::pair<std::vector<bool>,std::vector<bool>>> atom_in_2D_list;
		std::vector<size_t> send_size_list;
		size_t recv_size;
	};
	#endif

public:
	std::vector<ModuleBase::matrix>        HK_Gamma_m2D;
	std::vector<ModuleBase::ComplexMatrix> HK_K_m2D;
	std::vector<std::deque<ModuleBase::matrix>>        HK_Gamma_m2D_pulay_seq;
	std::vector<std::deque<ModuleBase::ComplexMatrix>> HK_K_m2D_pulay_seq;

	enum class Mixing_Mode{ No, Plain, Pulay };
	Mixing_Mode mixing_mode;
	double mixing_beta;

#ifdef __MPI	
	Allreduce2 allreduce2;
#endif
};

#endif
#endif