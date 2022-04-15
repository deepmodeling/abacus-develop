#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_ALLREDUCE_TEMPLATE_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_ALLREDUCE_TEMPLATE_H

#include "exx_abfs-parallel-communicate-allreduce.h"

#include <vector>

#ifdef __MPI
#include "mpi.h"
#endif
#include <boost/dynamic_bitset.hpp>
#include "../src_lcao/serialization_boost.h"

#include <thread>
#include <atomic>

#ifdef __MKL
#include <mkl_service.h>
#endif

#include "../src_external/src_test/src_ri/exx_lcao-test.h"


template< typename T >
T Exx_Abfs::Parallel::Communicate::Allreduce::allreduce(
	T & data_local,
	std::function<void(T&,T&)> & insert_function)
{
std::ofstream ofs_mpi( "allreduce_"+ModuleBase::GlobalFunc::TO_STRING(my_rank), std::ofstream::app );
timeval t_start;
double t_probe=0, t_recv=0;
	
	std::vector<std::thread> threads;
	atomic_flag insert_lock = ATOMIC_FLAG_INIT;
	atomic<int> rank_delta=1;
	for( int i=0; i<2; ++i )
		ask(rank_delta++);

gettimeofday(&t_start,NULL);
	boost::mpi::packed_oarchive oar(mpi_comm);
	oar << data_local;
ofs_mpi<<"sizeof_oar:\t"<<get_sizeof(data_local)<<"\t"<<oar.size()<<std::endl;
ofs_mpi<<"oar<<\t"<<time_during(t_start)<<std::endl;

	std::vector<MPI_Request> request_isends(comm_sz);
	boost::dynamic_bitset<> flag_isends(comm_sz,false);
	flag_isends[my_rank] = true;

gettimeofday(&t_start,NULL);
	T data_all;
	insert_function( data_all, data_local );
	boost::dynamic_bitset<> flag_irecvs(comm_sz,false);
	flag_irecvs[my_rank] = true;
ofs_mpi<<"insert\t"<<time_during(t_start)<<std::endl;

	std::vector<boost::mpi::packed_iarchive*> iarps( comm_sz );
	for( auto &iarp : iarps )
		iarp = new boost::mpi::packed_iarchive{MPI_COMM_WORLD};

	#ifdef __MKL
    const int mkl_threads = mkl_get_max_threads();
	mkl_set_num_threads(1);
	#endif
	
	while( !flag_irecvs.all() || !flag_isends.all() )
	{
		MPI_Status status;
gettimeofday(&t_start,NULL);
		MPI_Probe( MPI_ANY_SOURCE, MPI_ANY_TAG, mpi_comm, &status );
t_probe+=time_during(t_start);
		switch(status.MPI_TAG)
		{
			case tag_ask:
			{
				const int rank_send = status.MPI_SOURCE;
				int message_ignore;
gettimeofday(&t_start,NULL);
				MPI_Recv( &message_ignore, 1, MPI_INT, rank_send, tag_ask, mpi_comm, MPI_STATUS_IGNORE );
t_recv+=time_during(t_start);
				MPI_Isend( oar.address(), oar.size(), MPI_PACKED, rank_send, tag_data, mpi_comm, &request_isends[rank_send] );
				flag_isends[rank_send] = true;
				break;
			}
			case tag_data:
			{
				const int rank_recv = status.MPI_SOURCE;
				
				int data_size;
				MPI_Get_count( &status, MPI_PACKED, &data_size );
				boost::mpi::packed_iarchive & iar = *iarps[rank_recv];
				iar.resize(data_size);

gettimeofday(&t_start,NULL);
				MPI_Recv( iar.address(), data_size, MPI_PACKED, rank_recv, tag_data, mpi_comm, MPI_STATUS_IGNORE );
t_recv+=time_during(t_start);

				threads.push_back(std::thread(
					&Exx_Abfs::Parallel::Communicate::Allreduce::allreduce_thread<T>, this,
					std::ref(iar), std::ref(insert_function), std::ref(data_all),
					std::ref(rank_delta), std::ref(insert_lock) ));
				flag_irecvs[rank_recv] = true;
ofs_mpi<<rank_recv<<"\t"<<threads.back().get_id()<<std::endl;
				break;
			}
			default:
				throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
for( int i=0; i<flag_isends.size(); ++i )
	ofs_mpi<<flag_isends[i];	ofs_mpi<<std::endl;
for( int i=0; i<flag_irecvs.size(); ++i )
	ofs_mpi<<flag_irecvs[i];	ofs_mpi<<std::endl;
	}

gettimeofday(&t_start,NULL);
	for( std::thread & t : threads )
		t.join();
	for( int i_rank=0; i_rank<comm_sz; ++i_rank )
		if( i_rank != my_rank )
			MPI_Wait( &request_isends[i_rank], MPI_STATUS_IGNORE );
	for( auto &iarp : iarps )
		delete iarp;
ofs_mpi<<"final\t"<<time_during(t_start)<<std::endl;

ofs_mpi<<"t_probe\t"<<t_probe<<std::endl;
ofs_mpi<<"t_recv\t"<<t_recv<<std::endl;

ofs_mpi.close();

	#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
	#endif
	
	return data_all;
}


template< typename T >
void Exx_Abfs::Parallel::Communicate::Allreduce::allreduce_thread(
	boost::mpi::packed_iarchive &iar,
	std::function<void(T&,T&)> &insert_function,
	T &data_all,
	atomic<int> &rank_delta,
	atomic_flag &insert_lock )
{
std::ofstream ofs_mpi( "allreduce_"+ModuleBase::GlobalFunc::TO_STRING(my_rank)+"_"+ModuleBase::GlobalFunc::TO_STRING(this_thread::get_id()), std::ofstream::app );
timeval t_start;

ofs_mpi<<"sizeof_iar:\t"<<iar.size()<<"\t";

gettimeofday(&t_start,NULL);
	T data_rank;
	iar >> data_rank;
ofs_mpi<<iar.size()<<"\t"<<get_sizeof(data_rank)<<std::endl;
ofs_mpi<<"iar>>\t"<<time_during(t_start)<<std::endl;

gettimeofday(&t_start,NULL);
	while( insert_lock.test_and_set() );
ofs_mpi<<"while\t"<<time_during(t_start)<<std::endl;
gettimeofday(&t_start,NULL);
	insert_function( data_all, data_rank );
	insert_lock.clear();
ofs_mpi<<"insert\t"<<time_during(t_start)<<std::endl;

	ask(rank_delta++);
ofs_mpi.close();
}


#endif