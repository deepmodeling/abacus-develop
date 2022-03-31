#include "exx_abfs-parallel-communicate-hexx.h"
#include "../src_lcao/serialization_boost.h"

#include "../src_pw/global.h"
#include "../src_lcao/global_fp.h"

#include <thread>

#ifdef __MKL
#include <mkl_service.h>
#endif

Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::Allreduce(
	const MPI_Comm & mpi_comm_in,
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_local_in )
	:mpi_comm(mpi_comm_in),
	 data_local(data_local_in),
	 lock_insert(ATOMIC_FLAG_INIT)
{
	data_all.resize(GlobalV::NSPIN);

	MPI_Comm_size( mpi_comm, &comm_sz );
	MPI_Comm_rank( mpi_comm, &my_rank );

	rank_delta = 1;

	set_atom_in_2D();
	oarp_atom_in_2D = new boost::mpi::packed_oarchive{mpi_comm};
	*oarp_atom_in_2D << atom_in_2D;

	oarps_isend_data.resize( comm_sz );
	for( auto &oarp_isend_data : oarps_isend_data )
		oarp_isend_data = new boost::mpi::packed_oarchive{mpi_comm};
	iarps_recv_data.resize( comm_sz );
	for( auto &iarp_recv_data : iarps_recv_data )
		iarp_recv_data = new boost::mpi::packed_iarchive{mpi_comm};
	iarps_atom_asked.resize( comm_sz );
	for( auto &iarp_atom_asked : iarps_atom_asked )
		iarp_atom_asked = new boost::mpi::packed_iarchive{mpi_comm};

	// 0: undo		1: thread finish		2: MPI begin
	flags_isend_data.resize(comm_sz);
	for( int irank=0; irank!=comm_sz; ++irank )
		flags_isend_data[irank] = new atomic<int>(0);
	*flags_isend_data[my_rank] = 2;

	flags_ask_atom.resize(comm_sz);
	for( int irank=0; irank!=comm_sz; ++irank )
		flags_ask_atom[irank] = new atomic<int>(0);
	*flags_ask_atom[my_rank] = 2;

	flags_recv_data.resize(comm_sz,false);
	flags_recv_data[my_rank] = true;
}


Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::~Allreduce()
{
	delete oarp_atom_in_2D;
	for( auto &oarp_isend_data : oarps_isend_data )
		if(oarp_isend_data)
			delete oarp_isend_data;
	for( auto &iarp_recv_data : iarps_recv_data )
		if(iarp_recv_data)
			delete iarp_recv_data;
	for( auto &iarp_atom_asked : iarps_atom_asked )
		if(iarp_atom_asked)
			delete iarp_atom_asked;
	for( int irank=0; irank!=comm_sz; ++irank )
		delete flags_isend_data[irank];
	for( int irank=0; irank!=comm_sz; ++irank )
		delete flags_ask_atom[irank];
}



std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>>
Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::exx_to_a2D()
{
	auto clear_oars = [&]( std::vector<MPI_Request> &requests, boost::dynamic_bitset<> &flags, std::vector<boost::mpi::packed_oarchive*> &oarps, const std::string &s )
	{
		if( flags.none() )	return;
		for( int irank=0; irank!=comm_sz; ++irank )
		{
			if( flags[irank] )
			{
				int flag_finish;
				if(MPI_SUCCESS!=MPI_Test( &requests[irank], &flag_finish, MPI_STATUS_IGNORE ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				if( flag_finish )
				{
					delete oarps[irank];
					oarps[irank] = nullptr;
				}
				flags[irank] = false;
			}
		}
	};

	auto if_finish = []( const std::vector<atomic<int>*> &flags ) -> bool
	{
		int sum=0;
		for( size_t i=0; i<flags.size(); ++i )
			sum += *flags[i];
		return sum == flags.size()*2;
	};

	for( int i=0; i<2; ++i )
		ask(rank_delta++);

	std::vector<std::thread> threads;

	std::vector<MPI_Request> requests_isend_data(comm_sz);
	std::vector<MPI_Request> requests_ask(comm_sz);
	boost::dynamic_bitset<> flags_request_isend_data(comm_sz,false);

	#ifdef __MKL
    const int mkl_threads = mkl_get_max_threads();
	mkl_set_num_threads(1);
	#endif

	while( !if_finish(flags_ask_atom) || !if_finish(flags_isend_data) || !flags_recv_data.all() )
	{
		MPI_Status status;
		int flag_message;
		if(MPI_SUCCESS!=MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, mpi_comm, &flag_message, &status ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		if(flag_message)
		{
			int message_size;
			if(MPI_SUCCESS!=MPI_Get_count( &status, MPI_PACKED, &message_size ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));

			switch(status.MPI_TAG)
			{
				case tag_ask:
				{
					const int rank_asked = status.MPI_SOURCE;

					iarps_atom_asked[rank_asked]->resize(message_size);
					if(MPI_SUCCESS!=MPI_Recv( iarps_atom_asked[rank_asked]->address(), message_size, MPI_PACKED, rank_asked, tag_ask, mpi_comm, MPI_STATUS_IGNORE ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));

					threads.push_back(std::thread(
						&Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::send_data_process, this,
						rank_asked ));
					break;
				}
				case tag_data:
				{
					const int rank_data = status.MPI_SOURCE;

					iarps_recv_data[rank_data]->resize(message_size);
					if(MPI_SUCCESS!=MPI_Recv( iarps_recv_data[rank_data]->address(), message_size, MPI_PACKED, rank_data, tag_data, mpi_comm, MPI_STATUS_IGNORE ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
					flags_recv_data[rank_data] = true;

					threads.push_back(std::thread(
						&Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::recv_data_process, this,
						rank_data ));
					break;
				}
				default:
					throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
			}
		}

		for( int rank_ask=0; rank_ask!=comm_sz; ++rank_ask )
			if( *flags_ask_atom[rank_ask] == 1 )
			{
				if(MPI_SUCCESS!=MPI_Isend( oarp_atom_in_2D->address(), oarp_atom_in_2D->size(), MPI_PACKED, rank_ask, tag_ask, mpi_comm, &requests_ask[rank_ask] ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				*flags_ask_atom[rank_ask] = 2;
			}
		for( int rank_asked=0; rank_asked!=comm_sz; ++rank_asked )
			if( *flags_isend_data[rank_asked] == 1 )
			{
				if(MPI_SUCCESS!=MPI_Isend( oarps_isend_data[rank_asked]->address(), oarps_isend_data[rank_asked]->size(), MPI_PACKED, rank_asked, tag_data, mpi_comm, &requests_isend_data[rank_asked] ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				flags_request_isend_data[rank_asked] = true;
				*flags_isend_data[rank_asked] = 2;
			}

		clear_oars( requests_isend_data, flags_request_isend_data, oarps_isend_data, "oarps_isend_data" );
	}

	while( lock_insert.test_and_set() );
	insert_data();
	lock_insert.clear();

	for( std::thread & t : threads )
		t.join();
	for( int i_rank=0; i_rank<comm_sz; ++i_rank )
		if( i_rank != my_rank )
		{
			if(MPI_SUCCESS!=MPI_Wait( &requests_isend_data[i_rank], MPI_STATUS_IGNORE ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
			if(MPI_SUCCESS!=MPI_Wait( &requests_ask[i_rank], MPI_STATUS_IGNORE ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

	#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
	#endif

	return data_all;
}



void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::set_atom_in_2D(Parallel_Orbitals &pv)
{
	atom_in_2D.first.resize(GlobalC::ucell.nat,false);
	atom_in_2D.second.resize(GlobalC::ucell.nat,false);
	for( size_t it=0; it!=GlobalC::ucell.ntype; ++it )
	{
		for( size_t ia=0; ia!=GlobalC::ucell.atoms[it].na; ++ia )
		{
			const size_t iat = GlobalC::ucell.itia2iat(it,ia);
			for( size_t iw=0; iw!=GlobalC::ucell.atoms[it].nw; ++iw )
			{
				const size_t iwt = GlobalC::ucell.itiaiw2iwt(it,ia,iw);
				if( pv.trace_loc_row[iwt]>=0 )
					atom_in_2D.first[iat] = true;
				if( pv.trace_loc_col[iwt]>=0 )
					atom_in_2D.second[iat] = true;
			}
		}
	}
}


void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::ask( const int rank_delta_now )
{
	if( rank_delta_now < comm_sz )
	{
		const int rank_ask = ( my_rank + rank_delta_now ) % comm_sz;
		*flags_ask_atom[rank_ask] = 1;
	}
}



void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::recv_data_process( const int rank_data )
{
	auto vector_empty = []( const std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> & v ) -> bool
	{
		for( const auto &i : v )
			if(!i.empty())	return false;
		return true;
	};

	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> data_rank;
	*iarps_recv_data[rank_data] >> data_rank;
	iarps_recv_data[rank_data]->resize(0);

	if(!vector_empty(data_rank))
	{
		while( lock_insert.test_and_set() );
		insert_data(data_rank);
		lock_insert.clear();
	}

	ask(rank_delta++);
	
	#ifdef MATRIX_WRAPPER_TIANHE2
	for( auto &data_rank_is : data_rank )
		for( auto &data_rank_A : data_rank_is )
			for( auto &data_rank_B : data_rank_A.second )
				for( auto &data_rank_C : data_rank_B.second )
					if(!data_rank_C.second.c)
						delete[] data_rank_C.second.c;
	#endif
}


void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::insert_data(
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> &data_rank )
{
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		auto &data_rank_is = data_rank[is];
		auto &data_all_is = data_all[is];
		for( auto &data_rank_A : data_rank_is )
		{
			const size_t iat1 = data_rank_A.first;
			auto &data_all_A = data_all_is[iat1];
			for( auto &data_rank_B : data_rank_A.second )
			{
				const size_t iat2 = data_rank_B.first;
				auto &data_all_B = data_all_A[iat2];
				for( auto &data_rank_C : data_rank_B.second )
				{
					const Abfs::Vector3_Order<int> &box2 = data_rank_C.first;
					auto &data_all_C = data_all_B[box2];
					if(data_all_C.c)
						data_all_C += data_rank_C.second.to_matrix();
					else
						data_all_C = data_rank_C.second.to_matrix();
				}
			}
		}
	}
}


void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::insert_data()
{
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		auto &data_local_is = data_local[is];
		auto &data_all_is = data_all[is];
		for( auto &data_local_A : data_local_is )
		{
			const size_t iat1 = data_local_A.first;
			if( !atom_in_2D.first[iat1] )	continue;
			auto &data_all_A = data_all_is[iat1];
			for( auto &data_local_B : data_local_A.second )
			{
				const size_t iat2 = data_local_B.first;
				if( !atom_in_2D.second[iat2] )	continue;
				auto &data_all_B = data_all_A[iat2];
				for( auto &data_local_C : data_local_B.second )
				{
					const Abfs::Vector3_Order<int> &box2 = data_local_C.first;
					auto &data_all_C = data_all_B[box2];
					if( data_all_C.c )
						data_all_C += data_local_C.second;
					else
						data_all_C = std::move(data_local_C.second);
				}
			}
		}
	}
}



void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::send_data_process( const int rank_asked )
{
	std::pair< std::vector<bool>, std::vector<bool> > atom_asked;
	*iarps_atom_asked[rank_asked] >> atom_asked;
	iarps_atom_asked[rank_asked]->resize(0);

	const std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> matrix_wrapped = get_data_local_wrapper(atom_asked);
	*oarps_isend_data[rank_asked] << matrix_wrapped;

	*flags_isend_data[rank_asked] = 1;
}


std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>>
Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::get_data_local_wrapper( const std::pair<std::vector<bool>,std::vector<bool>> &atom_asked ) const
{
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> mw(GlobalV::NSPIN);
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		auto &data_local_is = data_local[is];
		auto &mw_is = mw[is];
		for( auto &data_local_A : data_local_is )
		{
			const size_t iat1 = data_local_A.first;
			if( !atom_asked.first[iat1] )	continue;
			auto &mw_A = mw_is[iat1];
			for( auto &data_local_B : data_local_A.second )
			{
				const size_t iat2 = data_local_B.first;
				if( !atom_asked.second[iat2] )	continue;
				auto &mw_B = mw_A[iat2];
				for( auto &data_local_C : data_local_B.second )
				{
					const Abfs::Vector3_Order<int> &box2 = data_local_C.first;
					mw_B[box2] = data_local_C.second;
				}
			}
		}
	}
	return mw;
}
