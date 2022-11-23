#include "exx_abfs-parallel-communicate-dm.h"
#include "../src_lcao/serialization_boost.h"

#include <thread>

#ifdef __MKL
#include <mkl_service.h>
#endif

Exx_Abfs::Parallel::Communicate::DM::Allreduce::Allreduce( 
	const MPI_Comm & mpi_comm_in, 
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_local_in,
	const Abfs::Vector3_Order<int> &Born_von_Karman_period,
	const set<std::pair<size_t,size_t>> &H_atom_pairs_core)
	:mpi_comm(mpi_comm_in),
	 data_local(data_local_in),
	 lock_insert(ATOMIC_FLAG_INIT),
	 lock_atom_unset_read(0)
{
	data_all.resize(GlobalV::NSPIN);
	
	set<Abfs::Vector3_Order<int>> Born_Von_Karman_boxes;
	for( int ix=0; ix!=Born_von_Karman_period.x; ++ix )
		for( int iy=0; iy!=Born_von_Karman_period.y; ++iy )
			for( int iz=0; iz!=Born_von_Karman_period.z; ++iz )
				Born_Von_Karman_boxes.insert({ix,iy,iz});
	for( const auto std::pair : H_atom_pairs_core )
		atom_unset[std::pair.first][std::pair.second] = Born_Von_Karman_boxes;

	MPI_Comm_size( mpi_comm, &comm_sz );
	MPI_Comm_rank( mpi_comm, &my_rank );
	
	rank_delta = 1;
	
	oarps_isend_data.resize( comm_sz );
	for( auto &oarp_isend_data : oarps_isend_data )
		oarp_isend_data = new boost::mpi::packed_oarchive{mpi_comm};
	oarps_atom_unset.resize( comm_sz );
	for( auto &oarp_atom_unset : oarps_atom_unset )
		oarp_atom_unset = new boost::mpi::packed_oarchive{mpi_comm};
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


Exx_Abfs::Parallel::Communicate::DM::Allreduce::~Allreduce()
{
	for( auto &oarp_isend_data : oarps_isend_data )
		if(oarp_isend_data)
			delete oarp_isend_data;
	for( auto &oarp_atom_unset : oarps_atom_unset )
		if(oarp_atom_unset)
			delete oarp_atom_unset;
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
Exx_Abfs::Parallel::Communicate::DM::Allreduce::grid_to_exx()
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
	
	data_localw = get_data_local_wrapper();
		
	do{ while( lock_atom_unset_read ); } while( lock_insert.test_and_set() );
	insert_data( data_local );
	lock_insert.clear();
	
	if( atom_unset.empty() )
		while(rank_delta<comm_sz)
			ask(rank_delta++);
	else
		for( int i=0; i<2; ++i )
			ask(rank_delta++);
	
	std::vector<std::thread> threads;

	std::vector<MPI_Request> requests_isend_data(comm_sz);
	std::vector<MPI_Request> requests_ask(comm_sz);
	boost::dynamic_bitset<> flags_request_isend_data(comm_sz,false);
	boost::dynamic_bitset<> flags_request_ask(comm_sz,false);
		
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
						&Exx_Abfs::Parallel::Communicate::DM::Allreduce::send_data_process, this,
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
						&Exx_Abfs::Parallel::Communicate::DM::Allreduce::recv_data_process, this,
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
				if(MPI_SUCCESS!=MPI_Isend( oarps_atom_unset[rank_ask]->address(), oarps_atom_unset[rank_ask]->size(), MPI_PACKED, rank_ask, tag_ask, mpi_comm, &requests_ask[rank_ask] ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				flags_request_ask[rank_ask] = true;
				*flags_ask_atom[rank_ask] = 2;
			}
		for( int rank_asked=0; rank_asked!=comm_sz; ++rank_asked )
			if( *flags_isend_data[rank_asked] == 1 )
			{
				if(MPI_SUCCESS!=MPI_Isend( oarps_isend_data[rank_asked]->address(), oarps_isend_data[rank_asked]->size(), MPI_PACKED, rank_asked, tag_data, mpi_comm, &requests_isend_data[rank_asked] ))	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				flags_request_isend_data[rank_asked] = true;
				*flags_isend_data[rank_asked] = 2;
			}
		
		clear_oars( requests_ask, flags_request_ask, oarps_atom_unset, "oarps_atom_unset" );
		clear_oars( requests_isend_data, flags_request_isend_data, oarps_isend_data, "oarps_isend_data" );
//test_flags(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	}

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


void Exx_Abfs::Parallel::Communicate::DM::Allreduce::ask( const int rank_delta_now )
{
	if( rank_delta_now < comm_sz )
	{
		const int rank_ask = ( my_rank + rank_delta_now ) % comm_sz;
		while( lock_insert.test_and_set() );
		++lock_atom_unset_read;
		lock_insert.clear();
		*oarps_atom_unset[rank_ask] << atom_unset;
		--lock_atom_unset_read;
		*flags_ask_atom[rank_ask] = 1;
	}
}



void Exx_Abfs::Parallel::Communicate::DM::Allreduce::recv_data_process( const int rank_data )
{
	auto vector_empty = []( const std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> & v ) -> bool
	{
		for( const auto &i : v )
			if(!i.empty())	return false;
		return true;
	};

	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> data_rank;
	*iarps_recv_data[rank_data] >> data_rank;
//	delete iarps_recv_data[rank_data];		iarps_recv_data[rank_data]=nullptr;
	iarps_recv_data[rank_data]->resize(0);
	
	if(!vector_empty(data_rank))
	{
		do{ while( lock_atom_unset_read ); } while( lock_insert.test_and_set() );
		insert_data(data_rank);
		lock_insert.clear();
	}
	
	if( atom_unset.empty() )
		while(rank_delta<comm_sz)
			ask(rank_delta++);
	else
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


void Exx_Abfs::Parallel::Communicate::DM::Allreduce::insert_data( 
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_rank )
{
	std::vector<const std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>*> data_rank_Ap(GlobalV::NSPIN,nullptr);
	std::vector<const std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>*> data_rank_Bp(GlobalV::NSPIN,nullptr);
	std::vector<const ModuleBase::matrix*> data_rank_Cp(GlobalV::NSPIN,nullptr);
	
	for( auto atom_unset_Ap=atom_unset.begin(); atom_unset_Ap!=atom_unset.end(); )
	{
		const size_t iat1 = atom_unset_Ap->first;
		for( int is=0; is!=GlobalV::NSPIN; ++is )	data_rank_Ap[is] = static_cast<const std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>*>( ModuleBase::GlobalFunc::MAP_EXIST( data_rank[is], iat1 ) );
		if( !data_rank_Ap[0] ){ ++atom_unset_Ap; continue; }
		
		for( auto atom_unset_Bp=atom_unset_Ap->second.begin(); atom_unset_Bp!=atom_unset_Ap->second.end(); )
		{
			const size_t iat2 = atom_unset_Bp->first;
			for( int is=0; is!=GlobalV::NSPIN; ++is )	data_rank_Bp[is] = static_cast<const std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>*>( ModuleBase::GlobalFunc::MAP_EXIST( *data_rank_Ap[is], iat2 ) );
			if( !data_rank_Bp[0] ){ ++atom_unset_Bp; continue; }
		
			for( auto atom_unset_Cp=atom_unset_Bp->second.begin(); atom_unset_Cp!=atom_unset_Bp->second.end(); )
			{
				const Abfs::Vector3_Order<int> &box2 = *atom_unset_Cp;
				for( int is=0; is!=GlobalV::NSPIN; ++is )	data_rank_Cp[is] = static_cast<const ModuleBase::matrix*>( ModuleBase::GlobalFunc::MAP_EXIST( *data_rank_Bp[is], box2 ) );
				if( !data_rank_Cp[0] ){ ++atom_unset_Cp; continue; }
				
				for( int is=0; is!=GlobalV::NSPIN; ++is )
					if( data_rank_Cp[is]->c )
					{
						auto &data_all_C = data_all[is][iat1][iat2][box2];
						if(!data_all_C.c)
							data_all_C = std::move(*data_rank_Cp[is]);
					}
				
				atom_unset_Bp->second.erase( atom_unset_Cp++ );
			}
			if( atom_unset_Bp->second.empty() )
				atom_unset_Ap->second.erase( atom_unset_Bp++ );
			else
				++atom_unset_Bp;
		}
		if( atom_unset_Ap->second.empty() )
			atom_unset.erase( atom_unset_Ap++ );
		else
			++atom_unset_Ap;
	}
}

void Exx_Abfs::Parallel::Communicate::DM::Allreduce::insert_data( 
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> &data_rank )
{
	std::vector<std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>*> data_rank_Ap(GlobalV::NSPIN,nullptr);
	std::vector<std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>*> data_rank_Bp(GlobalV::NSPIN,nullptr);
	std::vector<ModuleBase::Matrix_Wrapper*> data_rank_Cp(GlobalV::NSPIN,nullptr);
	
	for( auto atom_unset_Ap=atom_unset.begin(); atom_unset_Ap!=atom_unset.end(); )
	{
		const size_t iat1 = atom_unset_Ap->first;
		for( int is=0; is!=GlobalV::NSPIN; ++is )	data_rank_Ap[is] = static_cast<std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>*>(ModuleBase::GlobalFunc::MAP_EXIST( data_rank[is], iat1 ));
		if( !data_rank_Ap[0] ){ ++atom_unset_Ap; continue; }
		
		for( auto atom_unset_Bp=atom_unset_Ap->second.begin(); atom_unset_Bp!=atom_unset_Ap->second.end(); )
		{
			const size_t iat2 = atom_unset_Bp->first;
			for( int is=0; is!=GlobalV::NSPIN; ++is )	data_rank_Bp[is] = static_cast<std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>*>(ModuleBase::GlobalFunc::MAP_EXIST( *data_rank_Ap[is], iat2 ));
			if( !data_rank_Bp[0] ){ ++atom_unset_Bp; continue; }
		
			for( auto atom_unset_Cp=atom_unset_Bp->second.begin(); atom_unset_Cp!=atom_unset_Bp->second.end(); )
			{
				const Abfs::Vector3_Order<int> &box2 = *atom_unset_Cp;
				for( int is=0; is!=GlobalV::NSPIN; ++is )	data_rank_Cp[is] = static_cast<ModuleBase::Matrix_Wrapper*>(ModuleBase::GlobalFunc::MAP_EXIST( *data_rank_Bp[is], box2 ));
				if( !data_rank_Cp[0] ){ ++atom_unset_Cp; continue; }
				
				for( int is=0; is!=GlobalV::NSPIN; ++is )
					if( data_rank_Cp[is]->c )
					{
						auto & data_all_C = data_all[is][iat1][iat2][box2];
						if( !data_all_C.c )
							data_all_C = data_rank_Cp[is]->to_matrix();
					}
				
				atom_unset_Bp->second.erase( atom_unset_Cp++ );
			}
			if( atom_unset_Bp->second.empty() )
				atom_unset_Ap->second.erase( atom_unset_Bp++ );
			else
				++atom_unset_Bp;
		}
		if( atom_unset_Ap->second.empty() )
			atom_unset.erase( atom_unset_Ap++ );
		else
			++atom_unset_Ap;
	}
}



void Exx_Abfs::Parallel::Communicate::DM::Allreduce::send_data_process( const int rank_asked )
{
	
	std::map<size_t,std::map<size_t,set<Abfs::Vector3_Order<int>>>> atom_asked;
	*iarps_atom_asked[rank_asked] >> atom_asked;
//	delete iarps_atom_asked[rank_asked];		iarps_atom_asked[rank_asked]=nullptr;
	iarps_atom_asked[rank_asked]->resize(0);
	
	const std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> matrix_wrapped = get_data_local_wrapper(atom_asked);
	*oarps_isend_data[rank_asked] << matrix_wrapped;	
	
	*flags_isend_data[rank_asked] = 1;
}

/*
std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>>
Exx_Abfs::Parallel::Communicate::DM::Allreduce::get_data_local_wrapper( 
	const std::map<size_t,std::map<size_t,set<Abfs::Vector3_Order<int>>>> & atom_asked ) const
{

	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> mw(GlobalV::NSPIN);
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		auto &mw_is = mw[is];
		const auto &data_local_is = data_local[is];
		for( const auto & atom_asked_A : atom_asked )
		{
			const size_t iat1 = atom_asked_A.first;
			if( auto data_local_A = static_cast<const std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>> * const>(ModuleBase::GlobalFunc::MAP_EXIST( data_local_is, iat1 )) )
			{
				auto &mw_A = mw_is[iat1];
				for( const auto & atom_asked_B : atom_asked_A.second )
				{
					const size_t iat2 = atom_asked_B.first;
					if( auto data_local_B = static_cast<const std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix> * const>(ModuleBase::GlobalFunc::MAP_EXIST( *data_local_A, iat2 )) )
					{
						auto &mw_B = mw_A[iat2];
						for( const auto & atom_asked_C : atom_asked_B.second )
						{
							const Abfs::Vector3_Order<int> &box2 = atom_asked_C;
							if( auto data_local_C = static_cast<const ModuleBase::matrix * const>(ModuleBase::GlobalFunc::MAP_EXIST( *data_local_B, box2 )) )
							{
								mw_B[box2] = *data_local_C;
							}
						}
					}
				}
			}
		}
	}
	return mw;
}
*/

std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>>
Exx_Abfs::Parallel::Communicate::DM::Allreduce::get_data_local_wrapper( 
	const std::map<size_t,std::map<size_t,set<Abfs::Vector3_Order<int>>>> & atom_asked ) const
{

	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> mw(GlobalV::NSPIN);
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		auto &mw_is = mw[is];
		const auto &data_localw_is = data_localw[is];
		for( const auto & atom_asked_A : atom_asked )
		{
			const size_t iat1 = atom_asked_A.first;
			if( auto data_localw_A = static_cast<const std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>> * const>(ModuleBase::GlobalFunc::MAP_EXIST( data_localw_is, iat1 )) )
			{
				auto &mw_A = mw_is[iat1];
				for( const auto & atom_asked_B : atom_asked_A.second )
				{
					const size_t iat2 = atom_asked_B.first;
					if( auto data_localw_B = static_cast<const std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper> * const>(ModuleBase::GlobalFunc::MAP_EXIST( *data_localw_A, iat2 )) )
					{
						auto &mw_B = mw_A[iat2];
						for( const auto & atom_asked_C : atom_asked_B.second )
						{
							const Abfs::Vector3_Order<int> &box2 = atom_asked_C;
							if( auto data_localw_C = static_cast<const ModuleBase::Matrix_Wrapper * const>(ModuleBase::GlobalFunc::MAP_EXIST( *data_localw_B, box2 )) )
							{
								mw_B[box2] = *data_localw_C;
							}
						}
					}
				}
			}
		}
	}
	return mw;
}

std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>>
Exx_Abfs::Parallel::Communicate::DM::Allreduce::get_data_local_wrapper() const
{

	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::Matrix_Wrapper>>>> mw(GlobalV::NSPIN);
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		const auto &data_local_is = data_local[is];
		auto &mw_is = mw[is];
		for( const auto &data_local_A : data_local_is )
		{
			const size_t iat1 = data_local_A.first;
			auto &mw_A = mw_is[iat1];
			for( const auto &data_local_B : data_local_A.second )
			{
				const size_t iat2 = data_local_B.first;
				auto &mw_B = mw_A[iat2];
				for( const auto &data_local_C : data_local_B.second )
				{
					const Abfs::Vector3_Order<int> &box2 = data_local_C.first;
					mw_B[box2] = data_local_C.second;
				}
			}
		}
	}
	return mw;
}
