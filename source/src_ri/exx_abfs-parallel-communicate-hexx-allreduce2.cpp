#ifdef __MPI
#include "exx_abfs-parallel-communicate-hexx.h"
#include "exx_abfs-parallel-communicate-function.h"

#include "../src_pw/global.h"
#include "../src_lcao/global_fp.h"

#include "mpi.h"
#include <thread>

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::init(
	const MPI_Comm &mpi_comm_in,
    const set<std::pair<size_t, size_t>>& H_atom_pairs_core,
    const Parallel_Orbitals &pv)
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::init");
	
	mpi_comm = mpi_comm_in;
	MPI_Comm_size(mpi_comm, &comm_sz);
	MPI_Comm_rank(mpi_comm, &my_rank);
	
	atom_in_2D_list = Exx_Abfs::Parallel::Communicate::Function::get_atom_in_2D_list(mpi_comm, pv);
	send_size_list = get_send_size_list(H_atom_pairs_core, atom_in_2D_list);
	recv_size = pv.nrow * pv.ncol * sizeof(double);
}

std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> 
Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::exx_to_a2D(
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_local) const
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::exx_to_a2D");
	
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> data_all(GlobalV::NSPIN);

	std::vector<atomic<Flag_Send>> flags_send(comm_sz);
	std::vector<atomic<Flag_Recv>> flags_recv(comm_sz);
	init_flags(flags_send, flags_recv);

	std::vector<std::valarray<double>> oarps_isend(comm_sz);
	std::vector<std::valarray<double>> iarps_irecv(comm_sz);

	int rank_send_now = my_rank;
	auto rank_send_next = [&]()->int{ return (rank_send_now+1)%comm_sz; };

	atomic_flag lock_insert = ATOMIC_FLAG_INIT;
	std::vector<thread> threads;
	std::vector<MPI_Request>requests_isend(comm_sz);
	std::vector<MPI_Request>requests_irecv(comm_sz);
	
	while(!finish_judge(flags_send, flags_recv))
	{
		if(rank_send_next()!=my_rank && memory_enough(rank_send_next(), flags_send))
		{
			rank_send_now = rank_send_next();
			flags_send[rank_send_now] = Flag_Send::begin_oar;
			threads.push_back(std::thread(
				&Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::send_data_process, this,
				rank_send_now, std::cref(data_local), std::ref(oarps_isend), std::ref(flags_send) ));
		}
		for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
		{
			if( flags_send[rank_send] == Flag_Send::finish_oar )
			{
				if(MPI_Isend( ModuleBase::GlobalFunc::VECTOR_TO_PTR(oarps_isend[rank_send]), oarps_isend[rank_send].size(), MPI_DOUBLE, rank_send, 0, mpi_comm, &requests_isend[rank_send] )!=MPI_SUCCESS)	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				flags_send[rank_send] = Flag_Send::begin_isend;
			}
		}
		for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
		{
			if( flags_send[rank_send] == Flag_Send::begin_isend )
			{
				int flag_finish;
				if(MPI_Test( &requests_isend[rank_send], &flag_finish, MPI_STATUS_IGNORE )!=MPI_SUCCESS)	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				if(flag_finish)
				{
					oarps_isend[rank_send].resize(0);
					flags_send[rank_send] = Flag_Send::finish_isend;
				}
			}
		}

		{
			MPI_Status status;
			int flag_message;
			if(MPI_Iprobe( MPI_ANY_SOURCE, 0, mpi_comm, &flag_message, &status )!=MPI_SUCCESS)	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
			if(flag_message)
			{
				int message_size;
				if(MPI_Get_count( &status, MPI_PACKED, &message_size )!=MPI_SUCCESS)	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				const int rank_recv = status.MPI_SOURCE;
				iarps_irecv[rank_recv].resize(message_size);
				if(MPI_Irecv( ModuleBase::GlobalFunc::VECTOR_TO_PTR(iarps_irecv[rank_recv]), message_size, MPI_DOUBLE, rank_recv, 0, mpi_comm, &requests_irecv[rank_recv] )!=MPI_SUCCESS)	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				flags_recv[rank_recv] = Flag_Recv::begin_irecv;
			}
		}
		for(int rank_recv=0; rank_recv!=comm_sz; ++rank_recv)
		{
			if(flags_recv[rank_recv] == Flag_Recv::begin_irecv)
			{
				int flag_finish;
				if(MPI_Test( &requests_irecv[rank_recv], &flag_finish, MPI_STATUS_IGNORE )!=MPI_SUCCESS)	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
				if(flag_finish)
				{
					flags_recv[rank_recv] = Flag_Recv::begin_iar;
					threads.push_back(std::thread(
						&Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::recv_data_process, this,
						rank_recv, std::ref(data_all), std::ref(iarps_irecv), std::ref(flags_recv), std::ref(lock_insert) ));
				}
			}
		}
	}

	while( lock_insert.test_and_set() );
	insert_data(data_local, data_all);
	lock_insert.clear();
	
	for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
	{
		if( flags_send[rank_send] == Flag_Send::begin_isend )
		{
			if(MPI_Wait( &requests_isend[rank_send], MPI_STATUS_IGNORE )!=MPI_SUCCESS)	throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
			oarps_isend[rank_send].resize(0);
			flags_send[rank_send] = Flag_Send::finish_isend;
		}
	}
	
	for(std::thread &t : threads)
		t.join();	
	return data_all;
}

// the upper limit of size sended to each process
std::vector<size_t> Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::get_send_size_list(
	const set<std::pair<size_t,size_t>> &H_atom_pairs_core,
	const std::vector<std::pair<std::vector<bool>,std::vector<bool>>> &atom_in_2D_list) const
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::get_send_size_list");
	std::vector<size_t> send_size_list(comm_sz,0);
	for(const auto &H_atom_pair : H_atom_pairs_core)
	{
		const size_t iat1 = H_atom_pair.first;
		const size_t iat2 = H_atom_pair.second;
		for(int rank_send=0; rank_send<comm_sz; ++rank_send)
		{
			if(atom_in_2D_list[rank_send].first[iat1] && atom_in_2D_list[rank_send].second[iat2])
				send_size_list[rank_send] += GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw * GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw * sizeof(double);
		}
	}
	return send_size_list;
}

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::init_flags(
	std::vector<atomic<Flag_Send>> &flags_send,
	std::vector<atomic<Flag_Recv>> &flags_recv) const
{
	for(atomic<Flag_Send> &flag_send : flags_send)
		flag_send = Flag_Send::undo;
	flags_send[my_rank] = Flag_Send::finish_isend;
	
	for(atomic<Flag_Recv> &flag_recv : flags_recv)
		flag_recv = Flag_Recv::undo;
	flags_recv[my_rank] = Flag_Recv::finish_iar;
}

bool Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::finish_judge(
	const std::vector<atomic<Flag_Send>> &flags_send,
	const std::vector<atomic<Flag_Recv>> &flags_recv) const
{
	for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
		if((flags_send[rank_send]!=Flag_Send::begin_isend) && (flags_send[rank_send]!=Flag_Send::finish_isend))
			return false;
	for(int rank_recv=0; rank_recv!=comm_sz; ++rank_recv)
		if((flags_recv[rank_recv]!=Flag_Recv::begin_iar) && (flags_recv[rank_recv]!=Flag_Recv::finish_iar))
			return false;
	return true;
}

bool Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::memory_enough(
	const int rank_send_next,
	const std::vector<atomic<Flag_Send>> &flags_send) const
{
	size_t memory_need = recv_size;
	for(int rank_send=0; rank_send<comm_sz; ++rank_send)
		if(flags_send[rank_send]==Flag_Send::begin_oar)
			memory_need += send_size_list[rank_send];
	const size_t memory_available = ModuleBase::GlobalFunc::MemAvailable()*1024;
	return (memory_available-memory_need>send_size_list[rank_send_next]) ? true : false;
}

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::send_data_process(
	const int rank_send_now,
	const std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_local,
	std::vector<std::valarray<double>> &oarps_isend,
	std::vector<atomic<Flag_Send>> &flags_send) const
{
	std::valarray<size_t> send_size(GlobalV::NSPIN);
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,const ModuleBase::matrix*>>>> mw(GlobalV::NSPIN);
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		send_size[is]=0;
		auto &data_local_is = data_local[is];
		auto &mw_is = mw[is];
		for( auto &data_local_A : data_local_is )
		{
			const size_t iat1 = data_local_A.first;
			if( !atom_in_2D_list[rank_send_now].first[iat1] )	continue;
			auto &mw_A = mw_is[iat1];
			for( auto &data_local_B : data_local_A.second )
			{
				const size_t iat2 = data_local_B.first;
				if( !atom_in_2D_list[rank_send_now].second[iat2] )	continue;
				auto &mw_B = mw_A[iat2];
				for( auto &data_local_C : data_local_B.second )
				{
					const Abfs::Vector3_Order<int> &box2 = data_local_C.first;
					mw_B[box2] = &data_local_C.second;
					send_size[is] += 7 + data_local_C.second.nr * data_local_C.second.nc;
				}
			}
		}
	}
	
	oarps_isend[rank_send_now].resize(send_size.sum()+GlobalV::NSPIN);
	double * ptr = ModuleBase::GlobalFunc::VECTOR_TO_PTR(oarps_isend[rank_send_now]);
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		ptr[0] = send_size[is];
		ptr += 1;
		for( const auto &mw_A : mw[is] )
		{
			for( const auto &mw_B : mw_A.second )
			{
				for( const auto &mw_C : mw_B.second )
				{
					ptr[0] = mw_A.first;
					ptr[1] = mw_B.first;
					ptr[2] = mw_C.first.x;    ptr[3] = mw_C.first.y;    ptr[4] = mw_C.first.z;
					const ModuleBase::matrix &m = *mw_C.second;
					ptr[5] = m.nr;   ptr[6] = m.nc;
					memcpy( ptr+7, m.c, m.nr*m.nc*sizeof(double) );
					ptr += 7+m.nr*m.nc;
				}
			}
		}
	}
	flags_send[rank_send_now] = Flag_Send::finish_oar;
}

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::recv_data_process(
	const int rank_recv,
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_all,
	std::vector<std::valarray<double>> &iarps_irecv,
	std::vector<atomic<Flag_Recv>> &flags_recv,
	atomic_flag &lock_insert) const
{
	auto vector_empty = []( const std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &v ) -> bool
	{
		for( const auto &i : v )
			if(!i.empty())	return false;
		return true;
	};

	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> data_rank(GlobalV::NSPIN);
	double *ptr = ModuleBase::GlobalFunc::VECTOR_TO_PTR(iarps_irecv[rank_recv]);
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		const size_t recv_size = ptr[0];
		ptr += 1;
		const double * const ptr_end = ptr + recv_size;
		while(ptr<ptr_end)
		{
			const size_t iat1=ptr[0], iat2=ptr[1];
			const Abfs::Vector3_Order<int> box2 = {ptr[2],ptr[3],ptr[4]};
			const int nr=ptr[5], nc=ptr[6];
			ModuleBase::matrix &m = data_rank[is][iat1][iat2][box2];
			m.create(nr,nc);
			memcpy( m.c, ptr+7, nr*nc*sizeof(double) );
			ptr += 7+nr*nc;
		}
	}
	iarps_irecv[rank_recv].resize(0);
	flags_recv[rank_recv] = Flag_Recv::finish_iar;

	if(!vector_empty(data_rank))
	{
		while( lock_insert.test_and_set() );
		insert_data(data_rank, data_all);
		lock_insert.clear();
	}
}

void Exx_Abfs::Parallel::Communicate::Hexx::Allreduce2::insert_data(
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_rank,
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &data_all) const
{
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		auto &data_rank_is = data_rank[is];
		auto &data_all_is = data_all[is];
		for( auto &data_rank_A : data_rank_is )
		{
			const size_t iat1 = data_rank_A.first;
			if( !atom_in_2D_list[my_rank].first[iat1] )	continue;
			auto &data_all_A = data_all_is[iat1];
			for( auto &data_rank_B : data_rank_A.second )
			{
				const size_t iat2 = data_rank_B.first;
				if( !atom_in_2D_list[my_rank].second[iat2] )	continue;
				auto &data_all_B = data_all_A[iat2];
				for( auto &data_rank_C : data_rank_B.second )
				{
					const Abfs::Vector3_Order<int> &box2 = data_rank_C.first;
					auto &data_all_C = data_all_B[box2];
					if( data_all_C.c )
						data_all_C += data_rank_C.second;
					else
						data_all_C = std::move(data_rank_C.second);
				}
			}
		}
	}
}
#endif
