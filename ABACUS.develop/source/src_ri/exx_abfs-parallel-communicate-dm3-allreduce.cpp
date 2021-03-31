#include "exx_abfs-parallel-communicate-dm3.h"
#include "exx_abfs-parallel-communicate-function.h"

#include "src_pw/global.h"
#include "src_global/global_function.h"

#ifdef USE_CEREAL_SERIALIZATION
	#include "src_global/serialization_cereal.h"
#endif

#include <mpi.h>
#include <thread>
#include <algorithm>


#include "src_external/src_test/test_function.h"


void Exx_Abfs::Parallel::Communicate::DM3::Allreduce::init(
	const MPI_Comm &mpi_comm_in,
	const map<set<size_t>,set<size_t>> &H_atom_pairs_group)
{
	TITLE("Exx_Abfs::Parallel::Communicate::DM3::Allreduce::init");
	
	mpi_comm = mpi_comm_in;
	MPI_Comm_size(mpi_comm, &comm_sz);
	MPI_Comm_rank(mpi_comm, &my_rank);

//	const vector<pair<bool,bool>> atom_in_2D = get_atom_in_2D();
	H_atom_pairs_group_rank = get_H_atom_pairs_group_rank(H_atom_pairs_group);
	get_send_recv_size(H_atom_pairs_group_rank, H_atom_pairs_group, send_size_list, recv_size);

ofstream ofs(exx_lcao.test_dir.process+"dm3_"+TO_STRING(my_rank));
//ofs<<H_atom_pairs_group<<endl;
//ofs<<atom_in_2D<<endl;
//for(int rank=0; rank!=comm_sz; ++rank)
//{
//	ofs<<rank<<endl;
//	for(const auto &i:H_atom_pairs_group_rank[rank])
//		ofs<<"\t"<<i.first<<"\t"<<*i.second<<endl;
//}
ofs<<send_size_list<<endl;
ofs<<recv_size<<endl;
}


vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>>
	Exx_Abfs::Parallel::Communicate::DM3::Allreduce::a2D_to_exx(
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local) const
{
ofstream ofs(exx_lcao.test_dir.process+"dm3_"+TO_STRING(my_rank), ofstream::app);
//ofs<<data_local<<endl<<"@@@"<<endl;;

	TITLE("Exx_Abfs::Parallel::Communicate::DM3::Allreduce::a2D_to_exx");
	
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> data_all(NSPIN);

	vector<atomic<Flag_Send>> flags_send(comm_sz);
    vector<atomic<Flag_Recv>> flags_recv(comm_sz);
	init_flags(flags_send, flags_recv);

	vector<vector<double>> oarps_isend(comm_sz);
	vector<vector<double>> iarps_irecv(comm_sz);

	int rank_send_now = my_rank;
	auto rank_send_next = [&]()->int{ return (rank_send_now+1)%comm_sz; };

	atomic_flag lock_insert = ATOMIC_FLAG_INIT;
	vector<thread> threads;
	vector<MPI_Request>requests_isend(comm_sz);
	vector<MPI_Request>requests_irecv(comm_sz);
ofs<<__LINE__<<endl;

	while(!finish_judge(flags_send, flags_recv))
	{
		if(rank_send_next()!=my_rank && memory_enough(rank_send_next(), flags_send))
		{
			rank_send_now = rank_send_next();
ofs<<__LINE__<<"\t"<<rank_send_now<<endl;
			flags_send[rank_send_now] = Flag_Send::begin_oar;
			threads.push_back(std::thread(
				&Exx_Abfs::Parallel::Communicate::DM3::Allreduce::send_data_process, this,
				rank_send_now, std::cref(data_local), std::ref(oarps_isend), std::ref(flags_send) ));
ofs<<__LINE__<<"\t"<<rank_send_now<<endl;
		}
		for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
		{
			if( flags_send[rank_send] == Flag_Send::finish_oar )
			{
ofs<<__LINE__<<"\t"<<rank_send<<endl;
				if(MPI_Isend( VECTOR_TO_PTR(oarps_isend[rank_send]), oarps_isend[rank_send].size(), MPI_DOUBLE, rank_send, 0, mpi_comm, &requests_isend[rank_send] )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
				flags_send[rank_send] = Flag_Send::begin_isend;
ofs<<__LINE__<<"\t"<<rank_send<<"\t"<<oarps_isend[rank_send].size()<<endl;
			}
		}
		for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
		{
			if( flags_send[rank_send] == Flag_Send::begin_isend )
			{
				int flag_finish;
				if(MPI_Test( &requests_isend[rank_send], &flag_finish, MPI_STATUS_IGNORE )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
				if(flag_finish)
				{
ofs<<__LINE__<<"\t"<<rank_send<<endl;
					oarps_isend[rank_send].resize(0);	oarps_isend[rank_send].shrink_to_fit();
					flags_send[rank_send] = Flag_Send::finish_isend;
				}
			}
		}

		{
			MPI_Status status;
			int flag_message;
			if(MPI_Iprobe( MPI_ANY_SOURCE, 0, mpi_comm, &flag_message, &status )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
			if(flag_message)
			{
				int message_size;
				if(MPI_Get_count( &status, MPI_PACKED, &message_size )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
				const int rank_recv = status.MPI_SOURCE;
ofs<<__LINE__<<"\t"<<rank_recv<<endl;
				iarps_irecv[rank_recv].resize(message_size);
				if(MPI_Irecv( VECTOR_TO_PTR(iarps_irecv[rank_recv]), message_size, MPI_DOUBLE, rank_recv, 0, mpi_comm, &requests_irecv[rank_recv] )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
				flags_recv[rank_recv] = Flag_Recv::begin_irecv;
ofs<<__LINE__<<"\t"<<rank_recv<<"\t"<<message_size<<endl;
			}
		}
		for(int rank_recv=0; rank_recv!=comm_sz; ++rank_recv)
		{
			if(flags_recv[rank_recv] == Flag_Recv::begin_irecv)
			{
				int flag_finish = 0;
ofs<<__LINE__<<"\t"<<rank_recv<<"\t"<<flag_finish<<endl;
				if(MPI_Test( &requests_irecv[rank_recv], &flag_finish, MPI_STATUS_IGNORE )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
ofs<<__LINE__<<"\t"<<rank_recv<<"\t"<<flag_finish<<endl;
				if(flag_finish)
				{
ofs<<__LINE__<<"\t"<<rank_recv<<"\t"<<flag_finish<<endl;
					flags_recv[rank_recv] = Flag_Recv::begin_iar;
					threads.push_back(std::thread(
						&Exx_Abfs::Parallel::Communicate::DM3::Allreduce::recv_data_process, this,
						rank_recv, std::ref(data_all), std::ref(iarps_irecv), std::ref(flags_recv), std::ref(lock_insert) ));
ofs<<__LINE__<<"\t"<<rank_recv<<endl;
				}
			}
		}
	}
ofs<<__LINE__<<endl;

	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix*>>>> data_intersection = get_intersection(my_rank, data_local);
	while( lock_insert.test_and_set() );
	insert_data(data_intersection, data_all);
	lock_insert.clear();
ofs<<__LINE__<<endl;

	for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
	{
		if( flags_send[rank_send] == Flag_Send::begin_isend )
		{
			if(MPI_Wait( &requests_isend[rank_send], MPI_STATUS_IGNORE )!=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
			oarps_isend[rank_send].resize(0);		oarps_isend[rank_send].shrink_to_fit();
			flags_send[rank_send] = Flag_Send::finish_isend;
		}
	}
ofs<<__LINE__<<endl;

	for(std::thread &t : threads)
		t.join();
ofs<<__LINE__<<endl;
//ofs<<"@@@"<<endl<<data_all<<endl;

	return data_all;
}


/*
vector<pair<bool,bool>> Exx_Abfs::Parallel::Communicate::DM3::Allreduce::get_atom_in_2D() const
{
	vector<pair<bool,bool>> atom_in_2D(ucell.nat,{false,false});
	for(int iwt=0; iwt<NLOCAL; ++iwt)
	{
		const int iat = ucell.iwt2iat[iwt];
		if(KS_SOLVER=="genelpa")
		{
			if(ParaO.trace_loc_col[iwt]>=0)
				atom_in_2D[iat].first = true;
			if(ParaO.trace_loc_row[iwt]>=0)
				atom_in_2D[iat].second = true;
		}
		else
			throw domain_error(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	}
	return atom_in_2D;
}
*/


/*
vector<map<size_t,shared_ptr<set<size_t>>>> Exx_Abfs::Parallel::Communicate::DM3::Allreduce::get_H_atom_pairs_group_rank(
	const map<set<size_t>,set<size_t>> &H_atom_pairs_group,
	const vector<pair<bool,bool>> &atom_in_2D) const
{
	vector<int> H_atom_pairs_group_vector;
	for(const auto &pairs : H_atom_pairs_group)
	{
		H_atom_pairs_group_vector.push_back(pairs.first.size());
		H_atom_pairs_group_vector.push_back(pairs.second.size());
		H_atom_pairs_group_vector.insert(H_atom_pairs_group_vector.end(), pairs.first.begin(), pairs.first.end());
		H_atom_pairs_group_vector.insert(H_atom_pairs_group_vector.end(), pairs.second.begin(), pairs.second.end());
	}

	vector<int> rank_size(comm_sz);
	vector<int> rank_index(comm_sz);
	for(size_t i=0; i<rank_index.size(); ++i)
		rank_index[i] = i;
	const int myrank_size = H_atom_pairs_group_vector.size();
	MPI_Allgatherv(
		&myrank_size, 1, MPI_INT,
		VECTOR_TO_PTR(rank_size), VECTOR_TO_PTR(vector<int>(comm_sz,1)), VECTOR_TO_PTR(rank_index), MPI_INT,
		mpi_comm );

	size_t index_accumulate = 0;
	for(size_t i=0; i<rank_size.size(); ++i)
	{
		rank_index[i] = index_accumulate;
		index_accumulate += rank_size[i];
	}
	vector<int> H_atom_pairs_group_rank_vector(index_accumulate);
	MPI_Allgatherv(
		VECTOR_TO_PTR(H_atom_pairs_group_vector), H_atom_pairs_group_vector.size(), MPI_INT,
		VECTOR_TO_PTR(H_atom_pairs_group_rank_vector), VECTOR_TO_PTR(rank_size), VECTOR_TO_PTR(rank_index), MPI_INT,
		mpi_comm);

	vector<map<size_t,shared_ptr<set<size_t>>>> H_atom_pairs_group_rank(comm_sz);
	auto ptr = H_atom_pairs_group_rank_vector.begin();
	for(int rank=0; rank!=comm_sz; ++rank)
	{
		const auto ptr_begin = ptr;
		while(ptr < ptr_begin+rank_size[rank])
		{
			const size_t atoms1_size = *ptr++;
			const size_t atoms2_size = *ptr++;
			if(atoms1_size&&atoms2_size)
			{
				shared_ptr<set<size_t>> atoms2 = make_shared<set<size_t>>();
				for(int i=0; i<atoms1_size; ++i)
				{
					const size_t atom_require = *ptr++;
					if(atom_in_2D[atom_require].first)
						H_atom_pairs_group_rank[rank][atom_require] = atoms2;
				}
				for(int i=0; i<atoms2_size; ++i)
				{
					const size_t atom_require = *ptr++;
					if(atom_in_2D[atom_require].second)
						atoms2->insert(atom_require);
				}
			}
			else
			{
				ptr += atoms1_size + atoms2_size;
			}
		}
	}
	return H_atom_pairs_group_rank;
}
*/


vector<map<size_t,set<size_t>>> Exx_Abfs::Parallel::Communicate::DM3::Allreduce::get_H_atom_pairs_group_rank(
	const map<set<size_t>,set<size_t>> &H_atom_pairs_group) const
{
	constexpr int tag = 0;

	const vector<pair<vector<bool>,vector<bool>>> atom_in_2D_list = Exx_Abfs::Parallel::Communicate::Function::get_atom_in_2D_list(mpi_comm);

	vector<MPI_Request> request(comm_sz);
	vector<string> atom_send_str(comm_sz);
	for(int rank_tmp=my_rank; rank_tmp!=comm_sz+my_rank; ++rank_tmp)
	{
		const int rank = rank_tmp%comm_sz;
		map<size_t,set<size_t>> atom_required;
		for(const auto & H_atom_pairs : H_atom_pairs_group)
			for(const size_t iat1 : H_atom_pairs.first)
				if(atom_in_2D_list[rank].first[iat1])
					for(const size_t iat2 : H_atom_pairs.second)
						if(atom_in_2D_list[rank].second[iat2])
							atom_required[iat1].insert(iat2);
							
		#ifdef USE_CEREAL_SERIALIZATION
		{
			stringstream atom_send_ss;
			cereal::BinaryOutputArchive ar(atom_send_ss);
			ar(atom_required);
			atom_send_str[rank] = atom_send_ss.str();
		}
		#else
			throw invalid_argument(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		#endif
		if( MPI_Isend( atom_send_str[rank].c_str(), atom_send_str[rank].size(), MPI_CHAR, rank, tag, mpi_comm, &request[rank] ) !=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
	}

	vector<map<size_t,set<size_t>>> H_atom_pairs_group_rank(comm_sz);
	for(int i=0; i!=comm_sz; ++i)
	{
		MPI_Status status;
		if( MPI_Probe( MPI_ANY_SOURCE, tag, mpi_comm, &status ) !=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
		int size_recv;
		if( MPI_Get_count( &status, MPI_CHAR, &size_recv) !=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));

		vector<char> ss_buffer(size_recv);
		if( MPI_Recv( ss_buffer.data(), size_recv, MPI_CHAR, status.MPI_SOURCE, tag, mpi_comm, MPI_STATUS_IGNORE ) !=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
		stringstream atom_recv_ss;  
		atom_recv_ss.rdbuf()->pubsetbuf(ss_buffer.data(),size_recv);
		#ifdef USE_CEREAL_SERIALIZATION
		{
			cereal::BinaryInputArchive ar(atom_recv_ss);
			ar(H_atom_pairs_group_rank[status.MPI_SOURCE]);
		}	
		#else
			throw invalid_argument(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		#endif
	}

	if( MPI_Waitall( comm_sz, request.data(), MPI_STATUSES_IGNORE ) !=MPI_SUCCESS)	throw runtime_error(TO_STRING(__FILE__)+TO_STRING(__LINE__));
	return H_atom_pairs_group_rank;
}

void Exx_Abfs::Parallel::Communicate::DM3::Allreduce::get_send_recv_size(
	const vector<map<size_t,set<size_t>>> &H_atom_pairs_group_rank,
	const map<set<size_t>,set<size_t>> &H_atom_pairs_group,
	vector<size_t> &send_size_list, size_t &recv_size) const
{
	send_size_list.resize(comm_sz,0);
	for(int rank=0; rank!=comm_sz; ++rank)
	{
		for(const auto &H_pairs_tmp : H_atom_pairs_group_rank[rank])
		{
			const size_t iat1 = H_pairs_tmp.first;
			for(const size_t iat2 : H_pairs_tmp.second)
				send_size_list[rank] += ucell.atoms[ucell.iat2it[iat1]].nw * ucell.atoms[ucell.iat2it[iat2]].nw;
		}
		send_size_list[rank] *= sizeof(double);
	}

	recv_size = 0;
	for(const auto &pairs : H_atom_pairs_group)
	{
		size_t nw1=0, nw2=0;
		for(size_t iat1 : pairs.first)
			nw1 += ucell.atoms[ucell.iat2it[iat1]].nw;
		for(size_t iat2 : pairs.second)
			nw2 += ucell.atoms[ucell.iat2it[iat2]].nw;
		recv_size += nw1*nw1;
	}
	recv_size *= sizeof(double);
}


void Exx_Abfs::Parallel::Communicate::DM3::Allreduce::init_flags(
	vector<atomic<Flag_Send>> &flags_send,
	vector<atomic<Flag_Recv>> &flags_recv) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::DM3::Allreduce::init_flags");
	
	for(atomic<Flag_Send> &flag_send : flags_send)
		flag_send = Flag_Send::undo;
	flags_send[my_rank] = Flag_Send::finish_isend;
	
	for(atomic<Flag_Recv> &flag_recv : flags_recv)
		flag_recv = Flag_Recv::undo;
	flags_recv[my_rank] = Flag_Recv::finish_iar;
}


bool Exx_Abfs::Parallel::Communicate::DM3::Allreduce::finish_judge(
	const vector<atomic<Flag_Send>> &flags_send,
	const vector<atomic<Flag_Recv>> &flags_recv) const
{
	for(int rank_send=0; rank_send!=comm_sz; ++rank_send)
		if((flags_send[rank_send]!=Flag_Send::begin_isend) && (flags_send[rank_send]!=Flag_Send::finish_isend))
			return false;
	for(int rank_recv=0; rank_recv!=comm_sz; ++rank_recv)
		if((flags_recv[rank_recv]!=Flag_Recv::begin_iar) && (flags_recv[rank_recv]!=Flag_Recv::finish_iar))
			return false;
	return true;
}


bool Exx_Abfs::Parallel::Communicate::DM3::Allreduce::memory_enough(
	const int rank_send_next,
	const vector<atomic<Flag_Send>> &flags_send) const
{
	size_t memory_need = recv_size;
	for(int rank_send=0; rank_send<comm_sz; ++rank_send)
		if(flags_send[rank_send]==Flag_Send::begin_oar)
			memory_need += send_size_list[rank_send];
	const size_t memory_available = MemAvailable()*1024;
	return (memory_available-memory_need>send_size_list[rank_send_next]) ? true : false;
}


/*
void Exx_Abfs::Parallel::Communicate::DM3::Allreduce::send_data_process(
	const int rank_send_now,
	const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local,
	vector<vector<double>> &oarps_isend,
	vector<atomic<Flag_Send>> &flags_send) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::DM3::Allreduce::send_data_process");
	
	vector<double> &oarp = oarps_isend[rank_send_now];
	for( int is=0; is!=NSPIN; ++is )
	{
		oarp.push_back(-999);		const size_t index_size = oarp.size()-1;
		auto ptr_data_A = data_local[is].begin();
		auto ptr_atom_A = H_atom_pairs_group_rank[rank_send_now].begin();
		while(ptr_data_A!=data_local[is].end() && ptr_atom_A!=H_atom_pairs_group_rank[rank_send_now].end())
		{
			if(ptr_data_A->first < ptr_atom_A->first)		++ptr_data_A;
			else if(ptr_data_A->first > ptr_atom_A->first)	++ptr_atom_A;
			else
			{
				const size_t iat1 = ptr_data_A->first;
				auto ptr_data_B = ptr_data_A->second.begin();
				auto ptr_atom_B = ptr_atom_A->second->begin();
				while(ptr_data_B!=ptr_data_A->second.end() && ptr_atom_B!=ptr_atom_A->second->end())
				{
					const size_t iat2 = ptr_data_B->first;
					if(ptr_data_B->first < *ptr_atom_B)			++ptr_data_B;
					else if(ptr_data_B->first > *ptr_atom_B)	++ptr_atom_B;
					else
					{
						for(const auto &data_C : ptr_data_B->second)
						{
							const Abfs::Vector3_Order<int> &box2 = data_C.first;
							const matrix &m = data_C.second;
							oarp.push_back(iat1);
							oarp.push_back(iat2);
							oarp.push_back(box2.x);	oarp.push_back(box2.y);	oarp.push_back(box2.z);
							oarp.push_back(m.nr);	oarp.push_back(m.nc);
							oarp.insert(oarp.end(), m.c, m.c+m.nr*m.nc);
						}
						++ptr_data_B;	++ptr_atom_B;
					}
				}
				++ptr_data_A;	++ptr_atom_A;
			}
		}
		oarp[index_size] = oarp.size() - (index_size+1);
	}
	flags_send[rank_send_now] = Flag_Send::finish_oar;
}
*/

void Exx_Abfs::Parallel::Communicate::DM3::Allreduce::send_data_process(
	const int rank_send_now,
	const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local,
	vector<vector<double>> &oarps_isend,
	vector<atomic<Flag_Send>> &flags_send) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::DM3::Allreduce::send_data_process");
	
	const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix*>>>> data_intersection = get_intersection(
		rank_send_now,
		const_cast<vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>>&>(data_local));

	vector<double> &oarp = oarps_isend[rank_send_now];
	for( int is=0; is!=NSPIN; ++is )
	{
		oarp.push_back(-999);		const size_t index_size = oarp.size()-1;
		for(const auto data_A : data_intersection[is])
		{
			const size_t iat1 = data_A.first;
			for(const auto data_B : data_A.second)
			{
				const size_t iat2 = data_B.first;
				for(const auto data_C : data_B.second)
				{
					const Abfs::Vector3_Order<int> &box2 = data_C.first;
					const matrix &m = *data_C.second;
					oarp.push_back(iat1);
					oarp.push_back(iat2);
					oarp.push_back(box2.x);	oarp.push_back(box2.y);	oarp.push_back(box2.z);
					oarp.push_back(m.nr);	oarp.push_back(m.nc);
					oarp.insert(oarp.end(), m.c, m.c+m.nr*m.nc);					
				}
			}
		}
		oarp[index_size] = oarp.size() - (index_size+1);
	}
	flags_send[rank_send_now] = Flag_Send::finish_oar;
}

vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix*>>>> Exx_Abfs::Parallel::Communicate::DM3::Allreduce::get_intersection(
	const int rank_send_now,
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_local) const
{
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix*>>>> data_intersection(NSPIN);
	for( int is=0; is!=NSPIN; ++is )
	{
		auto ptr_data_A = data_local[is].begin();
		auto ptr_atom_A = H_atom_pairs_group_rank[rank_send_now].begin();
		while(ptr_data_A!=data_local[is].end() && ptr_atom_A!=H_atom_pairs_group_rank[rank_send_now].end())
		{
			if(ptr_data_A->first < ptr_atom_A->first)		++ptr_data_A;
			else if(ptr_data_A->first > ptr_atom_A->first)	++ptr_atom_A;
			else
			{
				const size_t iat1 = ptr_data_A->first;
				auto ptr_data_B = ptr_data_A->second.begin();
				auto ptr_atom_B = ptr_atom_A->second.begin();
				while(ptr_data_B!=ptr_data_A->second.end() && ptr_atom_B!=ptr_atom_A->second.end())
				{
					const size_t iat2 = ptr_data_B->first;
					if(ptr_data_B->first < *ptr_atom_B)			++ptr_data_B;
					else if(ptr_data_B->first > *ptr_atom_B)	++ptr_atom_B;
					else
					{
						for(auto &data_C : ptr_data_B->second)
						{
							const Abfs::Vector3_Order<int> &box2 = data_C.first;
							data_intersection[is][iat1][iat2][box2] = &data_C.second;
						}
						++ptr_data_B;	++ptr_atom_B;
					}
				}
				++ptr_data_A;	++ptr_atom_A;
			}
		}
	}
	return data_intersection;
}

void Exx_Abfs::Parallel::Communicate::DM3::Allreduce::recv_data_process(
	const int rank_recv,
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_all,
	vector<vector<double>> &iarps_irecv,
	vector<atomic<Flag_Recv>> &flags_recv,
	atomic_flag &lock_insert) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::DM3::Allreduce::recv_data_process");
	
	auto vector_empty = []( const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &v ) -> bool
	{
		for( const auto &i : v )
			if(!i.empty())	return false;
		return true;
	};

	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> data_rank(NSPIN);
	auto ptr = iarps_irecv[rank_recv].begin();
	for( int is=0; is!=NSPIN; ++is )
	{
		const size_t recv_size = ptr[0];
		ptr += 1;
		auto ptr_end = ptr + recv_size;
		while(ptr<ptr_end)
		{
			const size_t iat1=ptr[0], iat2=ptr[1];
			const Abfs::Vector3_Order<int> box2 = {ptr[2],ptr[3],ptr[4]};
			const int nr=ptr[5], nc=ptr[6];
			matrix &m = data_rank[is][iat1][iat2][box2];
			m.create(nr,nc);
			copy( ptr+7, ptr+7+nr*nc, m.c );
			ptr += 7+nr*nc;
		}
	}
	iarps_irecv[rank_recv].resize(0);		iarps_irecv[rank_recv].shrink_to_fit();
	flags_recv[rank_recv] = Flag_Recv::finish_iar;

	if(!vector_empty(data_rank))
	{
		while( lock_insert.test_and_set() );
		insert_data(data_rank, data_all);
		lock_insert.clear();
	}
}


matrix& Exx_Abfs::Parallel::Communicate::DM3::Allreduce::get_matrix(matrix&m)const{ return m; }
matrix& Exx_Abfs::Parallel::Communicate::DM3::Allreduce::get_matrix(matrix*pm)const{ return *pm; }


template<typename M>
void Exx_Abfs::Parallel::Communicate::DM3::Allreduce::insert_data(
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,M>>>> &data_rank,
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &data_all) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::DM3::Allreduce::insert_data");
	for( int is=0; is!=NSPIN; ++is )
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
					matrix &m = get_matrix(data_rank_C.second);
					if( data_all_C.c )
						data_all_C += m;
					else
						data_all_C = std::move(m);
				}
			}
		}
	}
}
