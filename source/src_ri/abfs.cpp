#include "abfs.h"
#include "abfs-vector3_order.h"
#include "abfs-template.h"

#include "exx_abfs-inverse_matrix_double.h"
#include "../src_pw/global.h"
#include "../module_base/global_function.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

#include <fstream>		// Peize Lin test
#include <iomanip>		// Peize Lin test
#include "../src_external/src_test/src_global/matrix-test.h"		// Peize Lin test
#include "../src_external/src_test/test_function.h"

std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> Abfs::cal_Cs(
	const set<size_t> &atom_centres,
	const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
	const Exx_Abfs::Matrix_Orbs21 &m_abfslcaos_lcaos,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_abfs,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_lcaos,
	const double threshold,
	std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Cws,
	std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Vws )
{
	ModuleBase::TITLE("Abfs","cal_Cs");
	pthread_rwlock_t rwlock_Cw;	pthread_rwlock_init(&rwlock_Cw,NULL);
	pthread_rwlock_t rwlock_Vw;	pthread_rwlock_init(&rwlock_Vw,NULL);
	
	const std::vector<size_t> atom_centres_vector( atom_centres.begin(), atom_centres.end() );
	const std::vector<std::map<size_t,std::vector<Abfs::Vector3_Order<int>>>> adjs = get_adjs();
	
	// pre-cal Vws on same atom, speed up DPcal_V() in DPcal_C()
	std::vector<std::shared_ptr<ModuleBase::matrix>> Vs_same_atom(GlobalC::ucell.ntype);

	#ifdef __MPI //liyuanbo 2022/2/23
	for(size_t it=0; it!=GlobalC::ucell.ntype; ++it)
		Vs_same_atom[it] = DPcal_V( it,it,{0,0,0}, m_abfs_abfs, index_abfs, 0,true, rwlock_Vw,Vws );
	#endif

#ifdef __MKL
    const int mkl_threads = mkl_get_max_threads();
	mkl_set_num_threads(std::max(1UL,mkl_threads/atom_centres_vector.size()));
#endif
	
	std::map<size_t,std::map<size_t,std::map<Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> Cs;
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for( int i_iat1=0; i_iat1<atom_centres_vector.size(); ++i_iat1 )
	{
		const size_t iat1 = atom_centres_vector[i_iat1];
		const size_t it1 = GlobalC::ucell.iat2it[iat1];
		const size_t ia1 = GlobalC::ucell.iat2ia[iat1];
		const Vector3_Order<double> tau1( GlobalC::ucell.atoms[it1].tau[ia1] );
		
		for( const auto & atom2 : adjs[iat1] )
		{
			const int iat2 = atom2.first;
			const int it2 = GlobalC::ucell.iat2it[iat2];
			const int ia2 = GlobalC::ucell.iat2ia[iat2];
			const Vector3_Order<double> tau2( GlobalC::ucell.atoms[it2].tau[ia2] );
			for( const ModuleBase::Vector3<int> &box2 : atom2.second )
			{
//				std::cout<<"cal_Cs\t"<<iat1<<"\t"<<iat2<<"\t"<<box2<<std::endl;
			#ifdef __MPI //liyuanbo 2022/2/23
				const std::shared_ptr<ModuleBase::matrix> C = DPcal_C( 
					it1, it2, -tau1+tau2+(box2*GlobalC::ucell.latvec), 
					m_abfs_abfs, m_abfslcaos_lcaos, index_abfs, index_lcaos, 
					threshold, true, rwlock_Cw, rwlock_Vw, Cws, Vws );
#ifdef _OPENMP
				#pragma omp critical(Abfs_cal_Cs)
#endif
				Cs[iat1][iat2][box2] = C;
			#endif
			}
		}
	}
#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
#endif
	Abfs::delete_threshold_ptrs( Cs, threshold );
	Vs_same_atom.clear();
	Abfs::delete_empty_ptrs( Vws );
	pthread_rwlock_destroy(&rwlock_Cw);
	pthread_rwlock_destroy(&rwlock_Vw);
	return Cs;
}

/*std::map<size_t,std::map<size_t,std::shared_ptr<matrix>>> 
	Abfs::cal_Vps(
		const set<std::pair<size_t,size_t>> &atom_pairs,
		const std::vector<Vector3_Order<int>> &Coulomb_potential_boxes,
		const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
		const ModuleBase::Element_Basis_Index::Index &index_abfs,
		std::map<size_t,std::map<size_t,std::map<Vector3_Order<int>,std::shared_ptr<matrix>>>> &Vs,
		std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<matrix>>>> &Vws )
{
	std::map<size_t,std::map<size_t,std::shared_ptr<matrix>>> Vps;
	
	for( const std::pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = GlobalC::ucell.iat2it[iat1];
		const size_t ia1 = GlobalC::ucell.iat2ia[iat1];
		const size_t it2 = GlobalC::ucell.iat2it[iat2];
		const size_t ia2 = GlobalC::ucell.iat2ia[iat2];
		const Vector3_Order<double> tau1 = GlobalC::ucell.atoms[it1].tau[ia1];
		const Vector3_Order<double> tau2 = GlobalC::ucell.atoms[it2].tau[ia2];
		
		for( const Vector3_Order<int> &box2 : Coulomb_potential_boxes )
		{
			const std::vector<std::shared_ptr<matrix>> Vs_tmp = DPcal_V( it1, it2, -tau1+tau2+(box2*GlobalC::ucell.latvec), m_abfs_abfs, index_abfs, Vws );
			Vs[iat1][iat2][box2] = Vs_tmp[0];	Vs[iat2][iat1][-box2] = Vs_tmp[1];
		}
		
		Vps[iat1][iat2] = make_shared<matrix>( Vs[iat1][iat2][{0,0,0}]->nr, Vs[iat1][iat2][{0,0,0}]->nc );
		for( const Vector3_Order<int> &box2 : Coulomb_potential_boxes )
			*Vps[iat1][iat2] += *Vs[iat1][iat2][box2];
	}
	return Vps;
}
*/

/*
std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::shared_ptr<matrix>>>> 
	Abfs::cal_Vps(
		const std::vector<std::pair<size_t,size_t>> &atom_pairs,
		const std::vector<Vector3_Order<int>> &Coulomb_potential_boxes,
		const std::vector<Vector3_Order<int>> &Born_von_Karman_boxes,
		const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_abfs,
		std::map<size_t,std::map<size_t,std::map<Vector3_Order<int>,std::shared_ptr<matrix>>>> &Vs,
		std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<matrix>>>> &Vws )
{
	ModuleBase::TITLE("Abfs","cal_Vps");	
	
	std::map<size_t,std::map<size_t,std::map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Vps;
	
	for( const std::pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = GlobalC::ucell.iat2it[iat1];
		const size_t ia1 = GlobalC::ucell.iat2ia[iat1];
		const size_t it2 = GlobalC::ucell.iat2it[iat2];
		const size_t ia2 = GlobalC::ucell.iat2ia[iat2];
		const Vector3_Order<double> tau1 = GlobalC::ucell.atoms[it1].tau[ia1];
		const Vector3_Order<double> tau2 = GlobalC::ucell.atoms[it2].tau[ia2];
		
		for( const Vector3_Order<int> &box2 : Coulomb_potential_boxes )
		{
std::cout<<"cal_Vs\t"<<iat1<<"\t"<<iat2<<"\t"<<box2<<std::endl;
			const std::vector<std::shared_ptr<matrix>> Vs_tmp = DPcal_V( it1, it2, -tau1+tau2+(box2*GlobalC::ucell.latvec), m_abfs_abfs, index_abfs, Vws );
			Vs[iat1][iat2][box2] = Vs_tmp[0];	Vs[iat2][iat1][-box2] = Vs_tmp[1];
		}
		
		std::vector<ModuleBase::ComplexMatrix> Vkps( GlobalC::kv.nks/GlobalV::NSPIN, {Vs[iat1][iat2][{0,0,0}]->nr, Vs[iat1][iat2][{0,0,0}]->nc} );
		for( size_t ik=0; ik!=GlobalC::kv.nks/GlobalV::NSPIN; ++ik )
			for( const Vector3_Order<int> &box2 : Coulomb_potential_boxes )
				Vkps[ik] += static_cast<ModuleBase::ComplexMatrix>(*Vs[iat1][iat2][box2]) * exp(-ModuleBase::TWO_PI*IMAG_UNIT*(GlobalC::kv.kvec_c[ik]*(box2*GlobalC::ucell.latvec)));
		
		for( const Vector3_Order<int> &box2 : Born_von_Karman_boxes )
		{
			ModuleBase::ComplexMatrix Vps_tmp ( Vkps[0].nr, Vkps[0].nc );
			for( size_t ik=0; ik!=GlobalC::kv.nks/GlobalV::NSPIN; ++ik )
				Vps_tmp += Vkps[ik] * GlobalC::kv.wk[ik]*(0.5*GlobalV::NSPIN) * exp(ModuleBase::TWO_PI*IMAG_UNIT*(GlobalC::kv.kvec_c[ik]*(box2*GlobalC::ucell.latvec)));
			Vps[iat1][iat2][box2] = make_shared<matrix>( Vps_tmp.real() );
		}
	}
	return Vps;
}
*/

std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> Abfs::cal_Vs(
	const std::vector<std::pair<size_t,size_t>> &atom_pairs,
	const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_abfs,
	const double rmesh_times,
	const double threshold,
	std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Vws )
{
	ModuleBase::TITLE("Abfs","cal_Vs");
	pthread_rwlock_t rwlock_Vw;	pthread_rwlock_init(&rwlock_Vw,NULL);
	std::vector<Abfs::Vector3_Order<int>> Coulomb_potential_boxes = get_Coulomb_potential_boxes(rmesh_times);

	std::map<size_t,std::map<size_t,std::map<Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> Vs;
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for( int i_atom_pair=0; i_atom_pair<atom_pairs.size(); ++i_atom_pair )
	{
		const size_t iat1 = atom_pairs[i_atom_pair].first;
		const size_t iat2 = atom_pairs[i_atom_pair].second;
		const size_t it1 = GlobalC::ucell.iat2it[iat1];
		const size_t ia1 = GlobalC::ucell.iat2ia[iat1];
		const size_t it2 = GlobalC::ucell.iat2it[iat2];
		const size_t ia2 = GlobalC::ucell.iat2ia[iat2];
		const Vector3_Order<double> tau1 = GlobalC::ucell.atoms[it1].tau[ia1];
		const Vector3_Order<double> tau2 = GlobalC::ucell.atoms[it2].tau[ia2];
		const double Rcut = std::min( GlobalC::ORB.Phi[it1].getRcut()*rmesh_times+GlobalC::ORB.Phi[it2].getRcut(), GlobalC::ORB.Phi[it1].getRcut()+GlobalC::ORB.Phi[it2].getRcut()*rmesh_times );
		
		for( const Vector3_Order<int> &box2 : Coulomb_potential_boxes )
		{
			const Vector3_Order<double> delta_R = -tau1+tau2+(box2*GlobalC::ucell.latvec);
			if( delta_R.norm()*GlobalC::ucell.lat0 < Rcut )
			{
			#ifdef __MPI //liyuanbo 2022/2/23
//				std::cout<<"cal_Vs\t"<<iat1<<"\t"<<iat2<<"\t"<<box2<<"\t"<<delta_R<<"\t"<<delta_R.norm()<<"\t"<<delta_R.norm()*GlobalC::ucell.lat0<<"\t"<<GlobalC::ORB.Phi[it1].getRcut()*rmesh_times+GlobalC::ORB.Phi[it2].getRcut()<<std::endl;
				const std::shared_ptr<ModuleBase::matrix> V = DPcal_V( 
					it1, it2, delta_R, 
					m_abfs_abfs, index_abfs, 
					threshold, true, rwlock_Vw, Vws );
#ifdef _OPENMP					
				#pragma omp critical(Abfs_cal_Vs)
#endif
				Vs[iat1][iat2][box2] = V;
			#endif
			}
		}
	}
	Abfs::delete_threshold_ptrs( Vs, threshold );
	pthread_rwlock_destroy(&rwlock_Vw);
	return Vs;
}

std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>> Abfs::cal_mps(
	const Abfs::Vector3_Order<int> &Born_von_Karman_period,
	const std::map<Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>> &ms )
{
	std::map< Vector3_Order<int>, std::vector<Vector3_Order<int>> > indexs;
	for( const auto & m : ms )
	{
		const Vector3_Order<int> & box = m.first;
		indexs[box%Born_von_Karman_period].push_back(box);
	}
	
	#if TEST_EXX_LCAO==1
	{
		static int istep=0;
		std::ofstream ofs( "mps_index_"+ModuleBase::GlobalFunc::TO_STRING(iat1)+"_"+ModuleBase::GlobalFunc::TO_STRING(iat2)+"_"+ModuleBase::GlobalFunc::TO_STRING(istep++) );
		for( const auto index : indexs )
		{
			ofs<<index.first<<std::endl;
			for( const Vector3_Order<int> & box : index.second )
				ofs<<"\t"<<box<<std::endl;
		}
		ofs.close();
	}
	#elif TEST_EXX_LCAO==-1
		#error
	#endif
	
	std::map<Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>> mps;
	for( const auto & index : indexs )
	{
		const Vector3_Order<int> & boxp = index.first;
		assert(index.second.size()>=1);
		if( index.second.size()==1 )
			mps[boxp] = ms.at(index.second[0]);
		else
		{
			mps[boxp] = make_shared<ModuleBase::matrix>();
			*mps[boxp] = *ms.at(index.second[0]);
			for( size_t i=1; i<index.second.size(); ++i )
				*mps[boxp] += *ms.at(index.second[i]);
		}
	}
	return mps;
}

std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> Abfs::cal_mps(
	const Abfs::Vector3_Order<int> &Born_von_Karman_period,
	const std::map<size_t,std::map<size_t,std::map<Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> &ms )
{
	ModuleBase::TITLE("Abfs","cal_mps");
	std::map<size_t,std::map<size_t,std::map<Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> mps;
	for( const auto & m1s : ms )
		for( const auto & m12s : m1s.second )
			mps[m1s.first][m12s.first] = cal_mps( Born_von_Karman_period, m12s.second );
	return mps;
}

#ifdef __MPI //liyuanbo 2022/2/23
std::shared_ptr<ModuleBase::matrix> Abfs::DPcal_C( 
	const size_t &it1, 
	const size_t &it2, 
	const Vector3_Order<double> &R, 
	const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
	const Exx_Abfs::Matrix_Orbs21 &m_abfslcaos_lcaos,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_abfs,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_lcaos,
	const double threshold,
	const bool writable,
	pthread_rwlock_t &rwlock_Cw,
	pthread_rwlock_t &rwlock_Vw,
	std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Cws,
	std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Vws )
{
	// Attention: 
	// 在计算 C[it1,it1,0][it2,R] 后，额外计算 C[it2,it2,0][it1,-R] 几乎不耗时。
	// if 仅保存 C[it1,it1,0][it2,R]，不保存 C[it2,it2,0][it1,-R]：
	//     if C[it2,it2,0][it1,-R] 在本核上亦需计算：
	//         浪费时间重算一遍
	//     if C[it2,it2,0][it1,-R] 在本核上不需计算：
	//         好
	// if 同时保存 C[it1,it1,0][it2,R]、C[it2,it2,0][it1,-R]：
	//     if C[it2,it2,0][it1,-R] 在本核上亦需计算：
	//         好
	//     if C[it2,it2,0][it1,-R] 在本核上不需计算：
	//         浪费空间（可在算完全部Cs后根据需要删除无需的Cs，释放空间？）
	// 根据后续情况选择权衡。
	
//	ModuleBase::TITLE("Abfs","DPcal_C");
	pthread_rwlock_rdlock(&rwlock_Cw);
	const std::weak_ptr<ModuleBase::matrix> * const Cws_ptr   = static_cast<const std::weak_ptr<ModuleBase::matrix> * const>( ModuleBase::GlobalFunc::MAP_EXIST( Cws, it1, it2, R ) );
	pthread_rwlock_unlock(&rwlock_Cw);
	
	if( Cws_ptr && !Cws_ptr->expired() )
		return Cws_ptr->lock();
	else
	{
//		std::cout<<"DPcal_C\t"<<it1<<"\t"<<it2<<"\t"<<R<<std::endl;
		if( (ModuleBase::Vector3<double>(0,0,0)==R) && (it1==it2) )
		{
			const std::shared_ptr<ModuleBase::matrix> A = 
				make_shared<ModuleBase::matrix>( m_abfslcaos_lcaos.cal_overlap_matrix( it1,it2,0,0,index_abfs,index_lcaos,index_lcaos,Exx_Abfs::Matrix_Orbs21::Matrix_Order::A2B_A1 ) );
			const std::shared_ptr<ModuleBase::matrix> V = DPcal_V(it1,it2,{0,0,0}, m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
			const std::shared_ptr<ModuleBase::matrix> L = cal_I(V);
			std::shared_ptr<ModuleBase::matrix> C = make_shared<ModuleBase::matrix>( 0.5 * *A * *L );
			if(C->absmax()<=threshold)	C->create(0,0);
			if(writable)
			{
				pthread_rwlock_wrlock(&rwlock_Cw);
				Cws[it1][it2][R] = C;
				pthread_rwlock_unlock(&rwlock_Cw);
			}
			return C;
		}
		else
		{
			const std::vector<std::shared_ptr<ModuleBase::matrix>> A = 
				{ make_shared<ModuleBase::matrix>( m_abfslcaos_lcaos.cal_overlap_matrix(it1,it2,0,R ,index_abfs,index_lcaos,index_lcaos,Exx_Abfs::Matrix_Orbs21::Matrix_Order::A2B_A1) ) ,
				  make_shared<ModuleBase::matrix>( m_abfslcaos_lcaos.cal_overlap_matrix(it2,it1,0,-R,index_abfs,index_lcaos,index_lcaos,Exx_Abfs::Matrix_Orbs21::Matrix_Order::BA2_A1) ) };

			const std::shared_ptr<ModuleBase::matrix> V_00 = DPcal_V(it1,it1,{0,0,0}, m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
			const std::shared_ptr<ModuleBase::matrix> V_01 = DPcal_V(it1,it2,R,       m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
			const std::shared_ptr<ModuleBase::matrix> V_10 = DPcal_V(it2,it1,-R,      m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
			const std::shared_ptr<ModuleBase::matrix> V_11 = DPcal_V(it2,it2,{0,0,0}, m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
			const std::vector<std::vector<std::shared_ptr<ModuleBase::matrix>>> V = 
				{{ V_00, V_01 },
				 { V_10, V_11 }};

			const std::vector<std::vector<std::shared_ptr<ModuleBase::matrix>>> L = cal_I(V);

			std::shared_ptr<ModuleBase::matrix> C = make_shared<ModuleBase::matrix>( *A[0] * *L[0][0] + *A[1] * *L[1][0] );
			if(C->absmax()<=threshold)	C->create(0,0);
			if(writable)
			{
				pthread_rwlock_wrlock(&rwlock_Cw);
				Cws[it1][it2][R] = C;
				pthread_rwlock_unlock(&rwlock_Cw);
			}
			return C;
		}
	}
}

std::shared_ptr<ModuleBase::matrix> Abfs::DPcal_V( 
	const size_t &it1, 
	const size_t &it2, 
	const Vector3_Order<double> &R, 
	const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_abfs,
	const double threshold,
	const bool writable,
	pthread_rwlock_t &rwlock_Vw,
	std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Vws)
{
//	ModuleBase::TITLE("Abfs","DPcal_V");
	pthread_rwlock_rdlock(&rwlock_Vw);
	const std::weak_ptr<ModuleBase::matrix> * const Vws12_ptr = static_cast<const std::weak_ptr<ModuleBase::matrix> * const>( ModuleBase::GlobalFunc::MAP_EXIST( Vws, it1, it2, R ) );
	const std::weak_ptr<ModuleBase::matrix> * const Vws21_ptr = static_cast<const std::weak_ptr<ModuleBase::matrix> * const>( ModuleBase::GlobalFunc::MAP_EXIST( Vws, it2, it1, -R ) );
	pthread_rwlock_unlock(&rwlock_Vw);
	
	if( Vws12_ptr && !Vws12_ptr->expired() )
		return Vws12_ptr->lock();
	else if( Vws21_ptr && !Vws21_ptr->expired() )
	{
		std::shared_ptr<ModuleBase::matrix> VT = make_shared<ModuleBase::matrix>( transpose(*Vws21_ptr->lock()) );
		if(VT->absmax()<=threshold)	VT->create(0,0);
		if(writable)
		{
			pthread_rwlock_wrlock(&rwlock_Vw);
			Vws[it1][it2][R] = VT;
			pthread_rwlock_unlock(&rwlock_Vw);
		}
		return VT;
	}
	else
	{
//		std::cout<<"DPcal_V\t"<<it1<<"\t"<<it2<<"\t"<<R<<std::endl;
		std::shared_ptr<ModuleBase::matrix> V = make_shared<ModuleBase::matrix>( m_abfs_abfs.cal_overlap_matrix(it1,it2,0,R,index_abfs,index_abfs) );
		if(V->absmax()<=threshold)	V->create(0,0);
		if(writable)
		{
			pthread_rwlock_wrlock(&rwlock_Vw);
			Vws[it1][it2][R] = V;
			pthread_rwlock_unlock(&rwlock_Vw);
		}
		return V;
	}
}
#endif

std::map<size_t,std::vector<Abfs::Vector3_Order<int>>> Abfs::get_adjs( const size_t &iat )
{
//	ModuleBase::TITLE("Abfs","get_adjs");
	const int it = GlobalC::ucell.iat2it[iat];
	const int ia = GlobalC::ucell.iat2ia[iat];
	const ModuleBase::Vector3<double> &tau = GlobalC::ucell.atoms[it].tau[ia];
	
	std::map<size_t,std::vector<Vector3_Order<int>>> adjs;
	GlobalC::GridD.Find_atom(GlobalC::ucell,  tau, it, ia );
	for( int ad=0; ad<GlobalC::GridD.getAdjacentNum()+1; ++ad )
	{
		const size_t it_ad = GlobalC::GridD.getType(ad);
		const size_t ia_ad = GlobalC::GridD.getNatom(ad);
		const int iat_ad = GlobalC::ucell.itia2iat(it_ad,ia_ad);
		const ModuleBase::Vector3<int> box_ad = GlobalC::GridD.getBox(ad);
		const ModuleBase::Vector3<double> tau_ad = GlobalC::GridD.getAdjacentTau(ad);
		
		if( (tau-tau_ad).norm()*GlobalC::ucell.lat0 < GlobalC::ORB.Phi[it].getRcut()+GlobalC::ORB.Phi[it_ad].getRcut() )
			adjs[iat_ad].push_back(box_ad);
	}
	return adjs;
}

std::vector<std::map<size_t,std::vector<Abfs::Vector3_Order<int>>>> Abfs::get_adjs()
{
	std::vector<std::map<size_t,std::vector<Abfs::Vector3_Order<int>>>> adjs(GlobalC::ucell.nat);
	for( size_t iat=0; iat!=GlobalC::ucell.nat; ++iat )
		adjs[iat] = Abfs::get_adjs(iat);
	return adjs;
}

/*
set<std::pair<size_t,size_t>> Abfs::get_H_pairs_core( const std::vector<std::pair<size_t,size_t>> &atom_pairs )
{
	ModuleBase::TITLE("Exx_Lcao","allocate_Hexx");

	set<std::pair<size_t,size_t>> H_atom_pairs_core;
	for( const std::pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;

		const std::map<size_t,std::vector<Abfs::Vector3_Order<int>>> adj1s = Abfs::get_adjs(iat1);
		const std::map<size_t,std::vector<Abfs::Vector3_Order<int>>> adj2s = Abfs::get_adjs(iat2);

		H_atom_pairs_core.insert({iat1,iat2});
		for( const auto & atom3 : adj1s )
		{
			const size_t iat3 = atom3.first;
			H_atom_pairs_core.insert({iat3,iat2});
			for( const auto & atom4 : adj2s )
			{
				const size_t iat4 = atom4.first;
				H_atom_pairs_core.insert({iat1,iat4});
				H_atom_pairs_core.insert({iat3,iat4});
			}
		}
	}
	return H_atom_pairs_core;
}
*/

std::map<set<size_t>,set<size_t>> Abfs::get_H_pairs_core_group( const std::vector<std::pair<size_t,size_t>> &atom_pairs )
{
	ModuleBase::TITLE("Abfs","get_H_pairs_core_group");
	
	const std::vector<std::map<size_t,std::vector<Abfs::Vector3_Order<int>>>> adjs = Abfs::get_adjs();
	
	auto get_set_adjs = [&adjs]( const set<size_t> &as ) -> set<size_t>
	{
		set<size_t> aRs;
		for( const size_t a : as )
			for( const auto &aR : adjs[a] )
				aRs.insert(aR.first);
		return aRs;
	};
	
	// {(0,1),(3,5),(0,4),(7,8),(6,5)}
	const std::vector<std::pair<size_t,size_t>> & a1_a2 = atom_pairs;
	
	// => {1:{0}, 5:{3,6}, 4:{0}, 8:{7}}
	std::map<size_t,set<size_t>> a2_a1s;
	for( const auto & a1_a2_i : a1_a2 )
		a2_a1s[ a1_a2_i.second ].insert( a1_a2_i.first );
	
	// => {{0}:{1,4}, {3,6}:{5}, {7}:{8}}
	std::map<set<size_t>,set<size_t>> a1s_a2s;
	for( const auto & a2_a1s_i : a2_a1s )
		a1s_a2s[ a2_a1s_i.second ].insert( a2_a1s_i.first );
	a2_a1s.clear();
	
	// => {R(0):R(1)UR(4), R(3)UR(6):R(5), R(7):R(8)}
	// imaging R(0)==R(3)UR(6):   
	// => {R(0):R(1)UR(4)UR(5), R(7):R(8)}
	std::map<set<size_t>,set<size_t>> a1Rs_a2s;
	for( const auto & a1s_a2s_i : a1s_a2s )
	{
		const set<size_t> a1Rs_i = get_set_adjs(a1s_a2s_i.first);
		const set<size_t> a2Rs_i = get_set_adjs(a1s_a2s_i.second);
		a1Rs_a2s[a1Rs_i].insert(a2Rs_i.begin(), a2Rs_i.end());
	}
	a1s_a2s.clear();
	
	return a1Rs_a2s;
}

set<std::pair<size_t,size_t>> Abfs::get_H_pairs_core( const std::vector<std::pair<size_t,size_t>> &atom_pairs )
{
	ModuleBase::TITLE("Exx_Lcao","get_H_pairs_core");
	
	const std::vector<std::map<size_t,std::vector<Abfs::Vector3_Order<int>>>> adjs = Abfs::get_adjs();
	
	auto get_set_adjs = [&adjs]( const set<size_t> &as ) -> set<size_t>
	{
		set<size_t> aRs;
		for( const size_t a : as )
			for( const auto &aR : adjs[a] )
				aRs.insert(aR.first);
		return aRs;
	};
	
	// {(0,1),(3,5),(0,4),(7,8),(6,5)}
	const std::vector<std::pair<size_t,size_t>> & a1_a2 = atom_pairs;
	
	// => {1:{0}, 5:{3,6}, 4:{0}, 8:{7}}
	std::map<size_t,set<size_t>> a2_a1s;
	for( const auto & a1_a2_i : a1_a2 )
		a2_a1s[ a1_a2_i.second ].insert( a1_a2_i.first );
	
	// => {{0}:{1,4}, {3,6}:{5}, {7}:{8}}
	std::map<set<size_t>,set<size_t>> a1s_a2s;
	for( const auto & a2_a1s_i : a2_a1s )
		a1s_a2s[ a2_a1s_i.second ].insert( a2_a1s_i.first );
	a2_a1s.clear();
	
	// => {R(0):{1,4}, R(3)UR(6):{5}, R(7):{8}}
	// imaging R(0)==R(3)UR(6):   
	// => {R(0):{1,4,5}, R(7):{8}}
	std::map<set<size_t>,set<size_t>> a1Rs_a2s;
	for( const auto & a1s_a2s_i : a1s_a2s )
	{
		const set<size_t> a1Rs_i = get_set_adjs(a1s_a2s_i.first);
		a1Rs_a2s[a1Rs_i].insert(a1s_a2s_i.second.begin(), a1s_a2s_i.second.end());
	}
	a1s_a2s.clear();
	
	// => {R(0):R(1)UR(4)UR(5), R(7):R(8)}
	set<std::pair<size_t,size_t>> a1Rs_a2Rs;
	for( const auto & a1Rs_a2s_i : a1Rs_a2s )
	{
		const set<size_t> a2Rs_i = get_set_adjs(a1Rs_a2s_i.second);
		for( const size_t a1R_i : a1Rs_a2s_i.first )
			for( const size_t a2R_i : a2Rs_i )
				a1Rs_a2Rs.insert({a1R_i,a2R_i});
	}
	a1Rs_a2s.clear();
	
	return a1Rs_a2Rs;
}

std::vector<Abfs::Vector3_Order<int>> Abfs::get_Coulomb_potential_boxes( const double rmesh_times )
{
	std::vector<Vector3_Order<int>> Coulomb_potential_boxes;
	const double Rcut = GlobalC::ORB.get_Rmax() * rmesh_times;
	const int nx = std::ceil(Rcut/std::abs((GlobalC::ucell.a2^GlobalC::ucell.a3).normalize()*GlobalC::ucell.a1*GlobalC::ucell.lat0));
	const int ny = std::ceil(Rcut/std::abs((GlobalC::ucell.a1^GlobalC::ucell.a3).normalize()*GlobalC::ucell.a2*GlobalC::ucell.lat0));
	const int nz = std::ceil(Rcut/std::abs((GlobalC::ucell.a1^GlobalC::ucell.a2).normalize()*GlobalC::ucell.a3*GlobalC::ucell.lat0));
	for(int x=-nx; x<=nx; ++x)
		for(int y=-ny; y<=ny; ++y)
			for(int z=-nz; z<=nz; ++z)
				Coulomb_potential_boxes.push_back({x,y,z});
	return Coulomb_potential_boxes;
}

std::vector<Abfs::Vector3_Order<int>> Abfs::get_Born_von_Karmen_boxes( const Abfs::Vector3_Order<int> &Born_von_Karman_period )
{
	std::vector<Abfs::Vector3_Order<int>> Born_von_Karman_boxes;
	for( int ix=0; ix<Born_von_Karman_period.x; ++ix )
		for( int iy=0; iy<Born_von_Karman_period.y; ++iy )
			for( int iz=0; iz<Born_von_Karman_period.z; ++iz )
				Born_von_Karman_boxes.push_back( Abfs::Vector3_Order<int>{ix,iy,iz} % Born_von_Karman_period );
	return Born_von_Karman_boxes;
}


#ifdef __MPI //liyuanbo 2022/2/23
std::shared_ptr<ModuleBase::matrix> Abfs::cal_I( const std::shared_ptr<ModuleBase::matrix> &m )
{
//	ModuleBase::TITLE("Abfs","cal_I1");
	Exx_Abfs::Inverse_Matrix_Double I;
	I.init( m->nc );
	I.input( m );
	
	#if TEST_EXX_LCAO==1
		std::ofstream ofs("I.dat",std::ofstream::app);
		ofs<<"@"<<std::endl<<I.A<<std::endl;
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
	
	I.cal_inverse( Exx_Abfs::Inverse_Matrix_Double::Method::dsyev );

	#if TEST_EXX_LCAO==1
		ofs<<"#"<<std::endl<<I.A<<std::endl;	
		ofs.close();
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif

	std::shared_ptr<ModuleBase::matrix> m_new = make_shared<ModuleBase::matrix>(m->nr,m->nc);
	I.output( m_new );
	return m_new;
}

std::vector<std::vector<std::shared_ptr<ModuleBase::matrix>>> Abfs::cal_I( const std::vector<std::vector<std::shared_ptr<ModuleBase::matrix>>> &ms )
{
//	ModuleBase::TITLE("Abfs","cal_I4");
	Exx_Abfs::Inverse_Matrix_Double I;
	I.init( ms[0][0]->nc + ms[0][1]->nc );
	I.input( ms[0][0], ms[0][1], ms[1][0], ms[1][1] );

	#if TEST_EXX_LCAO==1
		std::ofstream ofs("I.dat",std::ofstream::app);
		ofs<<"@"<<std::endl<<I.A<<std::endl;
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
	
	I.cal_inverse( Exx_Abfs::Inverse_Matrix_Double::Method::dsyev );

	#if TEST_EXX_LCAO==1
		ofs<<"#"<<std::endl<<I.A<<std::endl;	
		ofs.close();
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
	
	std::vector<std::vector<std::shared_ptr<ModuleBase::matrix>>> ms_new( 2, std::vector<std::shared_ptr<ModuleBase::matrix>>(2) );
	for( size_t i=0; i!=2; ++i )
		for(size_t j=0; j!=2; ++j )
			ms_new[i][j] = make_shared<ModuleBase::matrix>( ms[i][j]->nr, ms[i][j]->nc );
	I.output( ms_new[0][0], ms_new[0][1], ms_new[1][0], ms_new[1][1] );
	return ms_new;
}
#endif
