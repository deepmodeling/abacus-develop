#include "abfs.h"
#include "abfs-vector3_order.h"
#include "abfs-template.h"

#include "exx_abfs-inverse_matrix_double.h"
#include "src_pw/global.h"
#include "src_global/global_function.h"
#include <omp.h>

#include <fstream>		// Peize Lin test
#include <iomanip>		// Peize Lin test
#include "src_external/src_test/src_global/matrix-test.h"		// Peize Lin test
#include"src_external/src_test/test_function.h"

map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> Abfs::cal_Cs(
	const set<size_t> &atom_centres,
	const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
	const Exx_Abfs::Matrix_Orbs21 &m_abfslcaos_lcaos,
	const Element_Basis_Index::IndexLNM &index_abfs,
	const Element_Basis_Index::IndexLNM &index_lcaos,
	const double threshold,
	map<size_t,map<size_t,map<Vector3_Order<double>,weak_ptr<matrix>>>> &Cws,
	map<size_t,map<size_t,map<Vector3_Order<double>,weak_ptr<matrix>>>> &Vws )
{
	TITLE("Abfs","cal_Cs");
	pthread_rwlock_t rwlock_Cw;	pthread_rwlock_init(&rwlock_Cw,NULL);
	pthread_rwlock_t rwlock_Vw;	pthread_rwlock_init(&rwlock_Vw,NULL);
	
	const vector<size_t> atom_centres_vector( atom_centres.begin(), atom_centres.end() );
	const vector<map<size_t,vector<Abfs::Vector3_Order<int>>>> adjs = get_adjs();
	
	// pre-cal Vws on same atom, speed up DPcal_V() in DPcal_C()
	vector<shared_ptr<matrix>> Vs_same_atom(ucell.ntype);
	for(size_t it=0; it!=ucell.ntype; ++it)
		Vs_same_atom[it] = DPcal_V( it,it,{0,0,0}, m_abfs_abfs, index_abfs, 0,true, rwlock_Vw,Vws );
	
	const int omp_threads = omp_get_max_threads();
	omp_set_num_threads(std::max(1UL,omp_threads/atom_centres_vector.size()));
	
	map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> Cs;
	#pragma omp parallel for
	for( int i_iat1=0; i_iat1<atom_centres_vector.size(); ++i_iat1 )
	{
		const size_t iat1 = atom_centres_vector[i_iat1];
		const size_t it1 = ucell.iat2it[iat1];
		const size_t ia1 = ucell.iat2ia[iat1];
		const Vector3_Order<double> tau1( ucell.atoms[it1].tau[ia1] );
		
		for( const auto & atom2 : adjs[iat1] )
		{
			const int iat2 = atom2.first;
			const int it2 = ucell.iat2it[iat2];
			const int ia2 = ucell.iat2ia[iat2];
			const Vector3_Order<double> tau2( ucell.atoms[it2].tau[ia2] );
			for( const Vector3<int> &box2 : atom2.second )
			{
//				cout<<"cal_Cs\t"<<iat1<<"\t"<<iat2<<"\t"<<box2<<endl;
				const shared_ptr<matrix> C = DPcal_C( 
					it1, it2, -tau1+tau2+(box2*ucell.latvec), 
					m_abfs_abfs, m_abfslcaos_lcaos, index_abfs, index_lcaos, 
					threshold, true, rwlock_Cw, rwlock_Vw, Cws, Vws );
				#pragma omp critical(Abfs_cal_Cs)
				Cs[iat1][iat2][box2] = C;
			}
		}
	}
	omp_set_num_threads(omp_threads);
	Abfs::delete_threshold_ptrs( Cs, threshold );
	Vs_same_atom.clear();
	Abfs::delete_empty_ptrs( Vws );
	pthread_rwlock_destroy(&rwlock_Cw);
	pthread_rwlock_destroy(&rwlock_Vw);
	return Cs;
}

/*map<size_t,map<size_t,shared_ptr<matrix>>> 
	Abfs::cal_Vps(
		const set<pair<size_t,size_t>> &atom_pairs,
		const vector<Vector3_Order<int>> &Coulomb_potential_boxes,
		const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
		const Element_Basis_Index::Index &index_abfs,
		map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> &Vs,
		map<size_t,map<size_t,map<Vector3_Order<double>,weak_ptr<matrix>>>> &Vws )
{
	map<size_t,map<size_t,shared_ptr<matrix>>> Vps;
	
	for( const pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = ucell.iat2it[iat1];
		const size_t ia1 = ucell.iat2ia[iat1];
		const size_t it2 = ucell.iat2it[iat2];
		const size_t ia2 = ucell.iat2ia[iat2];
		const Vector3_Order<double> tau1 = ucell.atoms[it1].tau[ia1];
		const Vector3_Order<double> tau2 = ucell.atoms[it2].tau[ia2];
		
		for( const Vector3_Order<int> &box2 : Coulomb_potential_boxes )
		{
			const vector<shared_ptr<matrix>> Vs_tmp = DPcal_V( it1, it2, -tau1+tau2+(box2*ucell.latvec), m_abfs_abfs, index_abfs, Vws );
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
map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> 
	Abfs::cal_Vps(
		const vector<pair<size_t,size_t>> &atom_pairs,
		const vector<Vector3_Order<int>> &Coulomb_potential_boxes,
		const vector<Vector3_Order<int>> &Born_von_Karman_boxes,
		const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
		const Element_Basis_Index::IndexLNM &index_abfs,
		map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> &Vs,
		map<size_t,map<size_t,map<Vector3_Order<double>,weak_ptr<matrix>>>> &Vws )
{
	TITLE("Abfs","cal_Vps");	
	
	map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> Vps;
	
	for( const pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = ucell.iat2it[iat1];
		const size_t ia1 = ucell.iat2ia[iat1];
		const size_t it2 = ucell.iat2it[iat2];
		const size_t ia2 = ucell.iat2ia[iat2];
		const Vector3_Order<double> tau1 = ucell.atoms[it1].tau[ia1];
		const Vector3_Order<double> tau2 = ucell.atoms[it2].tau[ia2];
		
		for( const Vector3_Order<int> &box2 : Coulomb_potential_boxes )
		{
cout<<"cal_Vs\t"<<iat1<<"\t"<<iat2<<"\t"<<box2<<endl;
			const vector<shared_ptr<matrix>> Vs_tmp = DPcal_V( it1, it2, -tau1+tau2+(box2*ucell.latvec), m_abfs_abfs, index_abfs, Vws );
			Vs[iat1][iat2][box2] = Vs_tmp[0];	Vs[iat2][iat1][-box2] = Vs_tmp[1];
		}
		
		vector<ComplexMatrix> Vkps( kv.nks/NSPIN, {Vs[iat1][iat2][{0,0,0}]->nr, Vs[iat1][iat2][{0,0,0}]->nc} );
		for( size_t ik=0; ik!=kv.nks/NSPIN; ++ik )
			for( const Vector3_Order<int> &box2 : Coulomb_potential_boxes )
				Vkps[ik] += static_cast<ComplexMatrix>(*Vs[iat1][iat2][box2]) * exp(-TWO_PI*IMAG_UNIT*(kv.kvec_c[ik]*(box2*ucell.latvec)));
		
		for( const Vector3_Order<int> &box2 : Born_von_Karman_boxes )
		{
			ComplexMatrix Vps_tmp ( Vkps[0].nr, Vkps[0].nc );
			for( size_t ik=0; ik!=kv.nks/NSPIN; ++ik )
				Vps_tmp += Vkps[ik] * kv.wk[ik]*(0.5*NSPIN) * exp(TWO_PI*IMAG_UNIT*(kv.kvec_c[ik]*(box2*ucell.latvec)));
			Vps[iat1][iat2][box2] = make_shared<matrix>( Vps_tmp.real() );
		}
	}
	return Vps;
}
*/

map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> Abfs::cal_Vs(
	const vector<pair<size_t,size_t>> &atom_pairs,
	const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
	const Element_Basis_Index::IndexLNM &index_abfs,
	const double rmesh_times,
	const double threshold,
	map<size_t,map<size_t,map<Vector3_Order<double>,weak_ptr<matrix>>>> &Vws )
{
	TITLE("Abfs","cal_Vs");
	pthread_rwlock_t rwlock_Vw;	pthread_rwlock_init(&rwlock_Vw,NULL);
	vector<Abfs::Vector3_Order<int>> Coulomb_potential_boxes = get_Coulomb_potential_boxes(rmesh_times);

	map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> Vs;
	#pragma omp parallel for
	for( int i_atom_pair=0; i_atom_pair<atom_pairs.size(); ++i_atom_pair )
	{
		const size_t iat1 = atom_pairs[i_atom_pair].first;
		const size_t iat2 = atom_pairs[i_atom_pair].second;
		const size_t it1 = ucell.iat2it[iat1];
		const size_t ia1 = ucell.iat2ia[iat1];
		const size_t it2 = ucell.iat2it[iat2];
		const size_t ia2 = ucell.iat2ia[iat2];
		const Vector3_Order<double> tau1 = ucell.atoms[it1].tau[ia1];
		const Vector3_Order<double> tau2 = ucell.atoms[it2].tau[ia2];
		const double Rcut = std::min( ORB.Phi[it1].getRcut()*rmesh_times+ORB.Phi[it2].getRcut(), ORB.Phi[it1].getRcut()+ORB.Phi[it2].getRcut()*rmesh_times );
		
		for( const Vector3_Order<int> &box2 : Coulomb_potential_boxes )
		{
			const Vector3_Order<double> delta_R = -tau1+tau2+(box2*ucell.latvec);
			if( delta_R.norm()*ucell.lat0 < Rcut )
			{
//				cout<<"cal_Vs\t"<<iat1<<"\t"<<iat2<<"\t"<<box2<<"\t"<<delta_R<<"\t"<<delta_R.norm()<<"\t"<<delta_R.norm()*ucell.lat0<<"\t"<<ORB.Phi[it1].getRcut()*rmesh_times+ORB.Phi[it2].getRcut()<<endl;
				const shared_ptr<matrix> V = DPcal_V( 
					it1, it2, delta_R, 
					m_abfs_abfs, index_abfs, 
					threshold, true, rwlock_Vw, Vws );
				#pragma omp critical(Abfs_cal_Vs)
				Vs[iat1][iat2][box2] = V;
			}
		}
	}
	Abfs::delete_threshold_ptrs( Vs, threshold );
	pthread_rwlock_destroy(&rwlock_Vw);
	return Vs;
}

map<Abfs::Vector3_Order<int>,shared_ptr<matrix>> Abfs::cal_mps(
	const Abfs::Vector3_Order<int> &Born_von_Karman_period,
	const map<Vector3_Order<int>,shared_ptr<matrix>> &ms )
{
	map< Vector3_Order<int>, vector<Vector3_Order<int>> > indexs;
	for( const auto & m : ms )
	{
		const Vector3_Order<int> & box = m.first;
		indexs[box%Born_von_Karman_period].push_back(box);
	}
	
	#if TEST_EXX_LCAO==1
	{
		static int istep=0;
		ofstream ofs( "mps_index_"+TO_STRING(iat1)+"_"+TO_STRING(iat2)+"_"+TO_STRING(istep++) );
		for( const auto index : indexs )
		{
			ofs<<index.first<<endl;
			for( const Vector3_Order<int> & box : index.second )
				ofs<<"\t"<<box<<endl;
		}
		ofs.close();
	}
	#elif TEST_EXX_LCAO==-1
		#error
	#endif
	
	map<Vector3_Order<int>,shared_ptr<matrix>> mps;
	for( const auto & index : indexs )
	{
		const Vector3_Order<int> & boxp = index.first;
		assert(index.second.size()>=1);
		if( index.second.size()==1 )
			mps[boxp] = ms.at(index.second[0]);
		else
		{
			mps[boxp] = make_shared<matrix>();
			*mps[boxp] = *ms.at(index.second[0]);
			for( size_t i=1; i<index.second.size(); ++i )
				*mps[boxp] += *ms.at(index.second[i]);
		}
	}
	return mps;
}

map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> Abfs::cal_mps(
	const Abfs::Vector3_Order<int> &Born_von_Karman_period,
	const map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> &ms )
{
	TITLE("Abfs","cal_mps");
	map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> mps;
	for( const auto & m1s : ms )
		for( const auto & m12s : m1s.second )
			mps[m1s.first][m12s.first] = cal_mps( Born_von_Karman_period, m12s.second );
	return mps;
}


shared_ptr<matrix> Abfs::DPcal_C( 
	const size_t &it1, 
	const size_t &it2, 
	const Vector3_Order<double> &R, 
	const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
	const Exx_Abfs::Matrix_Orbs21 &m_abfslcaos_lcaos,
	const Element_Basis_Index::IndexLNM &index_abfs,
	const Element_Basis_Index::IndexLNM &index_lcaos,
	const double threshold,
	const bool writable,
	pthread_rwlock_t &rwlock_Cw,
	pthread_rwlock_t &rwlock_Vw,
	map<size_t,map<size_t,map<Vector3_Order<double>,weak_ptr<matrix>>>> &Cws,
	map<size_t,map<size_t,map<Vector3_Order<double>,weak_ptr<matrix>>>> &Vws )
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
	
//	TITLE("Abfs","DPcal_C");
	pthread_rwlock_rdlock(&rwlock_Cw);
	const weak_ptr<matrix> * const Cws_ptr   = static_cast<const weak_ptr<matrix> * const>( MAP_EXIST( Cws, it1, it2, R ) );
	pthread_rwlock_unlock(&rwlock_Cw);
	
	if( Cws_ptr && !Cws_ptr->expired() )
		return Cws_ptr->lock();
	else
	{
//		cout<<"DPcal_C\t"<<it1<<"\t"<<it2<<"\t"<<R<<endl;
		if( (Vector3<double>(0,0,0)==R) && (it1==it2) )
		{
			const shared_ptr<matrix> A = 
				make_shared<matrix>( m_abfslcaos_lcaos.cal_overlap_matrix( it1,it2,0,0,index_abfs,index_lcaos,index_lcaos,Exx_Abfs::Matrix_Orbs21::Matrix_Order::A2B_A1 ) );
			const shared_ptr<matrix> V = DPcal_V(it1,it2,{0,0,0}, m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
			const shared_ptr<matrix> L = cal_I(V);
			shared_ptr<matrix> C = make_shared<matrix>( 0.5 * *A * *L );
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
			const vector<shared_ptr<matrix>> A = 
				{ make_shared<matrix>( m_abfslcaos_lcaos.cal_overlap_matrix(it1,it2,0,R ,index_abfs,index_lcaos,index_lcaos,Exx_Abfs::Matrix_Orbs21::Matrix_Order::A2B_A1) ) ,
				  make_shared<matrix>( m_abfslcaos_lcaos.cal_overlap_matrix(it2,it1,0,-R,index_abfs,index_lcaos,index_lcaos,Exx_Abfs::Matrix_Orbs21::Matrix_Order::BA2_A1) ) };
			
			const shared_ptr<matrix> V_00 = DPcal_V(it1,it1,{0,0,0}, m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
			const shared_ptr<matrix> V_01 = DPcal_V(it1,it2,R,       m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
			const shared_ptr<matrix> V_10 = DPcal_V(it2,it1,-R,      m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
			const shared_ptr<matrix> V_11 = DPcal_V(it2,it2,{0,0,0}, m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
			const vector<vector<shared_ptr<matrix>>> V = 
				{{ V_00, V_01 },
				 { V_10, V_11 }};

			const vector<vector<shared_ptr<matrix>>> L = cal_I(V);

			shared_ptr<matrix> C = make_shared<matrix>( *A[0] * *L[0][0] + *A[1] * *L[1][0] );
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

shared_ptr<matrix> Abfs::DPcal_V( 
	const size_t &it1, 
	const size_t &it2, 
	const Vector3_Order<double> &R, 
	const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
	const Element_Basis_Index::IndexLNM &index_abfs,
	const double threshold,
	const bool writable,
	pthread_rwlock_t &rwlock_Vw,
	map<size_t,map<size_t,map<Vector3_Order<double>,weak_ptr<matrix>>>> &Vws)
{
//	TITLE("Abfs","DPcal_V");
	pthread_rwlock_rdlock(&rwlock_Vw);
	const weak_ptr<matrix> * const Vws12_ptr = static_cast<const weak_ptr<matrix> * const>( MAP_EXIST( Vws, it1, it2, R ) );
	const weak_ptr<matrix> * const Vws21_ptr = static_cast<const weak_ptr<matrix> * const>( MAP_EXIST( Vws, it2, it1, -R ) );
	pthread_rwlock_unlock(&rwlock_Vw);
	
	if( Vws12_ptr && !Vws12_ptr->expired() )
		return Vws12_ptr->lock();
	else if( Vws21_ptr && !Vws21_ptr->expired() )
	{
		shared_ptr<matrix> VT = make_shared<matrix>( transpose(*Vws21_ptr->lock()) );
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
//		cout<<"DPcal_V\t"<<it1<<"\t"<<it2<<"\t"<<R<<endl;
		shared_ptr<matrix> V = make_shared<matrix>( m_abfs_abfs.cal_overlap_matrix(it1,it2,0,R,index_abfs,index_abfs) );
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

map<size_t,vector<Abfs::Vector3_Order<int>>> Abfs::get_adjs( const size_t &iat )
{
//	TITLE("Abfs","get_adjs");
	const int it = ucell.iat2it[iat];
	const int ia = ucell.iat2ia[iat];
	const Vector3<double> &tau = ucell.atoms[it].tau[ia];
	
	map<size_t,vector<Vector3_Order<int>>> adjs;
	GridD.Find_atom( tau, it, ia );
	for( int ad=0; ad<GridD.getAdjacentNum()+1; ++ad )
	{
		const size_t it_ad = GridD.getType(ad);
		const size_t ia_ad = GridD.getNatom(ad);
		const int iat_ad = ucell.itia2iat(it_ad,ia_ad);
		const Vector3<int> box_ad = GridD.getBox(ad);
		const Vector3<double> tau_ad = GridD.getAdjacentTau(ad);
		
		if( (tau-tau_ad).norm()*ucell.lat0 < ORB.Phi[it].getRcut()+ORB.Phi[it_ad].getRcut() )
			adjs[iat_ad].push_back(box_ad);
	}
	return adjs;
}

vector<map<size_t,vector<Abfs::Vector3_Order<int>>>> Abfs::get_adjs()
{
	vector<map<size_t,vector<Abfs::Vector3_Order<int>>>> adjs(ucell.nat);
	for( size_t iat=0; iat!=ucell.nat; ++iat )
		adjs[iat] = Abfs::get_adjs(iat);
	return adjs;
}

/*
set<pair<size_t,size_t>> Abfs::get_H_pairs_core( const vector<pair<size_t,size_t>> &atom_pairs )
{
	TITLE("Exx_Lcao","allocate_Hexx");

	set<pair<size_t,size_t>> H_atom_pairs_core;
	for( const pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;

		const map<size_t,vector<Abfs::Vector3_Order<int>>> adj1s = Abfs::get_adjs(iat1);
		const map<size_t,vector<Abfs::Vector3_Order<int>>> adj2s = Abfs::get_adjs(iat2);

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

map<set<size_t>,set<size_t>> Abfs::get_H_pairs_core_group( const vector<pair<size_t,size_t>> &atom_pairs )
{
	TITLE("Abfs","get_H_pairs_core_group");
	
	const vector<map<size_t,vector<Abfs::Vector3_Order<int>>>> adjs = Abfs::get_adjs();
	
	auto get_set_adjs = [&adjs]( const set<size_t> &as ) -> set<size_t>
	{
		set<size_t> aRs;
		for( const size_t a : as )
			for( const auto &aR : adjs[a] )
				aRs.insert(aR.first);
		return aRs;
	};
	
	// {(0,1),(3,5),(0,4),(7,8),(6,5)}
	const vector<pair<size_t,size_t>> & a1_a2 = atom_pairs;
	
	// => {1:{0}, 5:{3,6}, 4:{0}, 8:{7}}
	map<size_t,set<size_t>> a2_a1s;
	for( const auto & a1_a2_i : a1_a2 )
		a2_a1s[ a1_a2_i.second ].insert( a1_a2_i.first );
	
	// => {{0}:{1,4}, {3,6}:{5}, {7}:{8}}
	map<set<size_t>,set<size_t>> a1s_a2s;
	for( const auto & a2_a1s_i : a2_a1s )
		a1s_a2s[ a2_a1s_i.second ].insert( a2_a1s_i.first );
	a2_a1s.clear();
	
	// => {R(0):R(1)UR(4), R(3)UR(6):R(5), R(7):R(8)}
	// imaging R(0)==R(3)UR(6):   
	// => {R(0):R(1)UR(4)UR(5), R(7):R(8)}
	map<set<size_t>,set<size_t>> a1Rs_a2s;
	for( const auto & a1s_a2s_i : a1s_a2s )
	{
		const set<size_t> a1Rs_i = get_set_adjs(a1s_a2s_i.first);
		const set<size_t> a2Rs_i = get_set_adjs(a1s_a2s_i.second);
		a1Rs_a2s[a1Rs_i].insert(a2Rs_i.begin(), a2Rs_i.end());
	}
	a1s_a2s.clear();
	
	return a1Rs_a2s;
}

set<pair<size_t,size_t>> Abfs::get_H_pairs_core( const vector<pair<size_t,size_t>> &atom_pairs )
{
	TITLE("Exx_Lcao","get_H_pairs_core");
	
	const vector<map<size_t,vector<Abfs::Vector3_Order<int>>>> adjs = Abfs::get_adjs();
	
	auto get_set_adjs = [&adjs]( const set<size_t> &as ) -> set<size_t>
	{
		set<size_t> aRs;
		for( const size_t a : as )
			for( const auto &aR : adjs[a] )
				aRs.insert(aR.first);
		return aRs;
	};
	
	// {(0,1),(3,5),(0,4),(7,8),(6,5)}
	const vector<pair<size_t,size_t>> & a1_a2 = atom_pairs;
	
	// => {1:{0}, 5:{3,6}, 4:{0}, 8:{7}}
	map<size_t,set<size_t>> a2_a1s;
	for( const auto & a1_a2_i : a1_a2 )
		a2_a1s[ a1_a2_i.second ].insert( a1_a2_i.first );
	
	// => {{0}:{1,4}, {3,6}:{5}, {7}:{8}}
	map<set<size_t>,set<size_t>> a1s_a2s;
	for( const auto & a2_a1s_i : a2_a1s )
		a1s_a2s[ a2_a1s_i.second ].insert( a2_a1s_i.first );
	a2_a1s.clear();
	
	// => {R(0):{1,4}, R(3)UR(6):{5}, R(7):{8}}
	// imaging R(0)==R(3)UR(6):   
	// => {R(0):{1,4,5}, R(7):{8}}
	map<set<size_t>,set<size_t>> a1Rs_a2s;
	for( const auto & a1s_a2s_i : a1s_a2s )
	{
		const set<size_t> a1Rs_i = get_set_adjs(a1s_a2s_i.first);
		a1Rs_a2s[a1Rs_i].insert(a1s_a2s_i.second.begin(), a1s_a2s_i.second.end());
	}
	a1s_a2s.clear();
	
	// => {R(0):R(1)UR(4)UR(5), R(7):R(8)}
	set<pair<size_t,size_t>> a1Rs_a2Rs;
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

vector<Abfs::Vector3_Order<int>> Abfs::get_Coulomb_potential_boxes( const double rmesh_times )
{
	vector<Vector3_Order<int>> Coulomb_potential_boxes;
	const double Rcut = ORB.get_Rmax() * rmesh_times;
	const int nx = std::ceil(Rcut/std::abs((ucell.a2^ucell.a3).normalize()*ucell.a1*ucell.lat0));
	const int ny = std::ceil(Rcut/std::abs((ucell.a1^ucell.a3).normalize()*ucell.a2*ucell.lat0));
	const int nz = std::ceil(Rcut/std::abs((ucell.a1^ucell.a2).normalize()*ucell.a3*ucell.lat0));
	for(int x=-nx; x<=nx; ++x)
		for(int y=-ny; y<=ny; ++y)
			for(int z=-nz; z<=nz; ++z)
				Coulomb_potential_boxes.push_back({x,y,z});
	return Coulomb_potential_boxes;
}

vector<Abfs::Vector3_Order<int>> Abfs::get_Born_von_Karmen_boxes( const Abfs::Vector3_Order<int> &Born_von_Karman_period )
{
	vector<Abfs::Vector3_Order<int>> Born_von_Karman_boxes;
	for( int ix=0; ix<Born_von_Karman_period.x; ++ix )
		for( int iy=0; iy<Born_von_Karman_period.y; ++iy )
			for( int iz=0; iz<Born_von_Karman_period.z; ++iz )
				Born_von_Karman_boxes.push_back( Abfs::Vector3_Order<int>{ix,iy,iz} % Born_von_Karman_period );
	return Born_von_Karman_boxes;
}

shared_ptr<matrix> Abfs::cal_I( const shared_ptr<matrix> &m )
{
//	TITLE("Abfs","cal_I1");
	Exx_Abfs::Inverse_Matrix_Double I;
	I.init( m->nc );
	I.input( m );
	
	#if TEST_EXX_LCAO==1
		ofstream ofs("I.dat",ofstream::app);
		ofs<<"@"<<endl<<I.A<<endl;
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
	
	I.cal_inverse( Exx_Abfs::Inverse_Matrix_Double::Method::dsyev );

	#if TEST_EXX_LCAO==1
		ofs<<"#"<<endl<<I.A<<endl;	
		ofs.close();
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif

	shared_ptr<matrix> m_new = make_shared<matrix>(m->nr,m->nc);
	I.output( m_new );
	return m_new;
}

vector<vector<shared_ptr<matrix>>> Abfs::cal_I( const vector<vector<shared_ptr<matrix>>> &ms )
{
//	TITLE("Abfs","cal_I4");
	Exx_Abfs::Inverse_Matrix_Double I;
	I.init( ms[0][0]->nc + ms[0][1]->nc );
	I.input( ms[0][0], ms[0][1], ms[1][0], ms[1][1] );

	#if TEST_EXX_LCAO==1
		ofstream ofs("I.dat",ofstream::app);
		ofs<<"@"<<endl<<I.A<<endl;
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
	
	I.cal_inverse( Exx_Abfs::Inverse_Matrix_Double::Method::dsyev );

	#if TEST_EXX_LCAO==1
		ofs<<"#"<<endl<<I.A<<endl;	
		ofs.close();
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
	
	vector<vector<shared_ptr<matrix>>> ms_new( 2, vector<shared_ptr<matrix>>(2) );
	for( size_t i=0; i!=2; ++i )
		for(size_t j=0; j!=2; ++j )
			ms_new[i][j] = make_shared<matrix>( ms[i][j]->nr, ms[i][j]->nc );
	I.output( ms_new[0][0], ms_new[0][1], ms_new[1][0], ms_new[1][1] );
	return ms_new;
}
