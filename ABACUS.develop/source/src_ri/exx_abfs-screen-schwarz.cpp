#include "exx_abfs-screen-schwarz.h"
#include "src_pw/global.h"
#include "src_global/lapack_connector.h"
#include <cassert>
#include <limits>
#include <thread>

#include "src_external/src_test/src_ri/exx_abfs-screen-test.h"
#include "src_external/src_test/src_ri/exx_lcao-test.h"
//double Exx_Abfs::Screen::Schwarz::num_screen = 0;
//double Exx_Abfs::Screen::Schwarz::num_cal = 0;
	
void Exx_Abfs::Screen::Schwarz::init(
	const bool flag_screen_schwarz_in,
	const double threshold_in )
{
	flag_screen_schwarz = flag_screen_schwarz_in;
	threshold = threshold_in * threshold_in;
}


void Exx_Abfs::Screen::Schwarz::cal_max_pair_fock(
	const set<size_t> &atom_centres,
	const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
	const Exx_Abfs::Matrix_Orbs21 &m_abfslcaos_lcaos,
	const Element_Basis_Index::IndexLNM &index_abfs,
	const Element_Basis_Index::IndexLNM &index_lcaos,
	const Abfs::Vector3_Order<int> &Born_von_Karman_period,
	map<size_t,map<size_t,map<Abfs::Vector3_Order<double>,weak_ptr<matrix>>>> &Cws,
	map<size_t,map<size_t,map<Abfs::Vector3_Order<double>,weak_ptr<matrix>>>> &Vws )
{
	if(!flag_screen_schwarz)	return;
	TITLE("Exx_Abfs::Screen::Schwarz::cal_max_pair_fock");
	pthread_rwlock_t rwlock_Cw;	pthread_rwlock_init(&rwlock_Cw,NULL);
	pthread_rwlock_t rwlock_Vw;	pthread_rwlock_init(&rwlock_Vw,NULL);
	
	// pre-cal Vws on same atom, speed up DPcal_V()
	vector<shared_ptr<matrix>> Vs_same_atom(ucell.ntype);
	for(size_t it=0; it!=ucell.ntype; ++it)
		Vs_same_atom[it] = Abfs::DPcal_V( it,it,{0,0,0}, m_abfs_abfs, index_abfs, 0,true, rwlock_Vw,Vws );	
	
	// m_out( i1, i2, i3 ) = m_in( i2, i1, i3 )
	auto change_matrix_order =[]( const matrix &m_in, const size_t n1 ) -> matrix
	{
		assert( m_in.nr%n1 == 0 );
		matrix m_out( m_in.nr, m_in.nc, false );
		const size_t n2  = m_in.nr/n1, 
					 n3  = m_in.nc, 
					 n13 = n1*n3;
		const auto length = n3*sizeof(double);
		double * m_out_ptr = m_out.c;
		for( size_t i1=0; i1!=n1; ++i1 )
		{
			const double * m_in_ptr = m_in.c + i1*n3;
			for( size_t i2=0; i2!=n2; ++i2, m_in_ptr+=n13, m_out_ptr+=n3 )
				memcpy( m_out_ptr, m_in_ptr, length );
		}
		return m_out;
	};
	
	// max_i m(i,i)
	auto max_diagnal = []( const matrix & m ) -> double
	{
		assert(m.nr==m.nc);
		double diagnal = -std::numeric_limits<double>::max();
		for( size_t i=0; i!=m.nr; ++i )
			diagnal = max( diagnal, m(i,i) );
		return diagnal;
	};
	
	// m * m.T
	auto m_mT = []( const matrix & m1, const matrix & m2 ) -> matrix
	{
		assert( m1.nc == m2.nc );
		matrix mm( m1.nr, m2.nr, false );
		LapackConnector::gemm( 'N','T', m1.nr,m2.nr,m1.nc, 1, m1.c,m1.nc, m2.c,m2.nc, 0, mm.c,mm.nc );
		return mm;
	};

	map<size_t,map<size_t,map<Abfs::Vector3_Order<double>,double>>> max_pair_fock_DP;
	for(const size_t iat1 : atom_centres)
	{	
		const size_t it1 = ucell.iat2it[iat1];
		const size_t ia1 = ucell.iat2ia[iat1];
		const Abfs::Vector3_Order<double> tau1( ucell.atoms[it1].tau[ia1] );
		
		const map<size_t,vector<Abfs::Vector3_Order<int>>> adj = Abfs::get_adjs(iat1);
		for( const auto & atom2 : adj )
		{
			const int iat2 = atom2.first;
			const int it2 = ucell.iat2it[iat2];
			const int ia2 = ucell.iat2ia[iat2];
			const Abfs::Vector3_Order<double> tau2( ucell.atoms[it2].tau[ia2] );
			
			map<Abfs::Vector3_Order<int>,shared_ptr<matrix>> pair_fock_s;
			for( const Abfs::Vector3_Order<int> &box2 : atom2.second )
			{
				const Abfs::Vector3_Order<int> box2p = box2%Born_von_Karman_period;
				if(const double*max_pair_fock_ptr=static_cast<const double*>(MAP_EXIST(max_pair_fock_DP,it1,it2,-tau1+tau2+box2p*ucell.latvec)))
					max_pair_fock[iat1][iat2][box2p] = *max_pair_fock_ptr;
				else
				{
					const Abfs::Vector3_Order<double> R = -tau1+tau2+box2*ucell.latvec;
					const matrix C_12 = *Abfs::DPcal_C( it1,it2,R,  m_abfs_abfs,m_abfslcaos_lcaos, index_abfs,index_lcaos, 0,false, rwlock_Cw,rwlock_Vw,Cws,Vws );
					const matrix C_21 = change_matrix_order( 
						*Abfs::DPcal_C( it2,it1,-R, m_abfs_abfs,m_abfslcaos_lcaos, index_abfs,index_lcaos, 0,false, rwlock_Cw,rwlock_Vw,Cws,Vws ),
						index_lcaos[it1].count_size);
					
					const matrix V_11 = *Abfs::DPcal_V(it1,it1,{0,0,0}, m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
					const matrix V_12 = *Abfs::DPcal_V(it1,it2,R,       m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
					const matrix V_21 = *Abfs::DPcal_V(it2,it1,-R,      m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
					const matrix V_22 = *Abfs::DPcal_V(it2,it2,{0,0,0}, m_abfs_abfs,index_abfs, 0,false, rwlock_Vw,Vws);
					
					pair_fock_s[box2] = make_shared<matrix>(
						  m_mT( C_12 * V_11, C_12 )
						+ m_mT( C_12 * V_12, C_21 )
						+ m_mT( C_21 * V_21, C_12 )
						+ m_mT( C_21 * V_22, C_21 ));
				}
			}
			
			const map<Abfs::Vector3_Order<int>,shared_ptr<matrix>> pair_fock_ps = Abfs::cal_mps( Born_von_Karman_period, pair_fock_s );

			for( const auto & pair_fock_p : pair_fock_ps )
			{
				const Abfs::Vector3_Order<int> & box2p = pair_fock_p.first;
				max_pair_fock[iat1][iat2][box2p] = max_pair_fock_DP[it1][it2][-tau1+tau2+box2p*ucell.latvec] = max_diagnal(*pair_fock_p.second);
			}
		}
	}
	
	Vs_same_atom.clear();
	Abfs::delete_empty_ptrs( Vws );
	pthread_rwlock_destroy(&rwlock_Cw);
	pthread_rwlock_destroy(&rwlock_Vw);
		
//test_screen("schwarz-max_pair_fock.dat",max_pair_fock);
}
