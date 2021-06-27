#include "exx_abfs-pca.h"

#include "exx_abfs-abfs_index.h"
#include "exx_abfs-matrix_orbs11.h"
#include "exx_abfs-matrix_orbs21.h"
#include "exx_abfs-inverse_matrix_double.h"

#include "module_base/lapack_connector.h"
#include "../module_base/global_function.h"

#include <cassert>
#include <limits>

#include "src_external/src_test/src_global/element_basis_index-test.h"		// Peize Lin test
#include "src_external/src_test/src_ri/exx_lcao-test.h"		// Peize Lin test
#include <sys/time.h>			// Peize Lin test
#include "src_lcao/global_fp.h"		// Peize Lin test

vector<vector<pair<vector<double>,matrix>>> Exx_Abfs::PCA::cal_PCA( 
	const vector<vector<vector<Numerical_Orbital_Lm>>> &lcaos, 
	const vector<vector<vector<Numerical_Orbital_Lm>>> &abfs,
	const double kmesh_times )
{
	TITLE("Exx_Abfs::PCA::cal_PCA");
ofstream ofs(exx_lcao.test_dir.process+"time_"+TO_STRING(MY_RANK),ofstream::app);
timeval t_start;
	
	const Element_Basis_Index::Range
		&& range_lcaos = Exx_Abfs::Abfs_Index::construct_range( lcaos );
	const Element_Basis_Index::IndexLNM
		&& index_lcaos = Element_Basis_Index::construct_index( range_lcaos );

	const Element_Basis_Index::Range
		&& range_abfs = Exx_Abfs::Abfs_Index::construct_range( abfs );
	const Element_Basis_Index::IndexLNM
		&& index_abfs = Element_Basis_Index::construct_index( range_abfs );
		
ofs<<range_lcaos<<endl;
ofs<<range_abfs<<endl;

	const int Lmax_bak = Exx_Abfs::Lmax;
	Exx_Abfs::Lmax = std::numeric_limits<int>::min();
	for( size_t T=0; T!=abfs.size(); ++T )
		Exx_Abfs::Lmax = std::max( Exx_Abfs::Lmax, static_cast<int>(abfs[T].size())-1 );

	Exx_Abfs::Matrix_Orbs21 m_abfslcaos_lcaos;
//gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init( 1, kmesh_times, 1 );
//ofs<<"TIME@m_abfslcaos_lcaos.init\t"<<time_during(t_start)<<endl;
//gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init_radial( abfs, lcaos, lcaos );
//ofs<<"TIME@m_abfslcaos_lcaos.init_radial\t"<<time_during(t_start)<<endl;

	map<size_t,map<size_t,set<double>>> delta_R;
	for( size_t it=0; it!=abfs.size(); ++it )
		delta_R[it][it] = {0.0};
//gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init_radial_table(delta_R);
//ofs<<"TIME@m_abfslcaos_lcaos.init_radial_table\t"<<time_during(t_start)<<endl;

	Exx_Abfs::Lmax = Lmax_bak;
	
	vector<vector<pair<vector<double>,matrix>>> eig(abfs.size());
	for( size_t T=0; T!=abfs.size(); ++T )
	{
		const matrix && A = m_abfslcaos_lcaos.cal_overlap_matrix(  
			T, 
			T, 
			Vector3<double>{0,0,0},
			Vector3<double>{0,0,0},  
			index_abfs, 
			index_lcaos,
			index_lcaos,
			Exx_Abfs::Matrix_Orbs21::Matrix_Order::A2B_A1);
//ofs<<"A:"<<endl<<A<<endl;
		
		eig[T].resize(abfs[T].size());
		for( size_t L=0; L!=abfs[T].size(); ++L )
		{
//ofs<<"get_sub_matrix:"<<endl<<get_sub_matrix( A, T, L, range_abfs, index_abfs )<<endl;
//			const matrix A_sub = get_column_mean0_matrix( get_sub_matrix( A, T, L, range_abfs, index_abfs ) );
			const matrix A_sub = get_sub_matrix( A, T, L, range_abfs, index_abfs );
//ofs<<"A_sub:"<<endl<<A_sub<<endl;
//ofs<<"transpose:"<<endl<<transpose(A_sub)<<endl;
//ofs<<"mul:"<<endl<<transpose(A_sub) * A_sub<<endl;
			matrix mm = transpose(A_sub) * A_sub;
//ofs<<"mm:"<<endl<<mm<<endl;
			vector<double> eig_value(mm.nr);
			
			int info;
gettimeofday( &t_start, NULL);
			LapackConnector::dsyev( 'V', 'U', mm, VECTOR_TO_PTR(eig_value), info );
ofs<<"TIME@LapackConnector::dsyev\t"<<time_during(t_start)<<endl;
			if( info )
			{
				cout<<endl<<"info_dsyev = "<<info<<endl;
				ofs_warning<<mm<<endl;
				cout<<"in file "<<__FILE__<<" line "<<__LINE__<<endl;
				QUIT();
			}
			eig[T][L] = make_pair( eig_value, mm );
//for( const auto v : eig_value )
//	ofs<<v<<"\t";
//ofs<<endl;
//ofs<<"mm:"<<endl<<mm<<endl;
		}
	}
ofs.close();
	
	return eig;
}

matrix Exx_Abfs::PCA::get_sub_matrix( 
	const matrix & m,
	const size_t & T,
	const size_t & L,
	const Element_Basis_Index::Range & range,
	const Element_Basis_Index::IndexLNM & index )
{
	TITLE("Exx_Abfs::PCA::get_sub_matrix");
	
	matrix m_sub( m.nr, range[T][L].N );
	for( size_t ir=0; ir!=m.nr; ++ir )
		for( size_t N=0; N!=range[T][L].N; ++N )
			m_sub( ir, N ) = m( ir, index[T][L][N][0] );
	return m_sub;
}


matrix Exx_Abfs::PCA::get_column_mean0_matrix( const matrix & m )
{
	matrix m_new( m.nr, m.nc );
	for( size_t ic=0; ic!=m.nc; ++ic )
	{
		double sum=0;
		for( size_t ir=0; ir!=m.nr; ++ir )
			sum += m(ir,ic);
		const double mean = sum/m.nr;
//ofs<<mean<<"\t";
		for( size_t ir=0; ir!=m.nr; ++ir )
			m_new(ir,ic) = m(ir,ic) - mean;
	}
//ofs<<endl;
	return m_new;
}
