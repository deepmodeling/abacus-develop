#include "exx_abfs-pca.h"

#include "exx_abfs-abfs_index.h"
#include "exx_abfs-matrix_orbs11.h"
#include "exx_abfs-matrix_orbs21.h"
#include "exx_abfs-inverse_matrix_double.h"

#include "../module_base/lapack_connector.h"
#include "../module_base/global_function.h"

#include <cassert>
#include <limits>

#include "../src_ri/test_code/element_basis_index-test.h"		// Peize Lin test
#include "../src_ri/test_code/exx_lcao-test.h"		// Peize Lin test
#include <sys/time.h>			// Peize Lin test
#include "../module_hamilt_lcao/hamilt_lcaodft/global_fp.h"		// Peize Lin test

std::vector<std::vector<std::pair<std::vector<double>,ModuleBase::matrix>>> Exx_Abfs::PCA::cal_PCA( 
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &lcaos, 
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs,
	const double kmesh_times )
{
	ModuleBase::TITLE("Exx_Abfs::PCA::cal_PCA");
//std::ofstream ofs(GlobalC::exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
//timeval t_start;
	
	const ModuleBase::Element_Basis_Index::Range
		&& range_lcaos = Exx_Abfs::Abfs_Index::construct_range( lcaos );
	const ModuleBase::Element_Basis_Index::IndexLNM
		&& index_lcaos = ModuleBase::Element_Basis_Index::construct_index( range_lcaos );

	const ModuleBase::Element_Basis_Index::Range
		&& range_abfs = Exx_Abfs::Abfs_Index::construct_range( abfs );
	const ModuleBase::Element_Basis_Index::IndexLNM
		&& index_abfs = ModuleBase::Element_Basis_Index::construct_index( range_abfs );
		
//ofs<<range_lcaos<<std::endl;
//ofs<<range_abfs<<std::endl;

	const int Lmax_bak = GlobalC::exx_info.info_ri.abfs_Lmax;
	GlobalC::exx_info.info_ri.abfs_Lmax = std::numeric_limits<int>::min();
	for( size_t T=0; T!=abfs.size(); ++T )
		GlobalC::exx_info.info_ri.abfs_Lmax = std::max( GlobalC::exx_info.info_ri.abfs_Lmax, static_cast<int>(abfs[T].size())-1 );

	Exx_Abfs::Matrix_Orbs21 m_abfslcaos_lcaos;
//gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init( 1, kmesh_times, 1 );
//ofs<<"TIME@m_abfslcaos_lcaos.init\t"<<time_during(t_start)<<std::endl;
//gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init_radial( abfs, lcaos, lcaos );
//ofs<<"TIME@m_abfslcaos_lcaos.init_radial\t"<<time_during(t_start)<<std::endl;

	std::map<size_t,std::map<size_t,set<double>>> delta_R;
	for( size_t it=0; it!=abfs.size(); ++it )
		delta_R[it][it] = {0.0};
//gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init_radial_table(delta_R);
//ofs<<"TIME@m_abfslcaos_lcaos.init_radial_table\t"<<time_during(t_start)<<std::endl;

	GlobalC::exx_info.info_ri.abfs_Lmax = Lmax_bak;
	
	std::vector<std::vector<std::pair<std::vector<double>,ModuleBase::matrix>>> eig(abfs.size());
	for( size_t T=0; T!=abfs.size(); ++T )
	{
		const ModuleBase::matrix && A = m_abfslcaos_lcaos.cal_overlap_matrix(  
			T, 
			T, 
			ModuleBase::Vector3<double>{0,0,0},
			ModuleBase::Vector3<double>{0,0,0},  
			index_abfs, 
			index_lcaos,
			index_lcaos,
			Exx_Abfs::Matrix_Orbs21::Matrix_Order::A2B_A1);
//ofs<<"A:"<<std::endl<<A<<std::endl;
		
		eig[T].resize(abfs[T].size());
		for( size_t L=0; L!=abfs[T].size(); ++L )
		{
//ofs<<"get_sub_matrix:"<<std::endl<<get_sub_matrix( A, T, L, range_abfs, index_abfs )<<std::endl;
//			const matrix A_sub = get_column_mean0_matrix( get_sub_matrix( A, T, L, range_abfs, index_abfs ) );
			const ModuleBase::matrix A_sub = get_sub_matrix( A, T, L, range_abfs, index_abfs );
//ofs<<"A_sub:"<<std::endl<<A_sub<<std::endl;
//ofs<<"transpose:"<<std::endl<<transpose(A_sub)<<std::endl;
//ofs<<"mul:"<<std::endl<<transpose(A_sub) * A_sub<<std::endl;
			ModuleBase::matrix mm = transpose(A_sub) * A_sub;
//ofs<<"mm:"<<std::endl<<mm<<std::endl;
			std::vector<double> eig_value(mm.nr);
			
			int info;
//gettimeofday( &t_start, NULL);
			LapackConnector::dsyev( 'V', 'U', mm, ModuleBase::GlobalFunc::VECTOR_TO_PTR(eig_value), info );
//ofs<<"TIME@LapackConnector::dsyev\t"<<time_during(t_start)<<std::endl;
			if( info )
			{
				std::cout<<std::endl<<"info_dsyev = "<<info<<std::endl;
				mm.print(GlobalV::ofs_warning)<<std::endl;
				std::cout<<"in file "<<__FILE__<<" line "<<__LINE__<<std::endl;
				ModuleBase::QUIT();
			}
			eig[T][L] = make_pair( eig_value, mm );
//for( const auto v : eig_value )
//	ofs<<v<<"\t";
//ofs<<std::endl;
//ofs<<"mm:"<<std::endl<<mm<<std::endl;
		}
	}
//ofs.close();
	
	return eig;
}

ModuleBase::matrix Exx_Abfs::PCA::get_sub_matrix( 
	const ModuleBase::matrix & m,
	const size_t & T,
	const size_t & L,
	const ModuleBase::Element_Basis_Index::Range & range,
	const ModuleBase::Element_Basis_Index::IndexLNM & index )
{
	ModuleBase::TITLE("Exx_Abfs::PCA::get_sub_matrix");
	
	ModuleBase::matrix m_sub( m.nr, range[T][L].N );
	for( size_t ir=0; ir!=m.nr; ++ir )
		for( size_t N=0; N!=range[T][L].N; ++N )
			m_sub( ir, N ) = m( ir, index[T][L][N][0] );
	return m_sub;
}


ModuleBase::matrix Exx_Abfs::PCA::get_column_mean0_matrix( const ModuleBase::matrix & m )
{
	ModuleBase::matrix m_new( m.nr, m.nc );
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
//ofs<<std::endl;
	return m_new;
}
