#include "exx_abfs-pca.h"

#include "exx_abfs-abfs_index.h"
#include "module_ri/Matrix_Orbs11.h"
#include "module_ri/Matrix_Orbs21.h"

#include "../module_base/lapack_connector.h"
#include "../module_base/global_function.h"

#include <cassert>
#include <limits>

#include "../module_ri/test_code/element_basis_index-test.h"		// Peize Lin test
#include <sys/time.h>			// Peize Lin test

static inline void tensor_dsyev(const char jobz, const char uplo, RI::Tensor<double>& a, double* w, int& info)
{
    // reference: dsyev in lapack_connector.h (for ModuleBase::matrix)
    assert(a.shape.size() == 2);
    assert(a.shape[0] == a.shape[1]);
    const int nr = a.shape[0];
    const int nc = a.shape[1];

    double work_tmp;
    constexpr int minus_one = -1;
    dsyev_(&jobz, &uplo, &nr, a.ptr(), &nc, w, &work_tmp, &minus_one, &info);		// get best lwork

    const int lwork = work_tmp;
    std::vector<double> work(std::max(1, lwork));
    dsyev_(&jobz, &uplo, &nr, a.ptr(), &nc, w, ModuleBase::GlobalFunc::VECTOR_TO_PTR(work), &lwork, &info);
};

std::vector<std::vector<std::pair<std::vector<double>, RI::Tensor<double>>>> Exx_Abfs::PCA::cal_PCA(
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

	Matrix_Orbs21 m_abfslcaos_lcaos;
//gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init( 1, kmesh_times, 1 );
//ofs<<"TIME@m_abfslcaos_lcaos.init\t"<<time_during(t_start)<<std::endl;
//gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init_radial( abfs, lcaos, lcaos );
//ofs<<"TIME@m_abfslcaos_lcaos.init_radial\t"<<time_during(t_start)<<std::endl;

	std::map<size_t,std::map<size_t,std::set<double>>> delta_R;
	for( size_t it=0; it!=abfs.size(); ++it )
		delta_R[it][it] = {0.0};
//gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init_radial_table(delta_R);
//ofs<<"TIME@m_abfslcaos_lcaos.init_radial_table\t"<<time_during(t_start)<<std::endl;

	GlobalC::exx_info.info_ri.abfs_Lmax = Lmax_bak;
	
	std::vector<std::vector<std::pair<std::vector<double>,RI::Tensor<double>>>> eig(abfs.size());
	for( size_t T=0; T!=abfs.size(); ++T )
	{
		const RI::Tensor<double> && A = m_abfslcaos_lcaos.cal_overlap_matrix<double>(  
			T, 
			T, 
			ModuleBase::Vector3<double>{0,0,0},
			ModuleBase::Vector3<double>{0,0,0},  
			index_abfs, 
			index_lcaos,
			index_lcaos,
			Matrix_Orbs21::Matrix_Order::A2BA1);
//ofs<<"A:"<<std::endl<<A<<std::endl;
		
		eig[T].resize(abfs[T].size());
		for( size_t L=0; L!=abfs[T].size(); ++L )
		{
//ofs<<"get_sub_matrix:"<<std::endl<<get_sub_matrix( A, T, L, range_abfs, index_abfs )<<std::endl;
//			const matrix A_sub = get_column_mean0_matrix( get_sub_matrix( A, T, L, range_abfs, index_abfs ) );
			const RI::Tensor<double> A_sub = get_sub_matrix( A, T, L, range_abfs, index_abfs );
//ofs<<"A_sub:"<<std::endl<<A_sub<<std::endl;
//ofs<<"transpose:"<<std::endl<<transpose(A_sub)<<std::endl;
//ofs<<"mul:"<<std::endl<<transpose(A_sub) * A_sub<<std::endl;
            RI::Tensor<double> mm = A_sub.transpose() * A_sub;
//ofs<<"mm:"<<std::endl<<mm<<std::endl;
			std::vector<double> eig_value(mm.shape[0]);
			
			int info;
            //gettimeofday( &t_start, NULL);

            tensor_dsyev('V', 'L', mm, ModuleBase::GlobalFunc::VECTOR_TO_PTR(eig_value), info);

            //ofs<<"TIME@LapackConnector::dsyev\t"<<time_during(t_start)<<std::endl;
			if( info )
			{
                std::cout << std::endl << "info_dsyev = " << info << std::endl;
                auto tensor_print = [](RI::Tensor<double>& m, std::ostream& os, const double threshold)
                {
                    for (int ir = 0; ir != m.shape[0]; ++ir)
                    {
                        for (int ic = 0; ic != m.shape[1]; ++ic)
                        {
                            if (std::abs(m(ir, ic)) > threshold)
                                os << m(ir, ic) << "\t";
                            else
                                os << 0 << "\t";
                        }
                        os << std::endl;
                    }
                    os << std::endl;
                };
                tensor_print(mm, GlobalV::ofs_warning, 0.0);
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

RI::Tensor<double> Exx_Abfs::PCA::get_sub_matrix( 
	const RI::Tensor<double> & m,
	const size_t & T,
	const size_t & L,
	const ModuleBase::Element_Basis_Index::Range & range,
	const ModuleBase::Element_Basis_Index::IndexLNM & index )
{
	ModuleBase::TITLE("Exx_Abfs::PCA::get_sub_matrix");
	
    assert(m.shape.size() == 3);
    RI::Tensor<double> m_sub(RI::Shape_Vector{m.shape[0] * m.shape[1], range[T][L].N });
    size_t count = 0;
    for (size_t ir = 0; ir != m.shape[0]; ++ir)
    {
        for (size_t jr = 0; jr != m.shape[1]; ++jr)
        {
            for (size_t N = 0; N != range[T][L].N; ++N)
            {
                m_sub(count, N) = m(ir, jr, index[T][L][N][0]);
            }
            ++count;
        }
    }

	return m_sub;
}


RI::Tensor<double> Exx_Abfs::PCA::get_column_mean0_matrix( const RI::Tensor<double> & m )
{
	RI::Tensor<double> m_new( m.shape);
	for( size_t ic=0; ic!=m.shape[1]; ++ic )
	{
		double sum=0;
		for( size_t ir=0; ir!=m.shape[0]; ++ir )
			sum += m(ir,ic);
		const double mean = sum/m.shape[0];
//ofs<<mean<<"\t";
		for( size_t ir=0; ir!=m.shape[0]; ++ir )
			m_new(ir,ic) = m(ir,ic) - mean;
	}
//ofs<<std::endl;
	return m_new;
}
