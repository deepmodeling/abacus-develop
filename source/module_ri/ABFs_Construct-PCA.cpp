#include "ABFs_Construct-PCA.h"

#include "exx_abfs-abfs_index.h"
#include "../module_base/lapack_connector.h"
#include "../module_base/global_function.h"
#include "../module_base/element_basis_index.h"
#include "../module_base/matrix.h"
#include "../module_ri/Matrix_Orbs11.h"
#include "../module_ri/Matrix_Orbs21.h"

#include <cassert>
#include <limits>

namespace ABFs_Construct::PCA
{
	void tensor_dsyev(const char jobz, const char uplo, RI::Tensor<double> & a, double*const w, int & info)
	{
		// reference: dsyev in lapack_connector.h (for ModuleBase::matrix)
		assert(a.shape.size() == 2);
		assert(a.shape[0] == a.shape[1]);
		const int nr = a.shape[0];
		const int nc = a.shape[1];

		double work_tmp=0.0;
		constexpr int minus_one = -1;
		dsyev_(&jobz, &uplo, &nr, a.ptr(), &nc, w, &work_tmp, &minus_one, &info);		// get best lwork

		const int lwork = work_tmp;
		std::vector<double> work(std::max(1, lwork));
		dsyev_(&jobz, &uplo, &nr, a.ptr(), &nc, w, work.data(), &lwork, &info);
	}

	RI::Tensor<double> get_sub_matrix( 
		const RI::Tensor<double> & m,		// size: (lcaos, lcaos, abfs)
		const std::size_t & T,
		const std::size_t & L,
		const ModuleBase::Element_Basis_Index::Range & range,
		const ModuleBase::Element_Basis_Index::IndexLNM & index )
	{
		ModuleBase::TITLE("ABFs_Construct::PCA::get_sub_matrix");		
		assert(m.shape.size() == 3);
		RI::Tensor<double> m_sub({ m.shape[0], m.shape[1], range[T][L].N });
		for (std::size_t ir=0; ir!=m.shape[0]; ++ir)
			for (std::size_t jr=0; jr!=m.shape[1]; ++jr)
				for (std::size_t N=0; N!=range[T][L].N; ++N)
					m_sub(ir, jr, N) = m(ir, jr, index[T][L][N][0]);
		m_sub = m_sub.reshape({ m.shape[0] * m.shape[1], range[T][L].N });
		return m_sub;
	}

	RI::Tensor<double> get_column_mean0_matrix( const RI::Tensor<double> & m )
	{
		ModuleBase::TITLE("ABFs_Construct::PCA::get_column_mean0_matrix");
		RI::Tensor<double> m_new( m.shape);
		for( std::size_t ic=0; ic!=m.shape[1]; ++ic )
		{
			double sum=0;
			for( std::size_t ir=0; ir!=m.shape[0]; ++ir )
				sum += m(ir,ic);
			const double mean = sum/m.shape[0];
			for( std::size_t ir=0; ir!=m.shape[0]; ++ir )
				m_new(ir,ic) = m(ir,ic) - mean;
		}
		return m_new;
	}

	std::vector<std::vector<std::pair<std::vector<double>, RI::Tensor<double>>>> cal_PCA(
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &lcaos, 
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs,
		const double kmesh_times )
	{
		ModuleBase::TITLE("ABFs_Construct::PCA::cal_PCA");
		
		const ModuleBase::Element_Basis_Index::Range
			range_lcaos = Exx_Abfs::Abfs_Index::construct_range( lcaos );
		const ModuleBase::Element_Basis_Index::IndexLNM
			index_lcaos = ModuleBase::Element_Basis_Index::construct_index( range_lcaos );

		const ModuleBase::Element_Basis_Index::Range
			range_abfs = Exx_Abfs::Abfs_Index::construct_range( abfs );
		const ModuleBase::Element_Basis_Index::IndexLNM
			index_abfs = ModuleBase::Element_Basis_Index::construct_index( range_abfs );

		const int Lmax_bak = GlobalC::exx_info.info_ri.abfs_Lmax;
		GlobalC::exx_info.info_ri.abfs_Lmax = std::numeric_limits<int>::min();
		for( std::size_t T=0; T!=abfs.size(); ++T )
			GlobalC::exx_info.info_ri.abfs_Lmax = std::max( GlobalC::exx_info.info_ri.abfs_Lmax, static_cast<int>(abfs[T].size())-1 );

		Matrix_Orbs21 m_abfslcaos_lcaos;
		m_abfslcaos_lcaos.init( 1, kmesh_times, 1 );
		m_abfslcaos_lcaos.init_radial( abfs, lcaos, lcaos );

		std::map<std::size_t,std::map<std::size_t,std::set<double>>> delta_R;
		for( std::size_t it=0; it!=abfs.size(); ++it )
			delta_R[it][it] = {0.0};
		m_abfslcaos_lcaos.init_radial_table(delta_R);

		GlobalC::exx_info.info_ri.abfs_Lmax = Lmax_bak;
		
		std::vector<std::vector<std::pair<std::vector<double>,RI::Tensor<double>>>> eig(abfs.size());
		for( std::size_t T=0; T!=abfs.size(); ++T )
		{
			const RI::Tensor<double> A = m_abfslcaos_lcaos.cal_overlap_matrix<double>(  
				T, 
				T, 
				ModuleBase::Vector3<double>{0,0,0},
				ModuleBase::Vector3<double>{0,0,0},  
				index_abfs, 
				index_lcaos,
				index_lcaos,
				Matrix_Orbs21::Matrix_Order::A2BA1);
			
			eig[T].resize(abfs[T].size());
			for( std::size_t L=0; L!=abfs[T].size(); ++L )
			{
				const RI::Tensor<double> A_sub = get_sub_matrix( A, T, L, range_abfs, index_abfs );
				RI::Tensor<double> mm = A_sub.transpose() * A_sub;
				std::vector<double> eig_value(mm.shape[0]);
				
				int info=1;

				tensor_dsyev('V', 'L', mm, eig_value.data(), info);

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
				eig[T][L] = std::make_pair( eig_value, mm );
			}
		}
		
		return eig;
	}

}	// namespace ABFs_Construct::PCA
