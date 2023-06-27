#include "exx_opt_orb.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "exx_abfs.h"
#include "exx_abfs-abfs_index.h"
#include "exx_abfs-construct_orbs.h"
#include "exx_abfs-jle.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_ri/Matrix_Orbs11.h"
#include "module_ri/Matrix_Orbs21.h"
#include "module_ri/Matrix_Orbs22.h"
#include "module_ri/LRI_CV_Tools.h"

#include "../module_ri/test_code/element_basis_index-test.h"
#include "../module_ri/test_code/test_function.h"

void Exx_Opt_Orb::generate_matrix(const K_Vectors &kv) const
{
// std::ofstream ofs_mpi(GlobalC::exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);

	ModuleBase::TITLE("Exx_Opt_Orb::generate_matrix");
// ofs_mpi<<"memory:\t"<<get_memory(10)<<std::endl;

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		lcaos = Exx_Abfs::Construct_Orbs::change_orbs( GlobalC::ORB, this->kmesh_times );

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		abfs = Exx_Abfs::Construct_Orbs::abfs_same_atom( lcaos, this->kmesh_times, GlobalC::exx_info.info_ri.pca_threshold );

// ofs_mpi<<"memory:\t"<<get_memory(10)<<std::endl;
	
	Exx_Abfs::Jle jle;
	jle.init_jle( this->kmesh_times );

// ofs_mpi<<"memory:\t"<<get_memory(10)<<std::endl;
	
	GlobalC::exx_info.info_ri.abfs_Lmax = Exx_Abfs::Jle::Lmax;
	for( size_t T=0; T!=abfs.size(); ++T )
		GlobalC::exx_info.info_ri.abfs_Lmax = std::max( GlobalC::exx_info.info_ri.abfs_Lmax, static_cast<int>(abfs[T].size())-1 );

	const ModuleBase::Element_Basis_Index::Range    range_lcaos = Exx_Abfs::Abfs_Index::construct_range( lcaos );
	const ModuleBase::Element_Basis_Index::IndexLNM index_lcaos = ModuleBase::Element_Basis_Index::construct_index( range_lcaos );

	const ModuleBase::Element_Basis_Index::Range    range_abfs = Exx_Abfs::Abfs_Index::construct_range( abfs );
	const ModuleBase::Element_Basis_Index::IndexLNM index_abfs = ModuleBase::Element_Basis_Index::construct_index( range_abfs );

	const ModuleBase::Element_Basis_Index::Range    range_jys = Exx_Abfs::Abfs_Index::construct_range( jle.jle );
	const ModuleBase::Element_Basis_Index::IndexLNM index_jys = ModuleBase::Element_Basis_Index::construct_index( range_jys );

// ofs_mpi<<range_lcaos<<std::endl;
// ofs_mpi<<range_abfs<<std::endl;
// ofs_mpi<<range_jys<<std::endl;

	std::map<size_t,std::map<size_t,set<double>>> radial_R = get_radial_R();
#if TEST_EXX_RADIAL==2
	{
		for(const auto & rA : radial_R)
			for(const auto & rB : rA.second)
			{
				ofs_mpi<<rA.first<<"\t"<<rB.first<<":\t";
				for(const auto & rC : rB.second)
					ofs_mpi<<rC<<"\t";
				ofs_mpi<<std::endl;
			}
	}
#endif

// ofs_mpi<<"memory:\t"<<get_memory(10)<<std::endl;

	// < lcaos lcaos | lcaos lcaos >
	const auto ms_lcaoslcaos_lcaoslcaos = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>> 
	{
		Matrix_Orbs22 m_lcaoslcaos_lcaoslcaos;
		m_lcaoslcaos_lcaoslcaos.init( 1, this->kmesh_times, 1 );
		m_lcaoslcaos_lcaoslcaos.init_radial( lcaos, lcaos, lcaos, lcaos );
		#if TEST_EXX_RADIAL>=1
		m_lcaoslcaos_lcaoslcaos.init_radial_table(radial_R);
		#else
		m_lcaoslcaos_lcaoslcaos.init_radial_table();
		#endif
		return m_lcaoslcaos_lcaoslcaos.cal_overlap_matrix_all<double>( index_lcaos, index_lcaos, index_lcaos, index_lcaos);
	}();
	
// ofs_mpi<<"memory:\t"<<get_memory(10)<<std::endl;

	// < lcaos lcaos | jys >
	const auto ms_lcaoslcaos_jys = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<RI::Tensor<double>>>>>>
	{
		Matrix_Orbs21 m_jyslcaos_lcaos;
		m_jyslcaos_lcaos.init( 1, this->kmesh_times, 1 );
		m_jyslcaos_lcaos.init_radial( jle.jle, lcaos, lcaos );
		#if TEST_EXX_RADIAL>=1
		m_jyslcaos_lcaos.init_radial_table(radial_R);
		#else
		m_jyslcaos_lcaos.init_radial_table();
		#endif
		return m_jyslcaos_lcaos.cal_overlap_matrix_all<double>( index_jys, index_lcaos, index_lcaos );
	}();

// ofs_mpi<<"memory:\t"<<get_memory(10)<<std::endl;

	// < jys | jys >
	const auto ms_jys_jys = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>>
	{
		Matrix_Orbs11 m_jys_jys;
		m_jys_jys.init( 2, this->kmesh_times, 1 );
		m_jys_jys.init_radial( jle.jle, jle.jle );
		#if TEST_EXX_RADIAL>=1
		m_jys_jys.init_radial_table(radial_R);
		#else
		m_jys_jys.init_radial_table();
		#endif
		return m_jys_jys.cal_overlap_matrix_all<double>( index_jys, index_jys );
	}();

// ofs_mpi<<"memory:\t"<<get_memory(10)<<std::endl;

	// < abfs | abfs >
	const auto ms_abfs_abfs = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>>
	{
		Matrix_Orbs11 m_abfs_abfs;
		m_abfs_abfs.init( 2, this->kmesh_times, 1 );
		m_abfs_abfs.init_radial( abfs, abfs );
		#if TEST_EXX_RADIAL>=1
		m_abfs_abfs.init_radial_table(radial_R);
		#else
		m_abfs_abfs.init_radial_table();
		#endif
		return m_abfs_abfs.cal_overlap_matrix_all<double>( index_abfs, index_abfs );
	}();

// ofs_mpi<<"memory:\t"<<get_memory(10)<<std::endl;

	// < lcaos lcaos | abfs >
	const auto ms_lcaoslcaos_abfs = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<RI::Tensor<double>>>>>>
	{
		Matrix_Orbs21 m_abfslcaos_lcaos;
		m_abfslcaos_lcaos.init( 1, this->kmesh_times, 1 );
		m_abfslcaos_lcaos.init_radial( abfs, lcaos, lcaos );
		#if TEST_EXX_RADIAL>=1
		m_abfslcaos_lcaos.init_radial_table(radial_R);
		#else
		m_abfslcaos_lcaos.init_radial_table();
		#endif
		return m_abfslcaos_lcaos.cal_overlap_matrix_all<double>( index_abfs, index_lcaos, index_lcaos );
	}();

// ofs_mpi<<"memory:\t"<<get_memory(10)<<std::endl;

	// < jys | abfs >
	const auto ms_jys_abfs = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>>
	{
		Matrix_Orbs11 m_jys_abfs;
		m_jys_abfs.init( 2, this->kmesh_times, 1 );
		m_jys_abfs.init_radial( jle.jle, abfs );
		#if TEST_EXX_RADIAL>=1
		m_jys_abfs.init_radial_table(radial_R);
		#else
		m_jys_abfs.init_radial_table();
		#endif
		return m_jys_abfs.cal_overlap_matrix_all<double>( index_jys, index_abfs );
	}();

// ofs_mpi<<"memory:\t"<<get_memory(10)<<std::endl;

// ofs_matrixes(GlobalC::exx_lcao.test_dir.matrix+"ms_jys_jys",ms_jys_jys);
// ofs_matrixes(GlobalC::exx_lcao.test_dir.matrix+"ms_lcaoslcaos_jys",ms_lcaoslcaos_jys);
// ofs_matrixes(GlobalC::exx_lcao.test_dir.matrix+"ms_lcaoslcaos_lcaoslcaos",ms_lcaoslcaos_lcaoslcaos);
// ofs_matrixes(GlobalC::exx_lcao.test_dir.matrix+"ms_abfs_abfs",ms_abfs_abfs);
// ofs_matrixes(GlobalC::exx_lcao.test_dir.matrix+"ms_lcaoslcaos_abfs",ms_lcaoslcaos_abfs);
// ofs_matrixes(GlobalC::exx_lcao.test_dir.matrix+"ms_jys_abfs",ms_jys_abfs);

	for( size_t TA=0; TA!=GlobalC::ucell.ntype; ++TA )
	{
		for( size_t IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA )
		{
			for( size_t TB=0; TB!=GlobalC::ucell.ntype; ++TB )
			{
				for( size_t IB=0; IB!=GlobalC::ucell.atoms[TB].na; ++IB )
				{
					if( TA==TB && IA==IB )
					{
						const size_t T=TA, I=IA;
						if(GlobalC::exx_info.info_ri.pca_threshold<=1)
						{
							// < abfs | abfs >.I
							const std::vector<std::vector<RI::Tensor<double>>> ms_abfs_abfs_I = cal_I( ms_abfs_abfs, T,I,T,I );
							// < lcaos lcaos | lcaos lcaos > - < lcaos lcaos | abfs > * < abfs | abfs >.I * < abfs | lcaos lcaos >
							const RI::Tensor<double> m_lcaoslcaos_lcaoslcaos_proj =
								cal_proj(
									ms_lcaoslcaos_lcaoslcaos.at(T).at(I).at(T).at(I),
									ms_lcaoslcaos_abfs.at(T).at(I).at(T).at(I),
									ms_abfs_abfs_I,
									ms_lcaoslcaos_abfs.at(T).at(I).at(T).at(I));
							// < lcaos lcaos | jys > - < lcaos lcaos | abfs > * < abfs | abfs >.I * < abfs | jys >
							const std::vector<RI::Tensor<double>> m_lcaoslcaos_jys_proj =
								{cal_proj(
									ms_lcaoslcaos_jys.at(T).at(I).at(T).at(I)[0],
									ms_lcaoslcaos_abfs.at(T).at(I).at(T).at(I),
									ms_abfs_abfs_I,
									{ms_jys_abfs.at(T).at(I).at(T).at(I)})};
							// < jys | jys > - < jys | abfs > * < abfs | abfs >.I * < abfs | jys >
							const std::vector<std::vector<RI::Tensor<double>>> m_jys_jys_proj =
								{{cal_proj(
									ms_jys_jys.at(T).at(I).at(T).at(I),
									{ms_jys_abfs.at(T).at(I).at(T).at(I)},
									ms_abfs_abfs_I,
									{ms_jys_abfs.at(T).at(I).at(T).at(I)})}};
							print_matrix(kv,
								"matrix",
								m_lcaoslcaos_jys_proj,
								m_jys_jys_proj,
								m_lcaoslcaos_lcaoslcaos_proj,
								TA, IA, TB, IB,
								range_jys, index_jys );
						}
						else
						{
							print_matrix(kv,
								"matrix",
								ms_lcaoslcaos_jys.at(T).at(I).at(T).at(I),
								{{ms_jys_jys.at(T).at(I).at(T).at(I)}},
								ms_lcaoslcaos_lcaoslcaos.at(T).at(I).at(T).at(I),
								TA, IA, TB, IB,
								range_jys, index_jys );
						}
					}
					else
					{
						if(GlobalC::exx_info.info_ri.pca_threshold<=1)
						{
							// < abfs | abfs >.I
							const std::vector<std::vector<RI::Tensor<double>>> ms_abfs_abfs_I = cal_I( ms_abfs_abfs, TA,IA,TB,IB );
							// < lcaos lcaos | lcaos lcaos > - < lcaos lcaos | abfs > * < abfs | abfs >.I * < abfs | lcaos lcaos >
							const RI::Tensor<double> m_lcaoslcaos_lcaoslcaos_proj =
								cal_proj(
									ms_lcaoslcaos_lcaoslcaos.at(TA).at(IA).at(TB).at(IB),
									ms_lcaoslcaos_abfs.at(TA).at(IA).at(TB).at(IB),
									ms_abfs_abfs_I,
									ms_lcaoslcaos_abfs.at(TA).at(IA).at(TB).at(IB));
							// < lcaos lcaos | jys > - < lcaos lcaos | abfs > * < abfs | abfs >.I * < abfs | jys >
							const std::vector<RI::Tensor<double>> m_lcaoslcaos_jys_proj =
								{cal_proj(
									ms_lcaoslcaos_jys.at(TA).at(IA).at(TB).at(IB)[0],
									ms_lcaoslcaos_abfs.at(TA).at(IA).at(TB).at(IB),
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TA).at(IA).at(TA).at(IA), ms_jys_abfs.at(TA).at(IA).at(TB).at(IB) }),
								 cal_proj(
									ms_lcaoslcaos_jys.at(TA).at(IA).at(TB).at(IB)[1],
									ms_lcaoslcaos_abfs.at(TA).at(IA).at(TB).at(IB),
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TB).at(IB).at(TA).at(IA), ms_jys_abfs.at(TB).at(IB).at(TB).at(IB) })};
							// < jys | jys > - < jys | abfs > * < abfs | abfs >.I * < abfs | jys >
							const std::vector<std::vector<RI::Tensor<double>>> m_jys_jys_proj =
								{{cal_proj(
									ms_jys_jys.at(TA).at(IA).at(TA).at(IA),
									{ ms_jys_abfs.at(TA).at(IA).at(TA).at(IA), ms_jys_abfs.at(TA).at(IA).at(TB).at(IB) },
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TA).at(IA).at(TA).at(IA), ms_jys_abfs.at(TA).at(IA).at(TB).at(IB) }),
								  cal_proj(
									ms_jys_jys.at(TA).at(IA).at(TB).at(IB),
									{ ms_jys_abfs.at(TA).at(IA).at(TA).at(IA), ms_jys_abfs.at(TA).at(IA).at(TB).at(IB) },
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TB).at(IB).at(TA).at(IA), ms_jys_abfs.at(TB).at(IB).at(TB).at(IB) }) },
								 {cal_proj(
									ms_jys_jys.at(TB).at(IB).at(TA).at(IA),
									{ ms_jys_abfs.at(TB).at(IB).at(TA).at(IA), ms_jys_abfs.at(TB).at(IB).at(TB).at(IB) },
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TA).at(IA).at(TA).at(IA), ms_jys_abfs.at(TA).at(IA).at(TB).at(IB) }),
								  cal_proj(
									ms_jys_jys.at(TB).at(IB).at(TB).at(IB),
									{ ms_jys_abfs.at(TB).at(IB).at(TA).at(IA), ms_jys_abfs.at(TB).at(IB).at(TB).at(IB) },
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TB).at(IB).at(TA).at(IA), ms_jys_abfs.at(TB).at(IB).at(TB).at(IB) }) }};
							print_matrix(kv,
								"matrix",
								m_lcaoslcaos_jys_proj,
								m_jys_jys_proj,
								m_lcaoslcaos_lcaoslcaos_proj,
								TA, IA, TB, IB,
								range_jys, index_jys );
						}
						else
						{
							print_matrix(kv,
								"matrix",
								ms_lcaoslcaos_jys.at(TA).at(IA).at(TB).at(IB),
								{{ms_jys_jys.at(TA).at(IA).at(TA).at(IA), ms_jys_jys.at(TA).at(IA).at(TB).at(IB)},
								 {ms_jys_jys.at(TB).at(IB).at(TA).at(IA), ms_jys_jys.at(TB).at(IB).at(TB).at(IB)}},
								ms_lcaoslcaos_lcaoslcaos.at(TA).at(IA).at(TB).at(IB),
								TA, IA, TB, IB,
								range_jys, index_jys );
						}
					}
				}
			}
		}
	}
}

// m_big - m_left * m_middle * m_right.T
RI::Tensor<double> Exx_Opt_Orb::cal_proj( 
	const RI::Tensor<double> & m_big, 
	const std::vector<RI::Tensor<double>> & m_left, 
	const std::vector<std::vector<RI::Tensor<double>>> & m_middle, 
	const std::vector<RI::Tensor<double>> & m_right ) const
{
    ModuleBase::TITLE("Exx_Opt_Orb::cal_proj");

//auto print_nrc = [](const matrix & m){ std::cout<<"\t"<<m.nr<<"\t"<<m.nc<<std::endl; };

	RI::Tensor<double> m_proj = m_big;
//print_nrc(m_proj);
	for( size_t il=0; il!=m_left.size(); ++il )
	{
		for( size_t ir=0; ir!=m_right.size(); ++ir )
		{
//std::cout<<il<<"\t"<<ir<<std::endl;
//print_nrc(m_left[il]);
//print_nrc(m_middle[il][ir]);
//print_nrc(m_right[ir]);
			m_proj = m_proj - m_left[il] *m_middle[il][ir] * m_right[ir].transpose();
		}
	}
	return m_proj;
}

std::vector<std::vector<RI::Tensor<double>>> Exx_Opt_Orb::cal_I(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>> &ms,
	const size_t TA, const size_t IA, const size_t TB, const size_t IB ) const
{
	ModuleBase::TITLE("Exx_Opt_Orb::cal_I");

	if( TA==TB && IA==IB )
	{
		std::vector<std::vector<RI::Tensor<double>>> m_I
        {{ RI::Tensor<double>(ms.at(TA).at(IA).at(TA).at(IA).shape) }};
		return LRI_CV_Tools::cal_I(m_I);
	}
	else
	{
		std::vector<std::vector<RI::Tensor<double>>> m_I
			{{ RI::Tensor<double>(ms.at(TA).at(IA).at(TA).at(IA).shape),
			   RI::Tensor<double>(ms.at(TA).at(IA).at(TB).at(IB).shape) },
			 { RI::Tensor<double>(ms.at(TB).at(IB).at(TA).at(IA).shape),
			   RI::Tensor<double>(ms.at(TB).at(IB).at(TB).at(IB).shape) }};
		return LRI_CV_Tools::cal_I(m_I);
	}
}

std::map<size_t,std::map<size_t,set<double>>> Exx_Opt_Orb::get_radial_R() const
{
	ModuleBase::TITLE("Exx_Opt_Orb::get_radial_R");
	std::map<size_t,std::map<size_t,set<double>>> radial_R;
	for( size_t TA=0; TA!=GlobalC::ucell.ntype; ++TA )
		for( size_t IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA )
			for( size_t TB=0; TB!=GlobalC::ucell.ntype; ++TB )
				for( size_t IB=0; IB!=GlobalC::ucell.atoms[TB].na; ++IB )
				{
					const ModuleBase::Vector3<double> &tauA = GlobalC::ucell.atoms[TA].tau[IA];
					const ModuleBase::Vector3<double> &tauB = GlobalC::ucell.atoms[TB].tau[IB];
					const double delta_R = (-tauA+tauB).norm();
					radial_R[TA][TB].insert( delta_R );
					radial_R[TB][TA].insert( delta_R );
				}
	return radial_R;
}