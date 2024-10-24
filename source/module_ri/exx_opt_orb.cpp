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
#include <RI/global/Tensor_Multiply.h>

void Exx_Opt_Orb::generate_matrix(
	const Exx_Info::Exx_Info_Opt_ABFs &info,
	const K_Vectors &kv,
	const LCAO_Orbitals &orb)
{
	ModuleBase::TITLE("Exx_Opt_Orb::generate_matrix");

	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		lcaos = Exx_Abfs::Construct_Orbs::change_orbs( orb, info.kmesh_times );
	Exx_Abfs::Construct_Orbs::filter_empty_orbs(lcaos);

	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		abfs = Exx_Abfs::Construct_Orbs::abfs_same_atom( orb, lcaos, info.kmesh_times, GlobalC::exx_info.info_ri.pca_threshold );
	Exx_Abfs::Construct_Orbs::filter_empty_orbs(abfs);

	std::vector< std::vector< std::vector< Numerical_Orbital_Lm>>>
		jle = Exx_Abfs::Jle::init_jle( info, info.kmesh_times, orb );
	Exx_Abfs::Construct_Orbs::filter_empty_orbs(jle);

	GlobalC::exx_info.info_ri.abfs_Lmax = info.abfs_Lmax;
	for( size_t T=0; T!=abfs.size(); ++T )
		{ GlobalC::exx_info.info_ri.abfs_Lmax = std::max( GlobalC::exx_info.info_ri.abfs_Lmax, static_cast<int>(abfs[T].size())-1 ); }

	const ModuleBase::Element_Basis_Index::Range    range_lcaos = Exx_Abfs::Abfs_Index::construct_range( lcaos );
	const ModuleBase::Element_Basis_Index::IndexLNM index_lcaos = ModuleBase::Element_Basis_Index::construct_index( range_lcaos );

	const ModuleBase::Element_Basis_Index::Range    range_abfs = Exx_Abfs::Abfs_Index::construct_range( abfs );
	const ModuleBase::Element_Basis_Index::IndexLNM index_abfs = ModuleBase::Element_Basis_Index::construct_index( range_abfs );

	const ModuleBase::Element_Basis_Index::Range    range_jys = Exx_Abfs::Abfs_Index::construct_range( jle );
	const ModuleBase::Element_Basis_Index::IndexLNM index_jys = ModuleBase::Element_Basis_Index::construct_index( range_jys );

	Exx_Abfs::Construct_Orbs::print_orbs_size(abfs, GlobalV::ofs_running);

	std::map<size_t,std::map<size_t,std::set<double>>> radial_R = get_radial_R();

	// < lcaos lcaos | lcaos lcaos >
	const auto ms_lcaoslcaos_lcaoslcaos = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>>
	{
		Matrix_Orbs22 m_lcaoslcaos_lcaoslcaos;
		m_lcaoslcaos_lcaoslcaos.init( 1, orb, info.kmesh_times, 1 );
		m_lcaoslcaos_lcaoslcaos.init_radial( lcaos, lcaos, lcaos, lcaos );
		#if TEST_EXX_RADIAL>=1
		m_lcaoslcaos_lcaoslcaos.init_radial_table(radial_R);
		#else
		m_lcaoslcaos_lcaoslcaos.init_radial_table();
		#endif
		return m_lcaoslcaos_lcaoslcaos.cal_overlap_matrix_all<double>( index_lcaos, index_lcaos, index_lcaos, index_lcaos);
	}();


	// < lcaos lcaos | jys >
	const auto ms_lcaoslcaos_jys = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<RI::Tensor<double>>>>>>
	{
		Matrix_Orbs21 m_jyslcaos_lcaos;
		m_jyslcaos_lcaos.init( 1, orb, info.kmesh_times, 1 );
		m_jyslcaos_lcaos.init_radial( jle, lcaos, lcaos );
		#if TEST_EXX_RADIAL>=1
		m_jyslcaos_lcaos.init_radial_table(radial_R);
		#else
		m_jyslcaos_lcaos.init_radial_table();
		#endif
		return m_jyslcaos_lcaos.cal_overlap_matrix_all<double>( index_jys, index_lcaos, index_lcaos );
	}();

	// < jys | jys >
	const auto ms_jys_jys = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>>
	{
		Matrix_Orbs11 m_jys_jys;
		m_jys_jys.init( 2, orb, info.kmesh_times, 1 );
		m_jys_jys.init_radial( jle, jle );
		#if TEST_EXX_RADIAL>=1
		m_jys_jys.init_radial_table(radial_R);
		#else
		m_jys_jys.init_radial_table();
		#endif
		return m_jys_jys.cal_overlap_matrix_all<double>( index_jys, index_jys );
	}();

	// < abfs | abfs >
	const auto ms_abfs_abfs = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>>
	{
		Matrix_Orbs11 m_abfs_abfs;
		m_abfs_abfs.init( 2, orb, info.kmesh_times, 1 );
		m_abfs_abfs.init_radial( abfs, abfs );
		#if TEST_EXX_RADIAL>=1
		m_abfs_abfs.init_radial_table(radial_R);
		#else
		m_abfs_abfs.init_radial_table();
		#endif
		return m_abfs_abfs.cal_overlap_matrix_all<double>( index_abfs, index_abfs );
	}();

	// < lcaos lcaos | abfs >
	const auto ms_lcaoslcaos_abfs = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<RI::Tensor<double>>>>>>
	{
		Matrix_Orbs21 m_abfslcaos_lcaos;
		m_abfslcaos_lcaos.init( 1, orb, info.kmesh_times, 1 );
		m_abfslcaos_lcaos.init_radial( abfs, lcaos, lcaos );
		#if TEST_EXX_RADIAL>=1
		m_abfslcaos_lcaos.init_radial_table(radial_R);
		#else
		m_abfslcaos_lcaos.init_radial_table();
		#endif
		return m_abfslcaos_lcaos.cal_overlap_matrix_all<double>( index_abfs, index_lcaos, index_lcaos );
	}();

	// < jys | abfs >
	const auto ms_jys_abfs = [&]() -> std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>>
	{
		Matrix_Orbs11 m_jys_abfs;
		m_jys_abfs.init( 2, orb, info.kmesh_times, 1 );
		m_jys_abfs.init_radial( jle, abfs );
		#if TEST_EXX_RADIAL>=1
		m_jys_abfs.init_radial_table(radial_R);
		#else
		m_jys_abfs.init_radial_table();
		#endif
		return m_jys_abfs.cal_overlap_matrix_all<double>( index_jys, index_abfs );
	}();

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
								cal_proj_22(
									ms_lcaoslcaos_lcaoslcaos.at(T).at(I).at(T).at(I),
									ms_lcaoslcaos_abfs.at(T).at(I).at(T).at(I),
									ms_abfs_abfs_I,
									ms_lcaoslcaos_abfs.at(T).at(I).at(T).at(I));
							// < lcaos lcaos | jys > - < lcaos lcaos | abfs > * < abfs | abfs >.I * < abfs | jys >
							const std::vector<RI::Tensor<double>> m_lcaoslcaos_jys_proj =
								{cal_proj_21(
									ms_lcaoslcaos_jys.at(T).at(I).at(T).at(I)[0],
									ms_lcaoslcaos_abfs.at(T).at(I).at(T).at(I),
									ms_abfs_abfs_I,
									{ms_jys_abfs.at(T).at(I).at(T).at(I)})};
							// < jys | jys > - < jys | abfs > * < abfs | abfs >.I * < abfs | jys >
							const std::vector<std::vector<RI::Tensor<double>>> m_jys_jys_proj =
								{{cal_proj_11(
									ms_jys_jys.at(T).at(I).at(T).at(I),
									{ms_jys_abfs.at(T).at(I).at(T).at(I)},
									ms_abfs_abfs_I,
									{ms_jys_abfs.at(T).at(I).at(T).at(I)})}};
							print_matrix(
								GlobalC::exx_info.info_opt_abfs,
								kv,
								PARAM.globalv.global_out_dir+"/matrix-opt-abfs",
								m_lcaoslcaos_jys_proj,
								m_jys_jys_proj,
								m_lcaoslcaos_lcaoslcaos_proj,
								TA, IA, TB, IB,
								orb.cutoffs(),
								range_jys, index_jys );
						}
						else
						{
							print_matrix(
								GlobalC::exx_info.info_opt_abfs,
								kv,
								PARAM.globalv.global_out_dir+"/matrix-opt-abfs",
								ms_lcaoslcaos_jys.at(T).at(I).at(T).at(I),
								{{ms_jys_jys.at(T).at(I).at(T).at(I)}},
								ms_lcaoslcaos_lcaoslcaos.at(T).at(I).at(T).at(I),
								TA, IA, TB, IB,
								orb.cutoffs(),
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
								cal_proj_22(
									ms_lcaoslcaos_lcaoslcaos.at(TA).at(IA).at(TB).at(IB),
									ms_lcaoslcaos_abfs.at(TA).at(IA).at(TB).at(IB),
									ms_abfs_abfs_I,
									ms_lcaoslcaos_abfs.at(TA).at(IA).at(TB).at(IB));
							// < lcaos lcaos | jys > - < lcaos lcaos | abfs > * < abfs | abfs >.I * < abfs | jys >
							const std::vector<RI::Tensor<double>> m_lcaoslcaos_jys_proj =
								{cal_proj_21(
									ms_lcaoslcaos_jys.at(TA).at(IA).at(TB).at(IB)[0],
									ms_lcaoslcaos_abfs.at(TA).at(IA).at(TB).at(IB),
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TA).at(IA).at(TA).at(IA), ms_jys_abfs.at(TA).at(IA).at(TB).at(IB) }),
								 cal_proj_21(
									ms_lcaoslcaos_jys.at(TA).at(IA).at(TB).at(IB)[1],
									ms_lcaoslcaos_abfs.at(TA).at(IA).at(TB).at(IB),
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TB).at(IB).at(TA).at(IA), ms_jys_abfs.at(TB).at(IB).at(TB).at(IB) })};
							// < jys | jys > - < jys | abfs > * < abfs | abfs >.I * < abfs | jys >
							const std::vector<std::vector<RI::Tensor<double>>> m_jys_jys_proj =
								{{cal_proj_11(
									ms_jys_jys.at(TA).at(IA).at(TA).at(IA),
									{ ms_jys_abfs.at(TA).at(IA).at(TA).at(IA), ms_jys_abfs.at(TA).at(IA).at(TB).at(IB) },
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TA).at(IA).at(TA).at(IA), ms_jys_abfs.at(TA).at(IA).at(TB).at(IB) }),
								  cal_proj_11(
									ms_jys_jys.at(TA).at(IA).at(TB).at(IB),
									{ ms_jys_abfs.at(TA).at(IA).at(TA).at(IA), ms_jys_abfs.at(TA).at(IA).at(TB).at(IB) },
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TB).at(IB).at(TA).at(IA), ms_jys_abfs.at(TB).at(IB).at(TB).at(IB) }) },
								 {cal_proj_11(
									ms_jys_jys.at(TB).at(IB).at(TA).at(IA),
									{ ms_jys_abfs.at(TB).at(IB).at(TA).at(IA), ms_jys_abfs.at(TB).at(IB).at(TB).at(IB) },
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TA).at(IA).at(TA).at(IA), ms_jys_abfs.at(TA).at(IA).at(TB).at(IB) }),
								  cal_proj_11(
									ms_jys_jys.at(TB).at(IB).at(TB).at(IB),
									{ ms_jys_abfs.at(TB).at(IB).at(TA).at(IA), ms_jys_abfs.at(TB).at(IB).at(TB).at(IB) },
									ms_abfs_abfs_I,
									{ ms_jys_abfs.at(TB).at(IB).at(TA).at(IA), ms_jys_abfs.at(TB).at(IB).at(TB).at(IB) }) }};
							print_matrix(
								GlobalC::exx_info.info_opt_abfs,
								kv,
								PARAM.globalv.global_out_dir+"/matrix-opt-abfs",
								m_lcaoslcaos_jys_proj,
								m_jys_jys_proj,
								m_lcaoslcaos_lcaoslcaos_proj,
								TA, IA, TB, IB,
								orb.cutoffs(),
								range_jys, index_jys );
						}
						else
						{
							print_matrix(
								GlobalC::exx_info.info_opt_abfs,
								kv,
								PARAM.globalv.global_out_dir+"/matrix-opt-abfs",
								ms_lcaoslcaos_jys.at(TA).at(IA).at(TB).at(IB),
								{{ms_jys_jys.at(TA).at(IA).at(TA).at(IA), ms_jys_jys.at(TA).at(IA).at(TB).at(IB)},
								 {ms_jys_jys.at(TB).at(IB).at(TA).at(IA), ms_jys_jys.at(TB).at(IB).at(TB).at(IB)}},
								ms_lcaoslcaos_lcaoslcaos.at(TA).at(IA).at(TB).at(IB),
								TA, IA, TB, IB,
								orb.cutoffs(),
								range_jys, index_jys );
						}
					}
				}
			}
		}
	}
}

// m_big - m_left * m_middle * m_right.T
RI::Tensor<double> Exx_Opt_Orb::cal_proj_22(
	const RI::Tensor<double> & m_big,
	const std::vector<RI::Tensor<double>> & m_left,
	const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
	const std::vector<RI::Tensor<double>> & m_right )
{
	ModuleBase::TITLE("Exx_Opt_Orb::cal_proj");
	RI::Tensor<double> m_proj = m_big.copy();
	for( size_t il=0; il!=m_left.size(); ++il )
	{
		for( size_t ir=0; ir!=m_right.size(); ++ir )
		{
			// m_proj = m_proj - m_left[il] * m_middle[il][ir] * m_right[ir].T;
			const RI::Tensor<double> m_lm = RI::Tensor_Multiply::x0x1y1_x0x1a_ay1(m_left[il], m_middle[il][ir]);
			const RI::Tensor<double> m_lmr = RI::Tensor_Multiply::x0x1y0y1_x0x1a_y0y1a(m_lm, m_right[ir]);
			m_proj -= m_lmr;
		}
	}
	return m_proj;
}
RI::Tensor<double> Exx_Opt_Orb::cal_proj_21(
	const RI::Tensor<double> & m_big,
	const std::vector<RI::Tensor<double>> & m_left,
	const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
	const std::vector<RI::Tensor<double>> & m_right )
{
	ModuleBase::TITLE("Exx_Opt_Orb::cal_proj");
	RI::Tensor<double> m_proj = m_big.copy();
	for( size_t il=0; il!=m_left.size(); ++il )
	{
		for( size_t ir=0; ir!=m_right.size(); ++ir )
		{
			// m_proj = m_proj - m_left[il] * m_middle[il][ir] * m_right[ir].T;
			const RI::Tensor<double> m_lm = RI::Tensor_Multiply::x0x1y1_x0x1a_ay1(m_left[il], m_middle[il][ir]);
			const RI::Tensor<double> m_lmr = RI::Tensor_Multiply::x0x1y0_x0x1a_y0a(m_lm, m_right[ir]);
			m_proj -= m_lmr;
		}
	}
	return m_proj;
}
/*RI::Tensor<double> Exx_Opt_Orb::cal_proj_12(
	const RI::Tensor<double> & m_big,
	const std::vector<RI::Tensor<double>> & m_left,
	const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
	const std::vector<RI::Tensor<double>> & m_right )
{
	ModuleBase::TITLE("Exx_Opt_Orb::cal_proj");
	RI::Tensor<double> m_proj = m_big.copy();
	for( size_t il=0; il!=m_left.size(); ++il )
	{
		for( size_t ir=0; ir!=m_right.size(); ++ir )
		{
			// m_proj = m_proj - m_left[il] * m_middle[il][ir] * m_right[ir].T;
			const RI::Tensor<double> m_lm = RI::Tensor_Multiply::x0y1_x0a_ay1(m_left[il], m_middle[il][ir]);
			const RI::Tensor<double> m_lmr = RI::Tensor_Multiply::x0y0y1_x0a_y0y1a(m_lm, m_right[ir]);
			m_proj -= m_lmr;
		}
	}
	return m_proj;
}*/
RI::Tensor<double> Exx_Opt_Orb::cal_proj_11(
	const RI::Tensor<double> & m_big,
	const std::vector<RI::Tensor<double>> & m_left,
	const std::vector<std::vector<RI::Tensor<double>>> & m_middle,
	const std::vector<RI::Tensor<double>> & m_right )
{
	ModuleBase::TITLE("Exx_Opt_Orb::cal_proj");
	RI::Tensor<double> m_proj = m_big.copy();
	for( size_t il=0; il!=m_left.size(); ++il )
	{
		for( size_t ir=0; ir!=m_right.size(); ++ir )
		{
			// m_proj = m_proj - m_left[il] * m_middle[il][ir] * m_right[ir].T;
			const RI::Tensor<double> m_lm = RI::Tensor_Multiply::x0y1_x0a_ay1(m_left[il], m_middle[il][ir]);
			const RI::Tensor<double> m_lmr = RI::Tensor_Multiply::x0y0_x0a_y0a(m_lm, m_right[ir]);
			m_proj -= m_lmr;
		}
	}
	return m_proj;
}

std::vector<std::vector<RI::Tensor<double>>> Exx_Opt_Orb::cal_I(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,RI::Tensor<double>>>>> &ms,
	const size_t TA, const size_t IA, const size_t TB, const size_t IB )
{
	ModuleBase::TITLE("Exx_Opt_Orb::cal_I");

	if( TA==TB && IA==IB )
	{
		return {{LRI_CV_Tools::cal_I(RI::Tensor<double>(ms.at(TA).at(IA).at(TA).at(IA)))}};
	}
	else
	{
		std::vector<std::vector<RI::Tensor<double>>> m_in
			{{ ms.at(TA).at(IA).at(TA).at(IA),
			   ms.at(TA).at(IA).at(TB).at(IB) },
			 { ms.at(TB).at(IB).at(TA).at(IA),
			   ms.at(TB).at(IB).at(TB).at(IB) }};
		return LRI_CV_Tools::cal_I(m_in);
	}
}

std::map<size_t,std::map<size_t,std::set<double>>> Exx_Opt_Orb::get_radial_R()
{
	ModuleBase::TITLE("Exx_Opt_Orb::get_radial_R");
	std::map<size_t,std::map<size_t,std::set<double>>> radial_R;
	for( size_t TA=0; TA!=GlobalC::ucell.ntype; ++TA )
	{
		for( size_t IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA )
		{
			for( size_t TB=0; TB!=GlobalC::ucell.ntype; ++TB )
			{
				for( size_t IB=0; IB!=GlobalC::ucell.atoms[TB].na; ++IB )
				{
					const ModuleBase::Vector3<double> &tauA = GlobalC::ucell.atoms[TA].tau[IA];
					const ModuleBase::Vector3<double> &tauB = GlobalC::ucell.atoms[TB].tau[IB];
					const double delta_R = (-tauA+tauB).norm();
					radial_R[TA][TB].insert( delta_R );
					radial_R[TB][TA].insert( delta_R );
				}
			}
		}
	}
	return radial_R;
}
