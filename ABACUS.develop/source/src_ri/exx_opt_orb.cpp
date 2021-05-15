#include "exx_opt_orb.h"
#include "src_pw/global.h"
#include "module_ORB/ORB_atomic_lm.h"
#include "src_ri/exx_abfs.h"
#include "src_ri/exx_abfs-abfs_index.h"
#include "src_ri/exx_abfs-construct_orbs.h"
#include "src_ri/exx_abfs-inverse_matrix_double.h"
#include "src_ri/exx_abfs-jle.h"
#include "module_ORB/ORB_read.h"
#include "src_ri/exx_abfs-matrix_orbs11.h"
#include "src_ri/exx_abfs-matrix_orbs21.h"
#include "src_ri/exx_abfs-matrix_orbs22.h"

#include "src_external/src_test/src_ri/exx_abfs-unittest.h"
#include "src_external/src_test/src_ri/exx_lcao-test.h"
#include "src_external/src_test/src_global/element_basis_index-test.h"
#include "src_external/src_test/test_function.h"

void Exx_Opt_Orb::generate_matrix() const
{
ofstream ofs_mpi(exx_lcao.test_dir.process+"time_"+TO_STRING(MY_RANK),ofstream::app);

	TITLE("Exx_Opt_Orb::generate_matrix");
ofs_mpi<<"memory:\t"<<get_memory(10)<<endl;

	const vector<vector<vector<Numerical_Orbital_Lm>>> lcaos = Exx_Abfs::Construct_Orbs::change_orbs( ORB, this->kmesh_times );

	const vector<vector<vector<Numerical_Orbital_Lm>>> abfs = Exx_Abfs::Construct_Orbs::abfs_same_atom( lcaos, this->kmesh_times, exx_lcao.info.pca_threshold );

ofs_mpi<<"memory:\t"<<get_memory(10)<<endl;
	
	Exx_Abfs::Jle jle;
	jle.init_jle( this->kmesh_times );

ofs_mpi<<"memory:\t"<<get_memory(10)<<endl;
	
	Exx_Abfs::Lmax = Exx_Abfs::Jle::Lmax;
	for( size_t T=0; T!=abfs.size(); ++T )
		Exx_Abfs::Lmax = std::max( Exx_Abfs::Lmax, static_cast<int>(abfs[T].size())-1 );

	const Element_Basis_Index::Range    range_lcaos = Exx_Abfs::Abfs_Index::construct_range( lcaos );
	const Element_Basis_Index::IndexLNM index_lcaos = Element_Basis_Index::construct_index( range_lcaos );

	const Element_Basis_Index::Range    range_abfs = Exx_Abfs::Abfs_Index::construct_range( abfs );
	const Element_Basis_Index::IndexLNM index_abfs = Element_Basis_Index::construct_index( range_abfs );

	const Element_Basis_Index::Range    range_jys = Exx_Abfs::Abfs_Index::construct_range( jle.jle );
	const Element_Basis_Index::IndexLNM index_jys = Element_Basis_Index::construct_index( range_jys );

ofs_mpi<<range_lcaos<<endl;
ofs_mpi<<range_abfs<<endl;
ofs_mpi<<range_jys<<endl;

	map<size_t,map<size_t,set<double>>> radial_R = get_radial_R();
#if TEST_EXX_RADIAL==2
	{
		for(const auto & rA : radial_R)
			for(const auto & rB : rA.second)
			{
				ofs_mpi<<rA.first<<"\t"<<rB.first<<":\t";
				for(const auto & rC : rB.second)
					ofs_mpi<<rC<<"\t";
				ofs_mpi<<endl;
			}
	}
#endif

ofs_mpi<<"memory:\t"<<get_memory(10)<<endl;

	// < lcaos lcaos | lcaos lcaos >
	const auto ms_lcaoslcaos_lcaoslcaos = [&]() -> map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> 
	{
		Exx_Abfs::Matrix_Orbs22 m_lcaoslcaos_lcaoslcaos;
		m_lcaoslcaos_lcaoslcaos.init( 1, this->kmesh_times, 1 );
		m_lcaoslcaos_lcaoslcaos.init_radial( lcaos, lcaos, lcaos, lcaos );
		#if TEST_EXX_RADIAL>=1
		m_lcaoslcaos_lcaoslcaos.init_radial_table(radial_R);
		#else
		m_lcaoslcaos_lcaoslcaos.init_radial_table();
		#endif
		return m_lcaoslcaos_lcaoslcaos.cal_overlap_matrix( index_lcaos, index_lcaos, index_lcaos, index_lcaos);
	}();
	
ofs_mpi<<"memory:\t"<<get_memory(10)<<endl;

	// < lcaos lcaos | jys >
	const auto ms_lcaoslcaos_jys = [&]() -> map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>>
	{
		Exx_Abfs::Matrix_Orbs21 m_jyslcaos_lcaos;
		m_jyslcaos_lcaos.init( 1, this->kmesh_times, 1 );
		m_jyslcaos_lcaos.init_radial( jle.jle, lcaos, lcaos );
		#if TEST_EXX_RADIAL>=1
		m_jyslcaos_lcaos.init_radial_table(radial_R);
		#else
		m_jyslcaos_lcaos.init_radial_table();
		#endif
		return m_jyslcaos_lcaos.cal_overlap_matrix( index_jys, index_lcaos, index_lcaos );
	}();

ofs_mpi<<"memory:\t"<<get_memory(10)<<endl;

	// < jys | jys >
	const auto ms_jys_jys = [&]() -> map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
	{
		Exx_Abfs::Matrix_Orbs11 m_jys_jys;
		m_jys_jys.init( 2, this->kmesh_times, 1 );
		m_jys_jys.init_radial( jle.jle, jle.jle );
		#if TEST_EXX_RADIAL>=1
		m_jys_jys.init_radial_table(radial_R);
		#else
		m_jys_jys.init_radial_table();
		#endif
		return m_jys_jys.cal_overlap_matrix( index_jys, index_jys );
	}();

ofs_mpi<<"memory:\t"<<get_memory(10)<<endl;

	// < abfs | abfs >
	const auto ms_abfs_abfs = [&]() -> map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
	{
		Exx_Abfs::Matrix_Orbs11 m_abfs_abfs;
		m_abfs_abfs.init( 2, this->kmesh_times, 1 );
		m_abfs_abfs.init_radial( abfs, abfs );
		#if TEST_EXX_RADIAL>=1
		m_abfs_abfs.init_radial_table(radial_R);
		#else
		m_abfs_abfs.init_radial_table();
		#endif
		return m_abfs_abfs.cal_overlap_matrix( index_abfs, index_abfs );
	}();

ofs_mpi<<"memory:\t"<<get_memory(10)<<endl;

	// < lcaos lcaos | abfs >
	const auto ms_lcaoslcaos_abfs = [&]() -> map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>>
	{
		Exx_Abfs::Matrix_Orbs21 m_abfslcaos_lcaos;
		m_abfslcaos_lcaos.init( 1, this->kmesh_times, 1 );
		m_abfslcaos_lcaos.init_radial( abfs, lcaos, lcaos );
		#if TEST_EXX_RADIAL>=1
		m_abfslcaos_lcaos.init_radial_table(radial_R);
		#else
		m_abfslcaos_lcaos.init_radial_table();
		#endif
		return m_abfslcaos_lcaos.cal_overlap_matrix( index_abfs, index_lcaos, index_lcaos );
	}();

ofs_mpi<<"memory:\t"<<get_memory(10)<<endl;

	// < jys | abfs >
	const auto ms_jys_abfs = [&]() -> map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
	{
		Exx_Abfs::Matrix_Orbs11 m_jys_abfs;
		m_jys_abfs.init( 2, this->kmesh_times, 1 );
		m_jys_abfs.init_radial( jle.jle, abfs );
		#if TEST_EXX_RADIAL>=1
		m_jys_abfs.init_radial_table(radial_R);
		#else
		m_jys_abfs.init_radial_table();
		#endif
		return m_jys_abfs.cal_overlap_matrix( index_jys, index_abfs );
	}();

ofs_mpi<<"memory:\t"<<get_memory(10)<<endl;

ofs_matrixes(exx_lcao.test_dir.matrix+"ms_jys_jys",ms_jys_jys);
ofs_matrixes(exx_lcao.test_dir.matrix+"ms_lcaoslcaos_jys",ms_lcaoslcaos_jys);
ofs_matrixes(exx_lcao.test_dir.matrix+"ms_lcaoslcaos_lcaoslcaos",ms_lcaoslcaos_lcaoslcaos);
ofs_matrixes(exx_lcao.test_dir.matrix+"ms_abfs_abfs",ms_abfs_abfs);
ofs_matrixes(exx_lcao.test_dir.matrix+"ms_lcaoslcaos_abfs",ms_lcaoslcaos_abfs);
ofs_matrixes(exx_lcao.test_dir.matrix+"ms_jys_abfs",ms_jys_abfs);

	for( size_t TA=0; TA!=ucell.ntype; ++TA )
	{
		for( size_t IA=0; IA!=ucell.atoms[TA].na; ++IA )
		{
			for( size_t TB=0; TB!=ucell.ntype; ++TB )
			{
				for( size_t IB=0; IB!=ucell.atoms[TB].na; ++IB )
				{
					if( TA==TB && IA==IB )
					{
						const size_t T=TA, I=IA;
						if(exx_lcao.info.pca_threshold<=1)
						{
							// < abfs | abfs >.I
							const vector<vector<matrix>> ms_abfs_abfs_I = cal_I( ms_abfs_abfs, T,I,T,I );
							// < lcaos lcaos | lcaos lcaos > - < lcaos lcaos | abfs > * < abfs | abfs >.I * < abfs | lcaos lcaos >
							const matrix m_lcaoslcaos_lcaoslcaos_proj =
								cal_proj(
									ms_lcaoslcaos_lcaoslcaos.at(T).at(I).at(T).at(I),
									ms_lcaoslcaos_abfs.at(T).at(I).at(T).at(I),
									ms_abfs_abfs_I,
									ms_lcaoslcaos_abfs.at(T).at(I).at(T).at(I));
							// < lcaos lcaos | jys > - < lcaos lcaos | abfs > * < abfs | abfs >.I * < abfs | jys >
							const vector<matrix> m_lcaoslcaos_jys_proj =
								{cal_proj(
									ms_lcaoslcaos_jys.at(T).at(I).at(T).at(I)[0],
									ms_lcaoslcaos_abfs.at(T).at(I).at(T).at(I),
									ms_abfs_abfs_I,
									{ms_jys_abfs.at(T).at(I).at(T).at(I)})};
							// < jys | jys > - < jys | abfs > * < abfs | abfs >.I * < abfs | jys >
							const vector<vector<matrix>> m_jys_jys_proj =
								{{cal_proj(
									ms_jys_jys.at(T).at(I).at(T).at(I),
									{ms_jys_abfs.at(T).at(I).at(T).at(I)},
									ms_abfs_abfs_I,
									{ms_jys_abfs.at(T).at(I).at(T).at(I)})}};
							print_matrix(
								"matrix",
								m_lcaoslcaos_jys_proj,
								m_jys_jys_proj,
								m_lcaoslcaos_lcaoslcaos_proj,
								TA, IA, TB, IB,
								range_jys, index_jys );
						}
						else
						{
							print_matrix(
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
						if(exx_lcao.info.pca_threshold<=1)
						{
							// < abfs | abfs >.I
							const vector<vector<matrix>> ms_abfs_abfs_I = cal_I( ms_abfs_abfs, TA,IA,TB,IB );
							// < lcaos lcaos | lcaos lcaos > - < lcaos lcaos | abfs > * < abfs | abfs >.I * < abfs | lcaos lcaos >
							const matrix m_lcaoslcaos_lcaoslcaos_proj =
								cal_proj(
									ms_lcaoslcaos_lcaoslcaos.at(TA).at(IA).at(TB).at(IB),
									ms_lcaoslcaos_abfs.at(TA).at(IA).at(TB).at(IB),
									ms_abfs_abfs_I,
									ms_lcaoslcaos_abfs.at(TA).at(IA).at(TB).at(IB));
							// < lcaos lcaos | jys > - < lcaos lcaos | abfs > * < abfs | abfs >.I * < abfs | jys >
							const vector<matrix> m_lcaoslcaos_jys_proj =
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
							const vector<vector<matrix>> m_jys_jys_proj =
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
							print_matrix(
								"matrix",
								m_lcaoslcaos_jys_proj,
								m_jys_jys_proj,
								m_lcaoslcaos_lcaoslcaos_proj,
								TA, IA, TB, IB,
								range_jys, index_jys );
						}
						else
						{
							print_matrix(
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
matrix Exx_Opt_Orb::cal_proj( 
	const matrix & m_big, 
	const vector<matrix> & m_left, 
	const vector<vector<matrix>> & m_middle, 
	const vector<matrix> & m_right ) const
{
	TITLE("Exx_Opt_Orb::cal_proj");

//auto print_nrc = [](const matrix & m){ cout<<"\t"<<m.nr<<"\t"<<m.nc<<endl; };

	matrix m_proj = m_big;
//print_nrc(m_proj);
	for( size_t il=0; il!=m_left.size(); ++il )
	{
		for( size_t ir=0; ir!=m_right.size(); ++ir )
		{
//cout<<il<<"\t"<<ir<<endl;
//print_nrc(m_left[il]);
//print_nrc(m_middle[il][ir]);
//print_nrc(m_right[ir]);
			m_proj -= m_left[il] * m_middle[il][ir] * transpose(m_right[ir]);
		}
	}
	return m_proj;
}

vector<vector<matrix>> Exx_Opt_Orb::cal_I(
	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &ms,
	const size_t TA, const size_t IA, const size_t TB, const size_t IB ) const
{
	TITLE("Exx_Opt_Orb::cal_I");

	if( TA==TB && IA==IB )
	{
		vector<vector<matrix>> m_I
			{{ {ms.at(TA).at(IA).at(TA).at(IA).nr, ms.at(TA).at(IA).at(TA).at(IA).nc} }};
		Exx_Abfs::Inverse_Matrix_Double inv;
		inv.init( ms.at(TA).at(IA).at(TA).at(IA).nr );
		inv.input( ms.at(TA).at(IA).at(TA).at(IA) );
		inv.cal_inverse( Exx_Abfs::Inverse_Matrix_Double::Method::dsyev );
		inv.output( m_I[0][0] );
		return m_I;
	}
	else
	{
		vector<vector<matrix>> m_I
			{{ {ms.at(TA).at(IA).at(TA).at(IA).nr, ms.at(TA).at(IA).at(TA).at(IA).nc},
			   {ms.at(TA).at(IA).at(TB).at(IB).nr, ms.at(TA).at(IA).at(TB).at(IB).nc} },
			 { {ms.at(TB).at(IB).at(TA).at(IA).nr, ms.at(TB).at(IB).at(TA).at(IA).nc},
			   {ms.at(TB).at(IB).at(TB).at(IB).nr, ms.at(TB).at(IB).at(TB).at(IB).nc} }};
		Exx_Abfs::Inverse_Matrix_Double inv;
		inv.init( ms.at(TA).at(IA).at(TA).at(IA).nr + ms.at(TB).at(IB).at(TA).at(IA).nr );
		inv.input( ms.at(TA).at(IA).at(TA).at(IA), ms.at(TA).at(IA).at(TB).at(IB), ms.at(TB).at(IB).at(TA).at(IA), ms.at(TB).at(IB).at(TB).at(IB) );
		inv.cal_inverse( Exx_Abfs::Inverse_Matrix_Double::Method::dsyev );
		inv.output( m_I[0][0], m_I[0][1], m_I[1][0], m_I[1][1] );
		return m_I;
	}
}

map<size_t,map<size_t,set<double>>> Exx_Opt_Orb::get_radial_R() const
{
	TITLE("Exx_Opt_Orb::get_radial_R");
	map<size_t,map<size_t,set<double>>> radial_R;
	for( size_t TA=0; TA!=ucell.ntype; ++TA )
		for( size_t IA=0; IA!=ucell.atoms[TA].na; ++IA )
			for( size_t TB=0; TB!=ucell.ntype; ++TB )
				for( size_t IB=0; IB!=ucell.atoms[TB].na; ++IB )
				{
					const Vector3<double> &tauA = ucell.atoms[TA].tau[IA];
					const Vector3<double> &tauB = ucell.atoms[TB].tau[IB];
					const double delta_R = (-tauA+tauB).norm();
					radial_R[TA][TB].insert( delta_R );
					radial_R[TB][TA].insert( delta_R );
				}
	return radial_R;
}