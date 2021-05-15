#include "exx_abfs.h"

#include "exx_abfs-abfs_index.h"
#include "exx_abfs-jle.h"
#include "exx_abfs-inverse_matrix_double.h"
#include "exx_abfs-io.h"
#include "exx_abfs-construct_orbs.h"

#include "exx_abfs-matrix_orbs11.h"
#include "exx_abfs-matrix_orbs21.h"
#include "exx_abfs-matrix_lcaoslcaos_lcaoslcaos.h"

#include "module_ORB/ORB_read.h"
#include "conv_coulomb_pot.h"
#include "conv_coulomb_pot-inl.h"

#include "src_global/global_function.h"

#include "src_external/src_test/test_function.h"				// Peize Lin test
#include "src_external/src_test/src_ri/exx_abfs-unittest.h"
#include "src_external/src_test/src_ri/make_gaunt_table-unittest.h"
#include "src_external/src_test/src_global/element_basis_index-test.h"			// Peize Lin test 2016-04-05
#include "src_pw/global.h"
#include<sys/time.h>				// Peize Lin test

int Exx_Abfs::Lmax = 0;		// Peize Lin test

void Exx_Abfs::test_all() const
{
	auto test_MGT = []()
	{
		const int Lmax = 3;
		ORB_gaunt_table MGT;
		MGT.init_Gaunt_CH( Lmax );
		MGT.init_Gaunt( Lmax );
		cout_MGT(MGT,Lmax);
	};

	auto test_abfs = []()
	{
		const vector<vector<vector<Numerical_Orbital_Lm>>>
			&&lcaos = Construct_Orbs::change_orbs( ORB, 1 );
		const vector<vector<vector<Numerical_Orbital_Lm>>>
			&&abfs = Construct_Orbs::abfs_same_atom( lcaos, 1 );

		for( size_t T=0; T!=abfs.size(); ++T )
			Lmax = std::max( Lmax, static_cast<int>(abfs[T].size())-1 );

		const Element_Basis_Index::Range
			&&range_abfs = Abfs_Index::construct_range( abfs );
		const Element_Basis_Index::IndexLNM
			&&index_abfs = Element_Basis_Index::construct_index( range_abfs );

		Matrix_Orbs11 m_abfs_abfs;
		cout<<"D1"<<endl;
		m_abfs_abfs.init( 2, 1, 1 );
		cout<<"D2"<<endl;
		m_abfs_abfs.init_radial( abfs, abfs );
		cout<<"D3"<<endl;
		m_abfs_abfs.init_radial_table();
		cout<<"D4"<<endl;
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
			&&ms_abfs_abfs = m_abfs_abfs.cal_overlap_matrix( index_abfs, index_abfs );
		ofs_ms("ms_abfs_abfs",ms_abfs_abfs);
	};

	auto test_svd = []()
	{
		const vector<vector<vector<Numerical_Orbital_Lm>>>
			&&lcaos = Construct_Orbs::change_orbs( ORB, 1 );
		const vector<vector<vector<Numerical_Orbital_Lm>>>
			&&abfs = Construct_Orbs::abfs_same_atom( lcaos, 1 );

		for( size_t T=0; T!=abfs.size(); ++T )
			Lmax = std::max( Lmax, static_cast<int>(abfs[T].size())-1 );

		const Element_Basis_Index::Range
			&&range_abfs = Abfs_Index::construct_range( abfs );
		const Element_Basis_Index::IndexLNM
			&&index_abfs = Element_Basis_Index::construct_index( range_abfs );

		const Element_Basis_Index::Range
			&&range_lcaos = Abfs_Index::construct_range( ORB );
		const Element_Basis_Index::IndexLNM
			&&index_lcaos = Element_Basis_Index::construct_index( range_lcaos );

		Matrix_Orbs21 m_abfslcaos_lcaos;
		m_abfslcaos_lcaos.init( 1, 1, 1 );
		m_abfslcaos_lcaos.init_radial( abfs, ORB, ORB );
		m_abfslcaos_lcaos.init_radial_table();
		const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>>
			&&ms_lcaos2_abfs = m_abfslcaos_lcaos.cal_overlap_matrix( index_abfs, index_lcaos, index_lcaos );
		ofs_ms("ms_lcaos2_abfs",ms_lcaos2_abfs);
	};

	auto test_GriD = []()
	{
		for (int T1 = 0; T1 < ucell.ntype; ++T1)
			for (int I1 = 0; I1 < ucell.atoms[T1].na; ++I1)
			{
				cout<<"@\t"<<T1<<"\t"<<I1<<endl;
				GridD.Find_atom( ucell.atoms[T1].tau[I1], T1, I1 );
				for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
					cout<<GridD.getBox(ad).x<<"\t"<<GridD.getBox(ad).y<<"\t"<<GridD.getBox(ad).z<<endl;
			}
	};

	auto test_k = []()
	{
		for( size_t ik=0; ik!=kv.nks; ++ik )
			cout<<kv.kvec_d[ik].x<<"\t"<<kv.kvec_d[ik].y<<"\t"<<kv.kvec_d[ik].z<<endl;
	};
}

void Exx_Abfs::generate_matrix() const
{

cout<<"A"<<endl;

	Jle jle;
	jle.init_jle(1);

cout<<"A1"<<endl;

//	const vector<vector<vector<Numerical_Orbital_Lm>>>
//		&&abfs_same_atom = Construct_Orbs::abfs_same_atom( ORB );

	Exx_Abfs::Lmax = Jle::Lmax;
//	for( size_t T=0; T!=abfs_same_atom.size(); ++T )
//		Exx_Abfs::Lmax = std::max( Exx_Abfs::Lmax, static_cast<int>(abfs_same_atom[T].size()) );

cout<<Exx_Abfs::Lmax<<endl;

/*
	// Peize Lin test
	ofstream ofsN("N_orbital.dat");
	for( size_t T=0; T!=abfs_same_atom.size(); ++T )
		for( size_t L=0; L!=abfs_same_atom[T].size(); ++L )
			for( size_t N=0; N!=abfs_same_atom[T][L].size(); ++N )
			{
				ofsN<<T<<"\t"<<L<<"\t"<<N<<endl;
				for( size_t ir=0; ir!=abfs_same_atom[T][L][N].getNr(); ++ir )
					ofsN<<abfs_same_atom[T][L][N].getPsi(ir)<<endl;
				ofsN<<endl;
			}
*/

cout<<"B"<<endl;

	Matrix_Orbs11 m_jys_jys;
	m_jys_jys.init(2);
	m_jys_jys.init_radial( jle.jle, jle.jle );
	m_jys_jys.init_radial_table();

cout<<"C"<<endl;

	Matrix_Orbs21 m_jyslcaos_lcaos;
	m_jyslcaos_lcaos.init(1);
	m_jyslcaos_lcaos.init_radial( jle.jle, ORB, ORB );
	m_jyslcaos_lcaos.init_radial_table();

cout<<"D"<<endl;

	Matrix_Lcaoslcaos_Lcaoslcaos mllll;
	mllll.init(1);
	mllll.init_radial( ORB, ORB );
	mllll.init_radial_table();

cout<<"D1"<<endl;

//	Matrix_Orbs21 m_asalcaos_lcaos;					// "asa" means "abfs_same_atom"
//cout<<"D11"<<endl;
//	m_asalcaos_lcaos.init(1);
//cout<<"D12"<<endl;
//	m_asalcaos_lcaos.init_radial( abfs_same_atom, ORB, ORB );
//cout<<"D13"<<endl;
//	m_asalcaos_lcaos.init_radial_table();

cout<<"D2"<<endl;

//	Matrix_Orbs11 m_asa_asa;
//cout<<"D21"<<endl;
//	m_asa_asa.init(2);
//cout<<"D22"<<endl;
//	m_asa_asa.init_radial( abfs_same_atom, abfs_same_atom );
//cout<<"D23"<<endl;
//	m_asa_asa.init_radial_table();

cout<<"D3"<<endl;

//	Matrix_Orbs11 m_asa_jys;
//cout<<"D31"<<endl;
//	m_asa_jys.init(2);
//cout<<"D32"<<endl;
//	m_asa_jys.init_radial( abfs_same_atom, jle.jle );
//cout<<"D33"<<endl;
//	m_asa_jys.init_radial_table();

cout<<"E"<<endl;

	const Element_Basis_Index::Range
		&&range_lcaos = Abfs_Index::construct_range( ORB );
	const Element_Basis_Index::IndexLNM
		&&index_lcaos = Element_Basis_Index::construct_index( range_lcaos );

cout<<"F"<<endl;

	const Element_Basis_Index::Range
		&&range_jys = Abfs_Index::construct_range( jle.jle );
	const Element_Basis_Index::IndexLNM
		&&index_jys = Element_Basis_Index::construct_index( range_jys );

//	const Element_Basis_Index::Range
//		&&range_asa = Abfs_Index::construct_range( abfs_same_atom );
//	const Element_Basis_Index::Index
//		&&index_asa = Element_Basis_Index::construct_index( range_asa );

cout<<"G"<<endl;

	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
		&&ms_jys_jys = m_jys_jys.cal_overlap_matrix( index_jys, index_jys );

ofs_ms("ms_jys_jys",ms_jys_jys);

	map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>>
		&&ms_lcaos2_jys = m_jyslcaos_lcaos.cal_overlap_matrix( index_jys, index_lcaos, index_lcaos );

ofs_ms("ms_lcaos2_jys",ms_lcaos2_jys);

cout<<"G2"<<endl;

	map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
		&&ms_lcaos2_lcaos2 = mllll.cal_overlap_matrix( index_lcaos, index_lcaos );

ofs_ms("ms_lcaos2_lcaos2",ms_lcaos2_lcaos2);

cout<<"G3"<<endl;

//	const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>>
//		&&ms_lcaos2_asa = m_asalcaos_lcaos.cal_overlap_matrix( index_asa, index_lcaos, index_lcaos );

//ofs_ms("ms_lcaos2_asa",ms_lcaos2_asa);

cout<<"G4"<<endl;

//	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
//		&&ms_asa_asa = m_asa_asa.cal_overlap_matrix( index_asa, index_asa );

//ofs_ms("ms_asa_asa",ms_asa_asa);

cout<<"G5"<<endl;

//	const map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>>
//		&&ms_asa_asa_I = cal_I( ms_asa_asa, index_asa );

//ofs_ms("ms_asa_asa_I",ms_asa_asa_I);

cout<<"G6"<<endl;

//	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
//		&&ms_asa_jys = m_asa_jys.cal_overlap_matrix( index_asa, index_jys );

//ofs_ms("ms_asa_jys",ms_asa_jys);

//	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
//		&&ms_lcaos2_lcaos2_proj_asa = cal_lcaos2_lcaos2_proj_asa( ms_lcaos2_asa, ms_asa_asa_I, range_lcaos, index_lcaos );

//ofs_ms("ms_lcaos2_lcaos2_proj_asa",ms_lcaos2_lcaos2_proj_asa);

//	const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>>
//		&&ms_lcaos2_jys_proj_asa = cal_lcaos2_jys_proj_asa( ms_lcaos2_asa, ms_asa_asa_I, ms_asa_jys );

//ofs_ms("ms_lcaos2_jys_proj_asa",ms_lcaos2_jys_proj_asa);

cout<<"G7"<<endl;

//	std::function< void( matrix &, const matrix & ) >
//		minus_matrix = []( matrix &mA, const matrix &mB ){ mA-=mB; };
//	// ms_lcaos2_lcaos2 -= ms_lcaos2_lcaos2_proj_asa
//	FUNC_EACH_2( ms_lcaos2_lcaos2, ms_lcaos2_lcaos2_proj_asa, minus_matrix );
//	// ms_lcaos2_jys -= ms_lcaos2_jys_proj_asa
//	FUNC_EACH_2( ms_lcaos2_jys, ms_lcaos2_jys_proj_asa, minus_matrix );

//ofs_ms("ms_lcaos2_lcaos2_new",ms_lcaos2_lcaos2);
//ofs_ms("ms_lcaos2_jys_new",ms_lcaos2_jys);

cout<<"H"<<endl;

	const string file_name_prefix = "";			// Peize Lin test
	IO::print_matrix(
		file_name_prefix,
		ms_lcaos2_jys,
		ms_jys_jys,
		ms_lcaos2_lcaos2,
		range_jys, index_jys,
		range_lcaos, index_lcaos );

cout<<"I"<<endl;

}


void Exx_Abfs::test_abfs1() const
{

cout<<"A"<<endl;

cout<<"A1"<<endl;

//for(const auto &f : this->files_abfs)
//	cout<<f<<endl;

	const vector<vector<vector<Numerical_Orbital_Lm>>>
		&&lcaos = Construct_Orbs::change_orbs( ORB, this->kmesh_times );
	const vector<vector<vector<Numerical_Orbital_Lm>>>
		&&abfs_same_atom = Construct_Orbs::abfs_same_atom( lcaos, this->kmesh_times );
	// Peize Lin test
//	vector<vector<vector<Numerical_Orbital_Lm>>>
//		&&abfs_same_atom = IO::Abfs_Same_Atom::construct( ORB );
//	for( size_t T=0; T!=abfs_same_atom.size(); ++T )
//		if( abfs_same_atom[T].size()>4 )
//			abfs_same_atom[T].resize(4);

const Element_Basis_Index::Range
	&&range_asa = Abfs_Index::construct_range( abfs_same_atom );
cout<<range_asa<<endl;
cout<<__FILE__<<__LINE__<<endl;
throw exception();

	const vector<vector<vector<Numerical_Orbital_Lm>>>
//		&&abfs = IO::construct_abfs( ORB, this->files_abfs, this->kmesh_times );						// Peize Lin test
//		&&abfs = IO::construct_abfs( abfs_same_atom, ORB, this->files_abfs, this->kmesh_times );						// Peize Lin test
		&abfs = abfs_same_atom;

	for( size_t T=0; T!=abfs.size(); ++T )
		this->Lmax = std::max( this->Lmax, static_cast<int>(abfs[T].size())-1 );

cout<<Exx_Abfs::Lmax<<endl;

	// Peize Lin test
	auto ofs_N_orbital = [&]()
	{
		ofstream ofsN("N_orbital.dat");
		for( size_t T=0; T!=abfs_same_atom.size(); ++T )
			for( size_t L=0; L!=abfs_same_atom[T].size(); ++L )
				for( size_t N=0; N!=abfs_same_atom[T][L].size(); ++N )
				{
					ofsN<<T<<"\t"<<L<<"\t"<<N<<endl;
					for( size_t ir=0; ir!=abfs_same_atom[T][L][N].getNr(); ++ir )
						ofsN<<abfs_same_atom[T][L][N].getPsi(ir)<<endl;
					ofsN<<endl;
				}
		ofsN.close();
	};

cout<<"B"<<endl;

cout<<"C"<<endl;

cout<<"D"<<endl;

	Matrix_Lcaoslcaos_Lcaoslcaos mllll;
	mllll.init(1);
	mllll.init_radial( ORB, ORB );
	mllll.init_radial_table();

cout<<"D1"<<endl;

	Matrix_Orbs21 m_abfslcaos_lcaos;					// "asa" means "abfs_same_atom"
cout<<"D11"<<endl;
	m_abfslcaos_lcaos.init(1,this->kmesh_times);
cout<<"D12"<<endl;
	m_abfslcaos_lcaos.init_radial( abfs, lcaos, lcaos );
cout<<"D13"<<endl;
	m_abfslcaos_lcaos.init_radial_table();

cout<<"D2"<<endl;

	Matrix_Orbs11 m_abfs_abfs;
cout<<"D21"<<endl;
	m_abfs_abfs.init(2,this->kmesh_times);
cout<<"D22"<<endl;
	m_abfs_abfs.init_radial( abfs, abfs );
cout<<"D23"<<endl;
	m_abfs_abfs.init_radial_table();

cout<<"D3"<<endl;

cout<<"E"<<endl;

	const Element_Basis_Index::Range
		&&range_lcaos = Abfs_Index::construct_range( ORB );
	const Element_Basis_Index::IndexLNM
		&&index_lcaos = Element_Basis_Index::construct_index( range_lcaos );

cout<<"F"<<endl;

	const Element_Basis_Index::Range
		&&range_abfs = Abfs_Index::construct_range( abfs );
	const Element_Basis_Index::IndexLNM
		&&index_abfs = Element_Basis_Index::construct_index( range_abfs );

//cout<<range_asa<<endl;
//cout<<index_asa<<endl;

cout<<"G"<<endl;

cout<<"G2"<<endl;

	map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
		&&ms_lcaos2_lcaos2 = mllll.cal_overlap_matrix( index_lcaos, index_lcaos );

ofs_ms("ms_lcaos2_lcaos2",ms_lcaos2_lcaos2);

cout<<"G3"<<endl;

	const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>>
		&&ms_lcaos2_abfs = m_abfslcaos_lcaos.cal_overlap_matrix( index_abfs, index_lcaos, index_lcaos );

ofs_ms("ms_lcaos2_abfs",ms_lcaos2_abfs);

cout<<"G4"<<endl;

	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
		&&ms_abfs_abfs = m_abfs_abfs.cal_overlap_matrix( index_abfs, index_abfs );

ofs_ms("ms_abfs_abfs",ms_abfs_abfs);

cout<<"G5"<<endl;

	const map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>>
		&&ms_abfs_abfs_I = cal_I( ms_abfs_abfs, index_abfs );

ofs_ms("ms_abfs_abfs_I",ms_abfs_abfs_I);

cout<<"G6"<<endl;

	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
		&&ms_lcaos2_lcaos2_proj_abfs = cal_lcaos2_lcaos2_proj_asa( ms_lcaos2_abfs, ms_abfs_abfs_I, range_lcaos, index_lcaos );

ofs_ms("ms_lcaos2_lcaos2_proj_abfs",ms_lcaos2_lcaos2_proj_abfs);

cout<<"G7"<<endl;

	auto test_ratio = [&]()
	{
		for(const auto m1 : ms_lcaos2_lcaos2_proj_abfs)
		{
			size_t TA = m1.first;
			for(const auto m2 : m1.second)
			{
				size_t IA = m2.first;
				for(const auto m3 : m2.second)
				{
					size_t TB = m3.first;
					for(const auto m4 : m3.second)
					{
						size_t IB = m4.first;
						matrix m = ms_lcaos2_lcaos2[TA][IA][TB][IB] - m4.second;
						for( int ir=0; ir!=m.nr; ++ir )
							for( size_t ic=0; ic!=m.nc; ++ic )
								m(ir,ic) /= ms_lcaos2_lcaos2[TA][IA][TB][IB](ir,ic);
						cout<<average(m)<<"\t"<<m.max()<<endl;
					}
				}
			}
		}
	};

test_ratio();

	std::function< void( matrix &, const matrix & ) >
		minus_matrix = []( matrix &mA, const matrix &mB ){ mA-=mB; };
	// ms_lcaos2_lcaos2 -= ms_lcaos2_lcaos2_proj_asa
	FUNC_EACH_2( ms_lcaos2_lcaos2, ms_lcaos2_lcaos2_proj_abfs, minus_matrix );

ofs_ms("ms_lcaos2_lcaos2_new",ms_lcaos2_lcaos2);

cout<<"H"<<endl;

for(const auto m1 : ms_lcaos2_lcaos2)
	for(const auto m2 : m1.second)
		for(const auto m3 : m2.second)
			for(const auto m4 : m3.second)
				cout<<average(m4.second)<<"\t"<<m4.second.max()<<endl;

cout<<"I"<<endl;

}

void Exx_Abfs::cal_exx() const
{
	// ȫ����ֻһ��
cout<<"A"<<endl;

	const vector<vector<vector<Numerical_Orbital_Lm>>>
		&&lcaos = Construct_Orbs::change_orbs( ORB, this->kmesh_times );
//		&&lcaos = Construct_Orbs::change_orbs( ORB, 1 );

cout<<"A1"<<endl;

	const vector<vector<vector<Numerical_Orbital_Lm>>>
		&&abfs_same_atom = Construct_Orbs::abfs_same_atom( lcaos, this->kmesh_times );
cout<<"A2"<<endl;
//	const vector<vector<vector<Numerical_Orbital_Lm>>>
//		&&abfs_origin = IO::construct_abfs( abfs_same_atom, ORB, this->files_abfs, this->kmesh_times );		// Peize Lin test
//		&&abfs = IO::construct_abfs( ORB, this->files_abfs, this->kmesh_times );						// Peize Lin test
cout<<"A3"<<endl;
	const vector<vector<vector<Numerical_Orbital_Lm>>>
//		&&abfs = Construct_Orbs::orth_orbs( abfs_origin );		// Peize Lin test
		&&abfs = Construct_Orbs::orth_orbs( abfs_same_atom );		// Peize Lin test
cout<<"A4"<<endl;

	vector<vector<vector<Numerical_Orbital_Lm>>> abfs_ccp;
	Conv_Coulomb_Pot::cal_orbs_ccp( abfs, abfs_ccp, this->rmesh_times, 1 );

	for( size_t T=0; T!=abfs.size(); ++T )
		this->Lmax = std::max( this->Lmax, static_cast<int>(abfs[T].size())-1 );

cout<<"B"<<endl;

	const Element_Basis_Index::Range
		&&range_lcaos = Abfs_Index::construct_range( lcaos );
	const Element_Basis_Index::IndexLNM
		&&index_lcaos = Element_Basis_Index::construct_index( range_lcaos );

cout<<"C"<<endl;

	const Element_Basis_Index::Range
		&&range_abfs = Abfs_Index::construct_range( abfs );
	const Element_Basis_Index::IndexLNM
		&&index_abfs = Element_Basis_Index::construct_index( range_abfs );

cout<<range_abfs<<endl;

cout<<"D"<<endl;

	// Peize Lin test
	auto test_matrix_lcaos_lcaos = [&]()
	{
		Matrix_Orbs11 m_lcaos_lcaos;
		m_lcaos_lcaos.init(1);
		m_lcaos_lcaos.init_radial( ORB, ORB );
		m_lcaos_lcaos.init_radial_table();
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
			&&matrix_V = m_lcaos_lcaos.cal_overlap_matrix(index_lcaos,index_lcaos);

		ofstream ofs("S_matrix.dat");
		for( const auto & m1 : matrix_V )
		{
			const size_t TA = m1.first;
			for( const auto & m2 : m1.second )
			{
				const size_t IA = m2.first;
				for( const auto & m3 : m2.second )
				{
					const size_t TB = m3.first;
					for( const auto & m4 : m3.second )
					{
						const size_t IB = m4.first;
						ofs<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<endl;
						ofs<<m4.second<<endl;
					}
				}
			}
		}
		ofs.close();
		WARNING_QUIT( TO_STRING(__FILE__), TO_STRING(__LINE__) );
	};

	// Peize Lin test
	auto test_matrix_lcaos_lcaos2 = [&]()
	{
		Matrix_Orbs11 m_lcaos_lcaos;
		m_lcaos_lcaos.init(1,this->kmesh_times);
		m_lcaos_lcaos.init_radial( lcaos, lcaos );
		m_lcaos_lcaos.init_radial_table();
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
			&&matrix_V = m_lcaos_lcaos.cal_overlap_matrix(index_lcaos,index_lcaos);

		ofstream ofs("S_matrix.dat");
		for( const auto & m1 : matrix_V )
		{
			const size_t TA = m1.first;
			for( const auto & m2 : m1.second )
			{
				const size_t IA = m2.first;
				for( const auto & m3 : m2.second )
				{
					const size_t TB = m3.first;
					for( const auto & m4 : m3.second )
					{
						const size_t IB = m4.first;
						ofs<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<endl;
						ofs<<m4.second<<endl;
					}
				}
			}
		}
		ofs.close();
		WARNING_QUIT( TO_STRING(__FILE__), "line "+TO_STRING(__LINE__) );
	};

	// Peize Lin test
/*
	auto test_matrix_standard_ll2_ll2 = [&]()
	{
		Matrix_Lcaoslcaos_Lcaoslcaos2<Center2_Orb::Orb22_Ccp> m_ll2_ll2;
cout<<"D0"<<endl;
		m_ll2_ll2.init( 1, this->kmesh_times, this->rmesh_times );
cout<<"D1"<<endl;
		m_ll2_ll2.init_radial( ORB, ORB );
cout<<"D2"<<endl;
		m_ll2_ll2.init_radial_table();
cout<<"D3"<<endl;
		const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
			&&ms_ll2_ll2 = m_ll2_ll2.cal_overlap_matrix(index_lcaos,index_lcaos);
cout<<"D4"<<endl;
		ofs_ms("ms_ll2_ll2",ms_ll2_ll2);
	};
	test_matrix_standard_ll2_ll2();
*/
cout<<"DD"<<endl;

	Matrix_Orbs11 m_abfs_abfs;
cout<<"D1"<<endl;
	m_abfs_abfs.init( 2, this->kmesh_times, this->rmesh_times );
cout<<"D2"<<endl;
	m_abfs_abfs.init_radial( abfs, abfs_ccp );
cout<<"D3"<<endl;
	m_abfs_abfs.init_radial_table();

cout<<"E"<<endl;

	Matrix_Orbs21 m_abfslcaos_lcaos;
	m_abfslcaos_lcaos.init( 1, this->kmesh_times, this->rmesh_times );
	m_abfslcaos_lcaos.init_radial( abfs_ccp, lcaos, lcaos );
	m_abfslcaos_lcaos.init_radial_table();

cout<<"F"<<endl;

	// ÿһ�����Ӳ�
	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>>
		&&ms_abfs_abfs = m_abfs_abfs.cal_overlap_matrix( index_abfs, index_abfs );
ofs_ms("ms_abfs_abfs",ms_abfs_abfs);

	const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>>
		&&ms_lcaos2_abfs = m_abfslcaos_lcaos.cal_overlap_matrix( index_abfs, index_lcaos, index_lcaos );
ofs_ms("ms_lcaos2_abfs",ms_lcaos2_abfs);

cout<<"G"<<endl;

	const map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>>
		&&ms_abfs_abfs_I = cal_I( ms_abfs_abfs, index_abfs );
ofs_ms("ms_abfs_abfs_I",ms_abfs_abfs_I);

cout<<"G1"<<endl;
	const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>>
		&&ms_C = cal_C( ms_lcaos2_abfs, ms_abfs_abfs_I );
ofs_ms("ms_C",ms_C);

cout<<"H"<<endl;

	// ÿһ�����Ӳ�
	timeval t_begin;
	gettimeofday( &t_begin, NULL);

	cal_CVC( ms_C, ms_abfs_abfs );

	timeval t_end;
	gettimeofday( &t_end, NULL);
	cout<<"time:\t"<<(double)(t_end.tv_sec-t_begin.tv_sec) + (double)(t_end.tv_usec-t_begin.tv_usec)/1000000.0<<endl;

cout<<"I"<<endl;

}

void Exx_Abfs::test_abfs2() const{}

// <a|a> .I
// &
// <a|a> <a|b>
// <b|a> <b|b> .I
map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>> Exx_Abfs::cal_I(
	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &ms,
	const Element_Basis_Index::IndexLNM &index )
{
	map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>> ms_I;

	for( const auto &m1 : ms )
	{
		const size_t TA = m1.first;
		for( const auto &m2 : m1.second )
		{
			const size_t IA = m2.first;
			for( const auto &m3 : m2.second )
			{
				const size_t TB = m3.first;
				for( const auto &m4 : m3.second )
				{
					const size_t IB = m4.first;

					Inverse_Matrix_Double ms_tmp;

					if( TA==TB && IA==IB )
					{
						const size_t size_A = index[TA].count_size;

						ms_tmp.init( size_A );
						ms_tmp.input(
							ms.at(TA).at(IA).at(TA).at(IA) );
					}
					else
					{
						const size_t size_A = index[TA].count_size;
						const size_t size_B = index[TB].count_size;

						ms_tmp.init( size_A + size_B );
						ms_tmp.input(
							ms.at(TA).at(IA).at(TA).at(IA),
							ms.at(TA).at(IA).at(TB).at(IB),
							ms.at(TB).at(IB).at(TA).at(IA),
							ms.at(TB).at(IB).at(TB).at(IB));
					}

					// Peize Lin test
//					cout<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<endl;
//					const matrix matrix_origin( ms_tmp.A );
//					cout<<matrix_origin<<endl;

					ms_tmp.cal_inverse( Exx_Abfs::Inverse_Matrix_Double::Method::dpotrf );

					// Peize Lin test
//					cout<<ms_tmp.A<<endl;
//					cout<<ms_tmp.A * matrix_origin <<endl;
//					cout<<matrix_origin * ms_tmp.A<<endl;

					if( TA==TB && IA==IB )
					{
						const size_t size_A = index[TA].count_size;

						ms_I[TA][IA][TB][IB].resize( 1, vector<matrix>(1) );
						ms_I[TA][IA][TB][IB][0][0].create( size_A, size_A );

						ms_tmp.output(
							ms_I[TA][IA][TB][IB][0][0]);
					}
					else
					{
						const size_t size_A = index[TA].count_size;
						const size_t size_B = index[TB].count_size;

						ms_I[TA][IA][TB][IB].resize( 2, vector<matrix>(2) );
						ms_I[TA][IA][TB][IB][0][0].create( size_A, size_A );
						ms_I[TA][IA][TB][IB][0][1].create( size_A, size_B );
						ms_I[TA][IA][TB][IB][1][0].create( size_B, size_A );
						ms_I[TA][IA][TB][IB][1][1].create( size_B, size_B );

						ms_tmp.output(
							ms_I[TA][IA][TB][IB][0][0],
							ms_I[TA][IA][TB][IB][0][1],
							ms_I[TA][IA][TB][IB][1][0],
							ms_I[TA][IA][TB][IB][1][1]);
					}
				}
			}
		}
	}
	return ms_I;
}

// <ij|P> * <P|P>.I
map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> Exx_Abfs::cal_C(
	const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &ms_lcaos2_abfs,
	const map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>> &ms_abfs_abfs_I )
{
	map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> ms_C;

	for( const auto & m1 : ms_lcaos2_abfs )
	{
		const size_t TA = m1.first;
		for( const auto & m2 : m1.second )
		{
			const size_t IA = m2.first;
			for( const auto &m3 : m2.second )
			{
				const size_t TB = m3.first;
				for( const auto &m4 : m3.second )
				{
					const size_t IB = m4.first;

					if( TA==TB && IA==IB )
					{
						ms_C[TA][IA][TB][IB].resize(1);
						ms_C[TA][IA][TB][IB][0].create( m4.second[0].nr, m4.second[0].nc );

						const auto &m_lcaos2_abfs = m4.second;
						const auto &m_abfs_abfs_I = ms_abfs_abfs_I.at(TA).at(IA).at(TB).at(IB);

						ms_C[TA][IA][TB][IB][0] = m_lcaos2_abfs[0] * m_abfs_abfs_I[0][0];
					}
					else
					{
						ms_C[TA][IA][TB][IB].resize(2);
						for( size_t i=0; i<2; ++i )
							ms_C[TA][IA][TB][IB][i].create( m4.second[i].nr, m4.second[i].nc );

						const auto &m_lcaos2_abfs = m4.second;
						const auto &m_abfs_abfs_I = ms_abfs_abfs_I.at(TA).at(IA).at(TB).at(IB);
						ms_C[TA][IA][TB][IB][0] = m_lcaos2_abfs[0] * m_abfs_abfs_I[0][0]
													+ m_lcaos2_abfs[1] * m_abfs_abfs_I[1][0];
						ms_C[TA][IA][TB][IB][1] = m_lcaos2_abfs[0] * m_abfs_abfs_I[0][1]
													+ m_lcaos2_abfs[1] * m_abfs_abfs_I[1][1];
					}
				}
			}
		}
	}
	return ms_C;
}

// (<ij|P>*<P|P>.I) * <P|P> * (<P|P>.I*<P|kl>)
void Exx_Abfs::cal_CVC(
	const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &ms_C,
	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &ms_abfs_abfs ) const
{
	ofstream ofs("ms_CVC");

	for( const auto & m11 : ms_C )
	{
		const size_t TA = m11.first;
		for( const auto & m12 : m11.second )
		{
			const size_t IA = m12.first;
			for( const auto & m13 : m12.second )
			{
				const size_t TB = m13.first;
				for( const auto & m14 : m13.second )
				{
					const size_t IB = m14.first;
					for( const auto & m21 : ms_C )
					{
						const size_t TC = m21.first;
						for( const auto & m22 : m21.second )
						{
							const size_t IC = m22.first;
							for( const auto & m23 : m22.second )
							{
								const size_t TD = m23.first;
								for( const auto & m24 : m23.second )
								{
									const size_t ID = m24.first;

									auto m_00 = [&](){ return m14.second[0] * ms_abfs_abfs.at(TA).at(IA).at(TC).at(IC) * transpose(m24.second[0]); };
									auto m_01 = [&](){ return m14.second[0] * ms_abfs_abfs.at(TA).at(IA).at(TD).at(ID) * transpose(m24.second[1]); };
									auto m_10 = [&](){ return m14.second[1] * ms_abfs_abfs.at(TB).at(IB).at(TC).at(IC) * transpose(m24.second[0]); };
									auto m_11 = [&](){ return m14.second[1] * ms_abfs_abfs.at(TB).at(IB).at(TD).at(ID) * transpose(m24.second[1]); };

									//matrix_CVC[TA][IA][TB][IB]|[TC][IC][TD][ID] =
									matrix mm;		// Peize Lin test

									if( TA==TB && IA==IB )
									{
										if( TC==TD && IC==ID )
										{
											mm = m_00();
										}
										else
										{
											mm = m_00()+m_01();
										}
									}
									else
									{
										if( TC==TD && IC==ID )
										{
											mm = m_00()+m_10();
										}
										else
										{
											mm = m_00()+m_01()+m_10()+m_11();
										}
									}

									// Peize Lin test
									ofs<<IA<<"\t"<<IB<<"\t"<<IC<<"\t"<<ID<<endl;
									ofs<<mm<<endl;
								}
							}
						}
					}
				}
			}
		}
	}
	ofs.close();
}

// <ij|P> * <P|P>.I * <P|ij>
map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> Exx_Abfs::cal_lcaos2_lcaos2_proj_asa(
	const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &ms_lcaos2_asa,
	const map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>> &ms_asa_asa_I,
	const Element_Basis_Index::Range &range,
	const Element_Basis_Index::IndexLNM &index)
{
	map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> ms_lcaos2_lcaos2_proj_asa;
	for( const auto &m1 : ms_lcaos2_asa )
	{
		const size_t TA = m1.first;
		for( const auto & m2 : m1.second )
		{
			const size_t IA = m2.first;
			for( const auto & m3 : m2.second )
			{
				const size_t TB = m3.first;
				for( const auto & m4 : m3.second )
				{
					const size_t IB = m4.first;
					const auto & m_lcaos2_asa = m4.second;
					const auto & m_abfs_abfs_I = ms_asa_asa_I.at(TA).at(IA).at(TB).at(IB);

					vector<matrix> mql(2);
					size_t matrix_num;
					if( TA==TB && IA==IB )
					{
						matrix_num = 1;
						mql[0] = m_lcaos2_asa[0] * m_abfs_abfs_I[0][0];
					}
					else
					{
						matrix_num = 2;
						mql[0] = m_lcaos2_asa[0] * m_abfs_abfs_I[0][0] + m_lcaos2_asa[1] * m_abfs_abfs_I[1][0];
						mql[1] = m_lcaos2_asa[0] * m_abfs_abfs_I[0][1] + m_lcaos2_asa[1] * m_abfs_abfs_I[1][1];
					}

					ms_lcaos2_lcaos2_proj_asa[TA][IA][TB][IB].create(
						index[TA].count_size,
						index[TB].count_size );
					for( size_t LA=0; LA!=range[TA].size(); ++LA )
						for( size_t NA=0; NA!=range[TA][LA].N; ++NA )
							for( size_t MA=0; MA!=range[TA][LA].M; ++MA )
								for( size_t LB=0; LB!=range[TB].size(); ++LB )
									for( size_t NB=0; NB!=range[TB][LB].N; ++NB )
										for( size_t MB=0; MB!=range[TB][LB].M; ++MB )
										{
											double mv_sas = 0.0;
											for( size_t im=0; im!=matrix_num; ++im )
											{
												assert( mql[im].nc==m_lcaos2_asa[im].nc );
												for( size_t ic=0; ic!=mql[im].nc; ++ic )
												{
													mv_sas
													+= mql[im](
														Exx_Abfs::Abfs_Index::get_index_index(
															index,TA,LA,NA,MA,
															index,TB,LB,NB,MB ),
														ic)
													* m_lcaos2_asa[im](
														Exx_Abfs::Abfs_Index::get_index_index(
															index,TA,LA,NA,MA,
															index,TB,LB,NB,MB ),
														ic);
												}
											}
											ms_lcaos2_lcaos2_proj_asa[TA][IA][TB][IB](
												index[TA][LA][NA][MA],
												index[TB][LB][NB][MB] )
											= mv_sas;
										}
				}
			}
		}
	}
	return ms_lcaos2_lcaos2_proj_asa;
}

// <ij|P> * <P|P>.I * <P|jY>
map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> Exx_Abfs::cal_lcaos2_jys_proj_asa(
	const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &ms_lcaos2_asa,
	const map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>> &ms_asa_asa_I,
	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &ms_asa_jys)
{
	map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> ms_lcaos2_jys_proj_asa;

	for( const auto & m1 : ms_lcaos2_asa )
	{
		const size_t TA = m1.first;
		for( const auto & m2 : m1.second )
		{
			const size_t IA = m2.first;
			for( const auto &m3 : m2.second )
			{
				const size_t TB = m3.first;
				for( const auto &m4 : m3.second )
				{
					const size_t IB = m4.first;

					if( TA==TB && IA==IB )
					{
						ms_lcaos2_jys_proj_asa[TA][IA][TB][IB].resize(1);

						ms_lcaos2_jys_proj_asa[TA][IA][TB][IB][0]
						= ms_lcaos2_asa.at(TA).at(IA).at(TB).at(IB)[0]
						* ms_asa_asa_I.at(TA).at(IA).at(TB).at(IB)[0][0]
						* ms_asa_jys.at(TA).at(IA).at(TA).at(IA);
					}
					else
					{
						vector<matrix> ms_tmp(2);
						for( size_t i=0; i<2; ++i )
							ms_tmp[i]
							= ms_lcaos2_asa.at(TA).at(IA).at(TB).at(IB)[0]
							* ms_asa_asa_I.at(TA).at(IA).at(TB).at(IB)[0][i]
							+ ms_lcaos2_asa.at(TA).at(IA).at(TB).at(IB)[1]
							* ms_asa_asa_I.at(TA).at(IA).at(TB).at(IB)[1][i];

						ms_lcaos2_jys_proj_asa[TA][IA][TB][IB].resize(2);

						ms_lcaos2_jys_proj_asa[TA][IA][TB][IB][0]
						= ms_tmp[0]
						* ms_asa_jys.at(TA).at(IA).at(TA).at(IA)
						+ ms_tmp[1]
						* ms_asa_jys.at(TB).at(IB).at(TA).at(IA);

						ms_lcaos2_jys_proj_asa[TA][IA][TB][IB][1]
						= ms_tmp[0]
						* ms_asa_jys.at(TA).at(IA).at(TB).at(IB)
						+ ms_tmp[1]
						* ms_asa_jys.at(TB).at(IB).at(TB).at(IB);
					}
				}
			}
		}
	}
	return ms_lcaos2_jys_proj_asa;
}

/*
void cal_R_supercell()
{
	vector<Vector3_Exx> R_supercell;
	for( size_t x=0; x!=kv.nmp[0]; ++x)
		for( size_t y=0; y!=kv.nmp[1]; ++y )
			for( size_t z=0; z!=kv.nmp[2]; ++z )
				R_supercell.push_back(Vector3_Exx(x,y,z));
}
*/
/*
void density_matrix()
{
	vector<matrix> DM_k(kv.nks);
	for( size_t ik=0; ik!=kv.nks; ++ik )
	{
		for( size_t ib=0; ib!=NBANDS; ++ib )
		{
			for( size_t iw1=0; iw1!=NLOCAL; ++iw1 )
			{
				for( size_t iw2=0; iw2!=NLOCAL; ++iw2 )
				{
					DM_k[ik](iw1,iw2) += wf.wg(ik,ib) * conj(LOWF.WFC_K[ik][ib][iw1]) * LOWF.WFC_K[ik][ib][iw2];
				}
			}
		}
	}

	vector<size_t,vector<matrix>> DM_R( NSPIN, vector<matrix>(R_supercell.size()) );
	for( size_t is=0; is!=NSPIN; ++is )
	{
		const size_t k_start = (NSPIN==1) ? 0 : ((is==0) ? 0 : (kv.nks/2));
		const size_t k_end = (NSPIN==1) ? kv.nks : ((is==0) ? (kv.nks/2) : kv.nks);
		for( size_it iR=0; iR!=R_supercell.size(); ++iR )
		{
			for( size_t ik=k_start; ik!=k_end; ++ik )
			{
				DM_R[is][iR] += exp(TWO_PI*IMAG_UNIT*dot(kv.kvec_d[ik],R_supercell[iR])) * DM[ik];
			}
		}
	}
}
*/
