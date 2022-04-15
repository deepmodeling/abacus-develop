#include "exx_abfs.h"

#include "exx_abfs-abfs_index.h"
#include "exx_abfs-jle.h"
#include "exx_abfs-inverse_matrix_double.h"
#include "exx_abfs-io.h"
#include "exx_abfs-construct_orbs.h"

#include "exx_abfs-matrix_orbs11.h"
#include "exx_abfs-matrix_orbs21.h"
#include "exx_abfs-matrix_lcaoslcaos_lcaoslcaos.h"

#include "../module_orbital/ORB_read.h"
#include "conv_coulomb_pot.h"
#include "conv_coulomb_pot-inl.h"

#include "../module_base/global_function.h"

#include "../src_external/src_test/test_function.h"				// Peize Lin test
#include "../src_external/src_test/src_ri/exx_abfs-unittest.h"
#include "../src_external/src_test/src_ri/make_gaunt_table-unittest.h"
#include "../src_external/src_test/src_global/element_basis_index-test.h"			// Peize Lin test 2016-04-05
#include "../src_pw/global.h"
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
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
			&&lcaos = Construct_Orbs::change_orbs( GlobalC::ORB, 1 );
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
			&&abfs = Construct_Orbs::abfs_same_atom( lcaos, 1 );

		for( size_t T=0; T!=abfs.size(); ++T )
			Lmax = std::max( Lmax, static_cast<int>(abfs[T].size())-1 );

		const ModuleBase::Element_Basis_Index::Range
			&&range_abfs = Abfs_Index::construct_range( abfs );
		const ModuleBase::Element_Basis_Index::IndexLNM
			&&index_abfs = ModuleBase::Element_Basis_Index::construct_index( range_abfs );

		Matrix_Orbs11 m_abfs_abfs;
		std::cout<<"D1"<<std::endl;
		m_abfs_abfs.init( 2, 1, 1 );
		std::cout<<"D2"<<std::endl;
		m_abfs_abfs.init_radial( abfs, abfs );
		std::cout<<"D3"<<std::endl;
		m_abfs_abfs.init_radial_table();
		std::cout<<"D4"<<std::endl;
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
			&&ms_abfs_abfs = m_abfs_abfs.cal_overlap_matrix( index_abfs, index_abfs );
		ofs_ms("ms_abfs_abfs",ms_abfs_abfs);
	};

	auto test_svd = []()
	{
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
			&&lcaos = Construct_Orbs::change_orbs( GlobalC::ORB, 1 );
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
			&&abfs = Construct_Orbs::abfs_same_atom( lcaos, 1 );

		for( size_t T=0; T!=abfs.size(); ++T )
			Lmax = std::max( Lmax, static_cast<int>(abfs[T].size())-1 );

		const ModuleBase::Element_Basis_Index::Range
			&&range_abfs = Abfs_Index::construct_range( abfs );
		const ModuleBase::Element_Basis_Index::IndexLNM
			&&index_abfs = ModuleBase::Element_Basis_Index::construct_index( range_abfs );

		const ModuleBase::Element_Basis_Index::Range
			&&range_lcaos = Abfs_Index::construct_range( GlobalC::ORB );
		const ModuleBase::Element_Basis_Index::IndexLNM
			&&index_lcaos = ModuleBase::Element_Basis_Index::construct_index( range_lcaos );

		Matrix_Orbs21 m_abfslcaos_lcaos;
		m_abfslcaos_lcaos.init( 1, 1, 1 );
		m_abfslcaos_lcaos.init_radial( abfs, GlobalC::ORB, GlobalC::ORB );
		m_abfslcaos_lcaos.init_radial_table();
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>>
			&&ms_lcaos2_abfs = m_abfslcaos_lcaos.cal_overlap_matrix( index_abfs, index_lcaos, index_lcaos );
		ofs_ms("ms_lcaos2_abfs",ms_lcaos2_abfs);
	};

	auto test_GriD = []()
	{
		for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
			for (int I1 = 0; I1 < GlobalC::ucell.atoms[T1].na; ++I1)
			{
				std::cout<<"@\t"<<T1<<"\t"<<I1<<std::endl;
				GlobalC::GridD.Find_atom(GlobalC::ucell,  GlobalC::ucell.atoms[T1].tau[I1], T1, I1 );
				for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
					std::cout<<GlobalC::GridD.getBox(ad).x<<"\t"<<GlobalC::GridD.getBox(ad).y<<"\t"<<GlobalC::GridD.getBox(ad).z<<std::endl;
			}
	};

	auto test_k = []()
	{
		for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
			std::cout<<GlobalC::kv.kvec_d[ik].x<<"\t"<<GlobalC::kv.kvec_d[ik].y<<"\t"<<GlobalC::kv.kvec_d[ik].z<<std::endl;
	};
}

void Exx_Abfs::generate_matrix() const
{

std::cout<<"A"<<std::endl;

	Jle jle;
	jle.init_jle(1);

std::cout<<"A1"<<std::endl;

//	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
//		&&abfs_same_atom = Construct_Orbs::abfs_same_atom( GlobalC::ORB );

	Exx_Abfs::Lmax = Jle::Lmax;
//	for( size_t T=0; T!=abfs_same_atom.size(); ++T )
//		Exx_Abfs::Lmax = std::max( Exx_Abfs::Lmax, static_cast<int>(abfs_same_atom[T].size()) );

std::cout<<Exx_Abfs::Lmax<<std::endl;

/*
	// Peize Lin test
	std::ofstream ofsN("N_orbital.dat");
	for( size_t T=0; T!=abfs_same_atom.size(); ++T )
		for( size_t L=0; L!=abfs_same_atom[T].size(); ++L )
			for( size_t N=0; N!=abfs_same_atom[T][L].size(); ++N )
			{
				ofsN<<T<<"\t"<<L<<"\t"<<N<<std::endl;
				for( size_t ir=0; ir!=abfs_same_atom[T][L][N].getNr(); ++ir )
					ofsN<<abfs_same_atom[T][L][N].getPsi(ir)<<std::endl;
				ofsN<<std::endl;
			}
*/

std::cout<<"B"<<std::endl;

	Matrix_Orbs11 m_jys_jys;
	m_jys_jys.init(2);
	m_jys_jys.init_radial( jle.jle, jle.jle );
	m_jys_jys.init_radial_table();

std::cout<<"C"<<std::endl;

	Matrix_Orbs21 m_jyslcaos_lcaos;
	m_jyslcaos_lcaos.init(1);
	m_jyslcaos_lcaos.init_radial( jle.jle, GlobalC::ORB, GlobalC::ORB );
	m_jyslcaos_lcaos.init_radial_table();

std::cout<<"D"<<std::endl;

	Matrix_Lcaoslcaos_Lcaoslcaos mllll;
	mllll.init(1);
	mllll.init_radial( GlobalC::ORB, GlobalC::ORB );
	mllll.init_radial_table();

std::cout<<"D1"<<std::endl;

//	Matrix_Orbs21 m_asalcaos_lcaos;					// "asa" means "abfs_same_atom"
//std::cout<<"D11"<<std::endl;
//	m_asalcaos_lcaos.init(1);
//std::cout<<"D12"<<std::endl;
//	m_asalcaos_lcaos.init_radial( abfs_same_atom, GlobalC::ORB, GlobalC::ORB );
//std::cout<<"D13"<<std::endl;
//	m_asalcaos_lcaos.init_radial_table();

std::cout<<"D2"<<std::endl;

//	Matrix_Orbs11 m_asa_asa;
//std::cout<<"D21"<<std::endl;
//	m_asa_asa.init(2);
//std::cout<<"D22"<<std::endl;
//	m_asa_asa.init_radial( abfs_same_atom, abfs_same_atom );
//std::cout<<"D23"<<std::endl;
//	m_asa_asa.init_radial_table();

std::cout<<"D3"<<std::endl;

//	Matrix_Orbs11 m_asa_jys;
//std::cout<<"D31"<<std::endl;
//	m_asa_jys.init(2);
//std::cout<<"D32"<<std::endl;
//	m_asa_jys.init_radial( abfs_same_atom, jle.jle );
//std::cout<<"D33"<<std::endl;
//	m_asa_jys.init_radial_table();

std::cout<<"E"<<std::endl;

	const ModuleBase::Element_Basis_Index::Range
		&&range_lcaos = Abfs_Index::construct_range( GlobalC::ORB );
	const ModuleBase::Element_Basis_Index::IndexLNM
		&&index_lcaos = ModuleBase::Element_Basis_Index::construct_index( range_lcaos );

std::cout<<"F"<<std::endl;

	const ModuleBase::Element_Basis_Index::Range
		&&range_jys = Abfs_Index::construct_range( jle.jle );
	const ModuleBase::Element_Basis_Index::IndexLNM
		&&index_jys = ModuleBase::Element_Basis_Index::construct_index( range_jys );

//	const ModuleBase::Element_Basis_Index::Range
//		&&range_asa = Abfs_Index::construct_range( abfs_same_atom );
//	const ModuleBase::Element_Basis_Index::Index
//		&&index_asa = ModuleBase::Element_Basis_Index::construct_index( range_asa );

std::cout<<"G"<<std::endl;

	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
		&&ms_jys_jys = m_jys_jys.cal_overlap_matrix( index_jys, index_jys );

ofs_ms("ms_jys_jys",ms_jys_jys);

	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>>
		&&ms_lcaos2_jys = m_jyslcaos_lcaos.cal_overlap_matrix( index_jys, index_lcaos, index_lcaos );

ofs_ms("ms_lcaos2_jys",ms_lcaos2_jys);

std::cout<<"G2"<<std::endl;

	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
		&&ms_lcaos2_lcaos2 = mllll.cal_overlap_matrix( index_lcaos, index_lcaos );

ofs_ms("ms_lcaos2_lcaos2",ms_lcaos2_lcaos2);

std::cout<<"G3"<<std::endl;

//	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>>
//		&&ms_lcaos2_asa = m_asalcaos_lcaos.cal_overlap_matrix( index_asa, index_lcaos, index_lcaos );

//ofs_ms("ms_lcaos2_asa",ms_lcaos2_asa);

std::cout<<"G4"<<std::endl;

//	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
//		&&ms_asa_asa = m_asa_asa.cal_overlap_matrix( index_asa, index_asa );

//ofs_ms("ms_asa_asa",ms_asa_asa);

std::cout<<"G5"<<std::endl;

//	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>>
//		&&ms_asa_asa_I = cal_I( ms_asa_asa, index_asa );

//ofs_ms("ms_asa_asa_I",ms_asa_asa_I);

std::cout<<"G6"<<std::endl;

//	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
//		&&ms_asa_jys = m_asa_jys.cal_overlap_matrix( index_asa, index_jys );

//ofs_ms("ms_asa_jys",ms_asa_jys);

//	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
//		&&ms_lcaos2_lcaos2_proj_asa = cal_lcaos2_lcaos2_proj_asa( ms_lcaos2_asa, ms_asa_asa_I, range_lcaos, index_lcaos );

//ofs_ms("ms_lcaos2_lcaos2_proj_asa",ms_lcaos2_lcaos2_proj_asa);

//	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>>
//		&&ms_lcaos2_jys_proj_asa = cal_lcaos2_jys_proj_asa( ms_lcaos2_asa, ms_asa_asa_I, ms_asa_jys );

//ofs_ms("ms_lcaos2_jys_proj_asa",ms_lcaos2_jys_proj_asa);

std::cout<<"G7"<<std::endl;

//	std::function< void( matrix &, const matrix & ) >
//		minus_matrix = []( matrix &mA, const matrix &mB ){ mA-=mB; };
//	// ms_lcaos2_lcaos2 -= ms_lcaos2_lcaos2_proj_asa
//	FUNC_EACH_2( ms_lcaos2_lcaos2, ms_lcaos2_lcaos2_proj_asa, minus_matrix );
//	// ms_lcaos2_jys -= ms_lcaos2_jys_proj_asa
//	FUNC_EACH_2( ms_lcaos2_jys, ms_lcaos2_jys_proj_asa, minus_matrix );

//ofs_ms("ms_lcaos2_lcaos2_new",ms_lcaos2_lcaos2);
//ofs_ms("ms_lcaos2_jys_new",ms_lcaos2_jys);

std::cout<<"H"<<std::endl;

	const std::string file_name_prefix = "";			// Peize Lin test
	IO::print_matrix(
		file_name_prefix,
		ms_lcaos2_jys,
		ms_jys_jys,
		ms_lcaos2_lcaos2,
		range_jys, index_jys,
		range_lcaos, index_lcaos );

std::cout<<"I"<<std::endl;

}


void Exx_Abfs::test_abfs1() const
{

std::cout<<"A"<<std::endl;

std::cout<<"A1"<<std::endl;

//for(const auto &f : this->files_abfs)
//	std::cout<<f<<std::endl;

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		&&lcaos = Construct_Orbs::change_orbs( GlobalC::ORB, this->kmesh_times );
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		&&abfs_same_atom = Construct_Orbs::abfs_same_atom( lcaos, this->kmesh_times );
	// Peize Lin test
//	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
//		&&abfs_same_atom = IO::Abfs_Same_Atom::construct( GlobalC::ORB );
//	for( size_t T=0; T!=abfs_same_atom.size(); ++T )
//		if( abfs_same_atom[T].size()>4 )
//			abfs_same_atom[T].resize(4);

const ModuleBase::Element_Basis_Index::Range
	&&range_asa = Abfs_Index::construct_range( abfs_same_atom );
std::cout<<range_asa<<std::endl;
std::cout<<__FILE__<<__LINE__<<std::endl;
throw exception();

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
//		&&abfs = IO::construct_abfs( GlobalC::ORB, this->files_abfs, this->kmesh_times );						// Peize Lin test
//		&&abfs = IO::construct_abfs( abfs_same_atom, GlobalC::ORB, this->files_abfs, this->kmesh_times );						// Peize Lin test
		&abfs = abfs_same_atom;

	for( size_t T=0; T!=abfs.size(); ++T )
		this->Lmax = std::max( this->Lmax, static_cast<int>(abfs[T].size())-1 );

std::cout<<Exx_Abfs::Lmax<<std::endl;

	// Peize Lin test
	auto ofs_N_orbital = [&]()
	{
		std::ofstream ofsN("N_orbital.dat");
		for( size_t T=0; T!=abfs_same_atom.size(); ++T )
			for( size_t L=0; L!=abfs_same_atom[T].size(); ++L )
				for( size_t N=0; N!=abfs_same_atom[T][L].size(); ++N )
				{
					ofsN<<T<<"\t"<<L<<"\t"<<N<<std::endl;
					for( size_t ir=0; ir!=abfs_same_atom[T][L][N].getNr(); ++ir )
						ofsN<<abfs_same_atom[T][L][N].getPsi(ir)<<std::endl;
					ofsN<<std::endl;
				}
		ofsN.close();
	};

std::cout<<"B"<<std::endl;

std::cout<<"C"<<std::endl;

std::cout<<"D"<<std::endl;

	Matrix_Lcaoslcaos_Lcaoslcaos mllll;
	mllll.init(1);
	mllll.init_radial( GlobalC::ORB, GlobalC::ORB );
	mllll.init_radial_table();

std::cout<<"D1"<<std::endl;

	Matrix_Orbs21 m_abfslcaos_lcaos;					// "asa" means "abfs_same_atom"
std::cout<<"D11"<<std::endl;
	m_abfslcaos_lcaos.init(1,this->kmesh_times);
std::cout<<"D12"<<std::endl;
	m_abfslcaos_lcaos.init_radial( abfs, lcaos, lcaos );
std::cout<<"D13"<<std::endl;
	m_abfslcaos_lcaos.init_radial_table();

std::cout<<"D2"<<std::endl;

	Matrix_Orbs11 m_abfs_abfs;
std::cout<<"D21"<<std::endl;
	m_abfs_abfs.init(2,this->kmesh_times);
std::cout<<"D22"<<std::endl;
	m_abfs_abfs.init_radial( abfs, abfs );
std::cout<<"D23"<<std::endl;
	m_abfs_abfs.init_radial_table();

std::cout<<"D3"<<std::endl;

std::cout<<"E"<<std::endl;

	const ModuleBase::Element_Basis_Index::Range
		&&range_lcaos = Abfs_Index::construct_range( GlobalC::ORB );
	const ModuleBase::Element_Basis_Index::IndexLNM
		&&index_lcaos = ModuleBase::Element_Basis_Index::construct_index( range_lcaos );

std::cout<<"F"<<std::endl;

	const ModuleBase::Element_Basis_Index::Range
		&&range_abfs = Abfs_Index::construct_range( abfs );
	const ModuleBase::Element_Basis_Index::IndexLNM
		&&index_abfs = ModuleBase::Element_Basis_Index::construct_index( range_abfs );

//std::cout<<range_asa<<std::endl;
//std::cout<<index_asa<<std::endl;

std::cout<<"G"<<std::endl;

std::cout<<"G2"<<std::endl;

	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
		&&ms_lcaos2_lcaos2 = mllll.cal_overlap_matrix( index_lcaos, index_lcaos );

ofs_ms("ms_lcaos2_lcaos2",ms_lcaos2_lcaos2);

std::cout<<"G3"<<std::endl;

	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>>
		&&ms_lcaos2_abfs = m_abfslcaos_lcaos.cal_overlap_matrix( index_abfs, index_lcaos, index_lcaos );

ofs_ms("ms_lcaos2_abfs",ms_lcaos2_abfs);

std::cout<<"G4"<<std::endl;

	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
		&&ms_abfs_abfs = m_abfs_abfs.cal_overlap_matrix( index_abfs, index_abfs );

ofs_ms("ms_abfs_abfs",ms_abfs_abfs);

std::cout<<"G5"<<std::endl;

	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>>
		&&ms_abfs_abfs_I = cal_I( ms_abfs_abfs, index_abfs );

ofs_ms("ms_abfs_abfs_I",ms_abfs_abfs_I);

std::cout<<"G6"<<std::endl;

	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
		&&ms_lcaos2_lcaos2_proj_abfs = cal_lcaos2_lcaos2_proj_asa( ms_lcaos2_abfs, ms_abfs_abfs_I, range_lcaos, index_lcaos );

ofs_ms("ms_lcaos2_lcaos2_proj_abfs",ms_lcaos2_lcaos2_proj_abfs);

std::cout<<"G7"<<std::endl;

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
						ModuleBase::matrix m = ms_lcaos2_lcaos2[TA][IA][TB][IB] - m4.second;
						for( int ir=0; ir!=m.nr; ++ir )
							for( size_t ic=0; ic!=m.nc; ++ic )
								m(ir,ic) /= ms_lcaos2_lcaos2[TA][IA][TB][IB](ir,ic);
						std::cout<<average(m)<<"\t"<<m.max()<<std::endl;
					}
				}
			}
		}
	};

test_ratio();

	std::function< void( ModuleBase::matrix &, const ModuleBase::matrix & ) >
		minus_matrix = []( ModuleBase::matrix &mA, const ModuleBase::matrix &mB ){ mA-=mB; };
	// ms_lcaos2_lcaos2 -= ms_lcaos2_lcaos2_proj_asa
	ModuleBase::GlobalFunc::FUNC_EACH_2( ms_lcaos2_lcaos2, ms_lcaos2_lcaos2_proj_abfs, minus_matrix );

ofs_ms("ms_lcaos2_lcaos2_new",ms_lcaos2_lcaos2);

std::cout<<"H"<<std::endl;

for(const auto m1 : ms_lcaos2_lcaos2)
	for(const auto m2 : m1.second)
		for(const auto m3 : m2.second)
			for(const auto m4 : m3.second)
				std::cout<<average(m4.second)<<"\t"<<m4.second.max()<<std::endl;

std::cout<<"I"<<std::endl;

}

void Exx_Abfs::cal_exx() const
{
	// ȫ����ֻһ��
std::cout<<"A"<<std::endl;

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		&&lcaos = Construct_Orbs::change_orbs( GlobalC::ORB, this->kmesh_times );
//		&&lcaos = Construct_Orbs::change_orbs( GlobalC::ORB, 1 );

std::cout<<"A1"<<std::endl;

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		&&abfs_same_atom = Construct_Orbs::abfs_same_atom( lcaos, this->kmesh_times );
std::cout<<"A2"<<std::endl;
//	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
//		&&abfs_origin = IO::construct_abfs( abfs_same_atom, GlobalC::ORB, this->files_abfs, this->kmesh_times );		// Peize Lin test
//		&&abfs = IO::construct_abfs( GlobalC::ORB, this->files_abfs, this->kmesh_times );						// Peize Lin test
std::cout<<"A3"<<std::endl;
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
//		&&abfs = Construct_Orbs::orth_orbs( abfs_origin );		// Peize Lin test
		&&abfs = Construct_Orbs::orth_orbs( abfs_same_atom );		// Peize Lin test
std::cout<<"A4"<<std::endl;

	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_ccp;
	Conv_Coulomb_Pot::cal_orbs_ccp( abfs, abfs_ccp, this->rmesh_times, 1 );

	for( size_t T=0; T!=abfs.size(); ++T )
		this->Lmax = std::max( this->Lmax, static_cast<int>(abfs[T].size())-1 );

std::cout<<"B"<<std::endl;

	const ModuleBase::Element_Basis_Index::Range
		&&range_lcaos = Abfs_Index::construct_range( lcaos );
	const ModuleBase::Element_Basis_Index::IndexLNM
		&&index_lcaos = ModuleBase::Element_Basis_Index::construct_index( range_lcaos );

std::cout<<"C"<<std::endl;

	const ModuleBase::Element_Basis_Index::Range
		&&range_abfs = Abfs_Index::construct_range( abfs );
	const ModuleBase::Element_Basis_Index::IndexLNM
		&&index_abfs = ModuleBase::Element_Basis_Index::construct_index( range_abfs );

std::cout<<range_abfs<<std::endl;

std::cout<<"D"<<std::endl;

	// Peize Lin test
	auto test_matrix_lcaos_lcaos = [&]()
	{
		Matrix_Orbs11 m_lcaos_lcaos;
		m_lcaos_lcaos.init(1);
		m_lcaos_lcaos.init_radial( GlobalC::ORB, GlobalC::ORB );
		m_lcaos_lcaos.init_radial_table();
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
			&&matrix_V = m_lcaos_lcaos.cal_overlap_matrix(index_lcaos,index_lcaos);

		std::ofstream ofs("S_matrix.dat");
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
						ofs<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<std::endl;
						m4.second.print(ofs, 1E-10)<<std::endl;
					}
				}
			}
		}
		ofs.close();
		ModuleBase::WARNING_QUIT( ModuleBase::GlobalFunc::TO_STRING(__FILE__), ModuleBase::GlobalFunc::TO_STRING(__LINE__) );
	};

	// Peize Lin test
	auto test_matrix_lcaos_lcaos2 = [&]()
	{
		Matrix_Orbs11 m_lcaos_lcaos;
		m_lcaos_lcaos.init(1,this->kmesh_times);
		m_lcaos_lcaos.init_radial( lcaos, lcaos );
		m_lcaos_lcaos.init_radial_table();
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
			&&matrix_V = m_lcaos_lcaos.cal_overlap_matrix(index_lcaos,index_lcaos);

		std::ofstream ofs("S_matrix.dat");
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
						ofs<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<std::endl;
						m4.second.print(ofs, 1E-10)<<std::endl;
					}
				}
			}
		}
		ofs.close();
		ModuleBase::WARNING_QUIT( ModuleBase::GlobalFunc::TO_STRING(__FILE__), "line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__) );
	};

	// Peize Lin test
/*
	auto test_matrix_standard_ll2_ll2 = [&]()
	{
		Matrix_Lcaoslcaos_Lcaoslcaos2<Center2_Orb::Orb22_Ccp> m_ll2_ll2;
std::cout<<"D0"<<std::endl;
		m_ll2_ll2.init( 1, this->kmesh_times, this->rmesh_times );
std::cout<<"D1"<<std::endl;
		m_ll2_ll2.init_radial( GlobalC::ORB, GlobalC::ORB );
std::cout<<"D2"<<std::endl;
		m_ll2_ll2.init_radial_table();
std::cout<<"D3"<<std::endl;
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
			&&ms_ll2_ll2 = m_ll2_ll2.cal_overlap_matrix(index_lcaos,index_lcaos);
std::cout<<"D4"<<std::endl;
		ofs_ms("ms_ll2_ll2",ms_ll2_ll2);
	};
	test_matrix_standard_ll2_ll2();
*/
std::cout<<"DD"<<std::endl;

	Matrix_Orbs11 m_abfs_abfs;
std::cout<<"D1"<<std::endl;
	m_abfs_abfs.init( 2, this->kmesh_times, this->rmesh_times );
std::cout<<"D2"<<std::endl;
	m_abfs_abfs.init_radial( abfs, abfs_ccp );
std::cout<<"D3"<<std::endl;
	m_abfs_abfs.init_radial_table();

std::cout<<"E"<<std::endl;

	Matrix_Orbs21 m_abfslcaos_lcaos;
	m_abfslcaos_lcaos.init( 1, this->kmesh_times, this->rmesh_times );
	m_abfslcaos_lcaos.init_radial( abfs_ccp, lcaos, lcaos );
	m_abfslcaos_lcaos.init_radial_table();

std::cout<<"F"<<std::endl;

	// ÿһ�����Ӳ�
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>>
		&&ms_abfs_abfs = m_abfs_abfs.cal_overlap_matrix( index_abfs, index_abfs );
ofs_ms("ms_abfs_abfs",ms_abfs_abfs);

	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>>
		&&ms_lcaos2_abfs = m_abfslcaos_lcaos.cal_overlap_matrix( index_abfs, index_lcaos, index_lcaos );
ofs_ms("ms_lcaos2_abfs",ms_lcaos2_abfs);

std::cout<<"G"<<std::endl;

	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>>
		&&ms_abfs_abfs_I = cal_I( ms_abfs_abfs, index_abfs );
ofs_ms("ms_abfs_abfs_I",ms_abfs_abfs_I);

std::cout<<"G1"<<std::endl;
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>>
		&&ms_C = cal_C( ms_lcaos2_abfs, ms_abfs_abfs_I );
ofs_ms("ms_C",ms_C);

std::cout<<"H"<<std::endl;

	// ÿһ�����Ӳ�
	timeval t_begin;
	gettimeofday( &t_begin, NULL);

	cal_CVC( ms_C, ms_abfs_abfs );

	timeval t_end;
	gettimeofday( &t_end, NULL);
	std::cout<<"time:\t"<<(double)(t_end.tv_sec-t_begin.tv_sec) + (double)(t_end.tv_usec-t_begin.tv_usec)/1000000.0<<std::endl;

std::cout<<"I"<<std::endl;

}

void Exx_Abfs::test_abfs2() const{}

// <a|a> .I
// &
// <a|a> <a|b>
// <b|a> <b|b> .I
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> Exx_Abfs::cal_I(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &ms,
	const ModuleBase::Element_Basis_Index::IndexLNM &index )
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> ms_I;

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
//					std::cout<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<std::endl;
//					const ModuleBase::matrix matrix_origin( ms_tmp.A );
//					std::cout<<matrix_origin<<std::endl;

					ms_tmp.cal_inverse( Exx_Abfs::Inverse_Matrix_Double::Method::dpotrf );

					// Peize Lin test
//					std::cout<<ms_tmp.A<<std::endl;
//					std::cout<<ms_tmp.A * matrix_origin <<std::endl;
//					std::cout<<matrix_origin * ms_tmp.A<<std::endl;

					if( TA==TB && IA==IB )
					{
						const size_t size_A = index[TA].count_size;

						ms_I[TA][IA][TB][IB].resize( 1, std::vector<ModuleBase::matrix>(1) );
						ms_I[TA][IA][TB][IB][0][0].create( size_A, size_A );

						ms_tmp.output(
							ms_I[TA][IA][TB][IB][0][0]);
					}
					else
					{
						const size_t size_A = index[TA].count_size;
						const size_t size_B = index[TB].count_size;

						ms_I[TA][IA][TB][IB].resize( 2, std::vector<ModuleBase::matrix>(2) );
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
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> Exx_Abfs::cal_C(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms_lcaos2_abfs,
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> &ms_abfs_abfs_I )
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> ms_C;

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
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms_C,
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &ms_abfs_abfs ) const
{
	std::ofstream ofs("ms_CVC");

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
									ModuleBase::matrix mm;		// Peize Lin test

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
									ofs<<IA<<"\t"<<IB<<"\t"<<IC<<"\t"<<ID<<std::endl;
									mm.print(ofs, 1E-10)<<std::endl;
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
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> Exx_Abfs::cal_lcaos2_lcaos2_proj_asa(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms_lcaos2_asa,
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> &ms_asa_asa_I,
	const ModuleBase::Element_Basis_Index::Range &range,
	const ModuleBase::Element_Basis_Index::IndexLNM &index)
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> ms_lcaos2_lcaos2_proj_asa;
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

					std::vector<ModuleBase::matrix> mql(2);
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
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> Exx_Abfs::cal_lcaos2_jys_proj_asa(
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms_lcaos2_asa,
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> &ms_asa_asa_I,
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &ms_asa_jys)
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> ms_lcaos2_jys_proj_asa;

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
						std::vector<ModuleBase::matrix> ms_tmp(2);
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
	std::vector<Vector3_Exx> R_supercell;
	for( size_t x=0; x!=GlobalC::kv.nmp[0]; ++x)
		for( size_t y=0; y!=GlobalC::kv.nmp[1]; ++y )
			for( size_t z=0; z!=GlobalC::kv.nmp[2]; ++z )
				R_supercell.push_back(Vector3_Exx(x,y,z));
}
*/
/*
void density_matrix()
{
	std::vector<ModuleBase::matrix> DM_k(GlobalC::kv.nks);
	for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
	{
		for( size_t ib=0; ib!=GlobalV::NBANDS; ++ib )
		{
			for( size_t iw1=0; iw1!=GlobalV::NLOCAL; ++iw1 )
			{
				for( size_t iw2=0; iw2!=GlobalV::NLOCAL; ++iw2 )
				{
					DM_k[ik](iw1,iw2) += GlobalC::wf.wg(ik,ib) * conj(GlobalC::LOWF.wfc_k_grid[ik][ib][iw1]) * GlobalC::LOWF.wfc_k_grid[ik][ib][iw2];
				}
			}
		}
	}

	std::vector<size_t,std::vector<ModuleBase::matrix>> DM_R( GlobalV::NSPIN, std::vector<ModuleBase::matrix>(R_supercell.size()) );
	for( size_t is=0; is!=GlobalV::NSPIN; ++is )
	{
		const size_t k_start = (GlobalV::NSPIN==1) ? 0 : ((is==0) ? 0 : (GlobalC::kv.nks/2));
		const size_t k_end = (GlobalV::NSPIN==1) ? GlobalC::kv.nks : ((is==0) ? (GlobalC::kv.nks/2) : GlobalC::kv.nks);
		for( size_it iR=0; iR!=R_supercell.size(); ++iR )
		{
			for( size_t ik=k_start; ik!=k_end; ++ik )
			{
				DM_R[is][iR] += exp(ModuleBase::TWO_PI*IMAG_UNIT*dot(GlobalC::kv.kvec_d[ik],R_supercell[iR])) * DM[ik];
			}
		}
	}
}
*/
