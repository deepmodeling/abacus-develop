#ifdef __MPI   
#include "exx_lcao.h"

#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../src_parallel/parallel_reduce.h"

#include "conv_coulomb_pot.h"
#include "conv_coulomb_pot-inl.h"
#include "conv_coulomb_pot_k.h"
#include "conv_coulomb_pot_k-template.h"

#include "abfs-template.h"
#include "exx_abfs.h"
#include "exx_abfs-construct_orbs.h"
#include "exx_abfs-abfs_index.h"
#include "exx_abfs-screen-schwarz.h"
#include "exx_abfs-io.h"
//#include "exx_abfs-io-template.h"
#include "exx_abfs-util.h"
#include "exx_abfs-parallel.h"
#include "exx_abfs-parallel-distribute-htime.h"
#include "exx_abfs-parallel-distribute-kmeans.h"
#include "exx_abfs-parallel-distribute-order.h"

#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

// Peize Lin test
#include<stdexcept>	
#include<sys/time.h>
#include "../src_external/src_test/src_global/element_basis_index-test.h"
#include "../src_external/src_test/src_global/matrix-test.h"
#include "../src_external/src_test/print_tmp.h"
#include "../src_external/src_test/src_global/sph_bessel-unittest.h"
#include "../src_external/src_test/src_ri/exx_lcao-test.h"
#include "../src_external/src_test/src_ri/abfs-test.h"

/*
// m_new( i2*n1+i1, i3 ) = m( i1*n2+i2, i3 )
static ModuleBase::matrix transform (
	const ModuleBase::matrix & m,
	const size_t n1, const size_t n2, const size_t n3 )
{
	assert( n1*n2*n3 == m.nr*m.nc );
	const auto length = sizeof(double)*n3;
	ModuleBase::matrix m_new( n1*n2, n3, false );
	for( size_t i1=0; i1!=n1; ++i1 )
	{
		for( size_t i2=0; i2!=n2; ++i2 )
			memcpy( m_new.c+(i2*n1+i1)*n3, m.c+(i1*n2+i2)*n3, length );
	}
	return m_new;
}
*/


// Peize Lin test
Exx_Lcao::Exx_Lcao( const Exx_Global::Exx_Info &info_global )
	:kmesh_times(4),
	 info(info_global)
{
	auto test_gemm = []()
	{
		auto init_matrix = [](const int nr, const int nc, const double add) -> ModuleBase::matrix
		{
			ModuleBase::matrix m(nr,nc);
			for(int i=0; i<m.nr*m.nc; ++i)
				m.c[i]=i+add;
			return m;
		};
		auto transpose = [](const ModuleBase::matrix &m) -> ModuleBase::matrix
		{
			ModuleBase::matrix mT(m.nc,m.nr);
			for(int ir=0; ir!=m.nr; ++ir)
				for(int ic=0; ic!=m.nc; ++ic)
					mT(ic,ir) = m(ir,ic);
			return mT;
		};
		auto print_matrix = [](const ModuleBase::matrix &m1,const ModuleBase::matrix &m2,const ModuleBase::matrix &m3)
		{
			m1.print(std::cout, 1E-10)<<std::endl;
			m2.print(std::cout, 1E-10)<<std::endl;
			m3.print(std::cout, 1E-10)<<std::endl;
			std::cout<<"============================="<<std::endl<<std::endl;
		};
		{
			const ModuleBase::matrix m1=init_matrix(1,3,10), m2=init_matrix(3,2,0);
			ModuleBase::matrix m3=m1*m2;
			print_matrix(m1,m2,m3);
		}
		{
			const ModuleBase::matrix m1=init_matrix(1,3,10), m2=init_matrix(3,2,0);
			ModuleBase::matrix m3(m1.nr,m2.nc);
			BlasConnector::gemm(
				'N', 'N', 
				m3.nr, m3.nc, m1.nc,
				1, m1.c, m1.nc, m2.c, m2.nc,
				0, m3.c, m3.nc);
			print_matrix(m1,m2,m3);
		}
		{
			const ModuleBase::matrix m1=transpose(init_matrix(1,3,10)), m2=init_matrix(3,2,0);
			ModuleBase::matrix m3(m1.nc,m2.nc);
			BlasConnector::gemm(
				'T', 'N', 
				m3.nr, m3.nc, m1.nr,
				1, m1.c, m1.nc, m2.c, m2.nc,
				0, m3.c, m3.nc);
			print_matrix(m1,m2,m3);
		}
		{
			const ModuleBase::matrix m1=transpose(init_matrix(1,3,10)), m2=init_matrix(3,2,0);
			ModuleBase::matrix m3(m1.nc,m2.nc);
			BlasConnector::gemm(
				'N', 'T', 
				m3.nr, m3.nc, m1.nr,
				1, m1.c, m1.nc, m2.c, m2.nc,
				0, m3.c, m3.nc);
			print_matrix(m1,m2,m3);
		}
		{
			const ModuleBase::matrix m1=init_matrix(1,3,10), m2=transpose(init_matrix(3,2,0));
			ModuleBase::matrix m3(m1.nr,m2.nr);
			BlasConnector::gemm(
				'N','T',
				m3.nr, m3.nc, m1.nc,
				1, m1.c, m1.nc, m2.c, m2.nc,
				0, m3.c, m3.nc);
			print_matrix(m1,m2,m3);
		}
		{
			const ModuleBase::matrix m1=init_matrix(1,3,10), m2=transpose(init_matrix(3,2,0));
			ModuleBase::matrix m3(m1.nr,m2.nr);
			BlasConnector::gemm(
				'T','N',
				m3.nr, m3.nc, m1.nc,
				1, m1.c, m1.nc, m2.c, m2.nc,
				0, m3.c, m3.nc);
			print_matrix(m1,m2,m3);
		}
	};

	auto test_gemm_2 = []()
	{
		constexpr int S=10000;
		constexpr int N=1000;
		ModuleBase::matrix m1(N,N), m2(N,N), m3(N,N);	
		timeval time;		
		gettimeofday(&time, NULL);
		
		for(int s=0; s<S; ++s)
			BlasConnector::gemm('N', 'N', N, N, N, 1, m1.c, N, m2.c, N, 0, m3.c, N);
		std::cout<<"NN\t"<<cut_time(time)<<std::endl;
		
		for(int s=0; s<S; ++s)
			BlasConnector::gemm('N', 'T', N, N, N, 1, m1.c, N, m2.c, N, 0, m3.c, N);
		std::cout<<"NT\t"<<cut_time(time)<<std::endl;
		
		for(int s=0; s<S; ++s)
			BlasConnector::gemm('T', 'N', N, N, N, 1, m1.c, N, m2.c, N, 0, m3.c, N);
		std::cout<<"TN\t"<<cut_time(time)<<std::endl;
		
		for(int s=0; s<S; ++s)
			BlasConnector::gemm('T', 'T', N, N, N, 1, m1.c, N, m2.c, N, 0, m3.c, N);
		std::cout<<"TT\t"<<cut_time(time)<<std::endl;
	};

	auto test_gemm_3 = []()
	{
		double s=0;
		const int M=1;
		for(int N=1; N<=4096; N*=2)
		{
			ModuleBase::matrix a(N,N), b(N,N);
			timeval t;	gettimeofday(&t,NULL);
			for(int m=0; m<M; ++m)
			{
				ModuleBase::matrix c = a * b;
				s+=c(0,0);
			}
			std::cout<<N<<"\t"<<cal_time(t)<<"\t"<<s<<std::endl;
		}
	};
}

Exx_Lcao::Exx_Info::Exx_Info( const Exx_Global::Exx_Info &info_global )
	:hybrid_type(info_global.hybrid_type),
	 hse_omega(info_global.hse_omega){} 

void Exx_Lcao::init()
{
	auto mkdir_test_dir = [&]()
	{
		auto mkdir_one = [](const std::string &dir)
		{
			const std::string command0 =  "test -d " + dir + " || mkdir " + dir;
			if(GlobalV::MY_RANK==0)
			{
				system( command0.c_str() );
			}
		};
		test_dir = {"test_exx/process/","test_exx/thread/","test_exx/matrix/"};
		mkdir_one("test_exx");
		mkdir_one(test_dir.process);
		mkdir_one(test_dir.thread);
		mkdir_one(test_dir.matrix);
	};	
	
	auto test_rk = []()
	{
		auto pr_v = []( const std::string & file, const std::vector<double> & v )
		{
			std::ofstream ofs(file);
			for( size_t i=0; i!=v.size(); ++i )
			{
				ofs<<v[i]<<std::endl;
			}
			ofs.close();
		};
		auto pr_orb = [&pr_v]( const std::string & file, const Numerical_Orbital_Lm & orb )
		{
			pr_v( file+"-psi", orb.get_psi() );
			pr_v( file+"-psi_f", orb.get_psif() );
			pr_v( file+"-psi_k", orb.get_psi_k() );
			pr_v( file+"-psi_k2", orb.get_psi_k2() );
		};
		auto pr_orb_all = [&pr_orb]( const std::string & file, const Numerical_Orbital_Lm & orb )
		{
			Numerical_Orbital_Lm orb_psi_T, orb_psi_F, orb_psif_T, orb_psik_T, orb_psik2_T;
			orb_psi_T.set_orbital_info
			(
				orb.getLabel(),
				orb.getType(),
				orb.getL(),
				orb.getChi(),
				orb.getNr(),
				orb.getRab(),
				orb.getRadial(),
				Numerical_Orbital_Lm::Psi_Type::Psi,
				orb.getPsi(),
				orb.getNk(),
				orb.getDk(),
				orb.getDruniform(),
				false,
				true, GlobalV::CAL_FORCE
			);
			orb_psi_F.set_orbital_info
			(
				orb.getLabel(),
				orb.getType(),
				orb.getL(),
				orb.getChi(),
				orb.getNr(),
				orb.getRab(),
				orb.getRadial(),
				Numerical_Orbital_Lm::Psi_Type::Psi,
				orb.getPsi(),
				orb.getNk(),
				orb.getDk(),
				orb.getDruniform(),
				false,
				false, GlobalV::CAL_FORCE
			);
			orb_psif_T.set_orbital_info
			(
				orb.getLabel(),
				orb.getType(),
				orb.getL(),
				orb.getChi(),
				orb.getNr(),
				orb.getRab(),
				orb.getRadial(),
				Numerical_Orbital_Lm::Psi_Type::Psif,
				orb.getPsif(),
				orb.getNk(),
				orb.getDk(),
				orb.getDruniform(),
				false,
				true, GlobalV::CAL_FORCE
			);
			orb_psik_T.set_orbital_info
			(
				orb.getLabel(),
				orb.getType(),
				orb.getL(),
				orb.getChi(),
				orb.getNr(),
				orb.getRab(),
				orb.getRadial(),
				Numerical_Orbital_Lm::Psi_Type::Psik,
				orb.getPsi_k(),
				orb.getNk(),
				orb.getDk(),
				orb.getDruniform(),
				false,
				true, GlobalV::CAL_FORCE
			);
			orb_psik2_T.set_orbital_info
			(
				orb.getLabel(),
				orb.getType(),
				orb.getL(),
				orb.getChi(),
				orb.getNr(),
				orb.getRab(),
				orb.getRadial(),
				Numerical_Orbital_Lm::Psi_Type::Psik2,
				orb.getPsi_k2(),
				orb.getNk(),
				orb.getDk(),
				orb.getDruniform(),
				false,
				true, GlobalV::CAL_FORCE
			);
			pr_orb( file+"-orb",orb );
			pr_orb( file+"-orb_psi_T",orb_psi_T );
			pr_orb( file+"-orb_psi_F",orb_psi_F );
			pr_orb( file+"-orb_psif_T",orb_psif_T );
			pr_orb( file+"-orb_psik_T",orb_psik_T );
			pr_orb( file+"-orb_psik2_T",orb_psik2_T );
		};
		auto pr_orb_all_kmesh = [&pr_orb]( const std::string & file, const Numerical_Orbital_Lm & orb, const int kmesh_times )
		{
			const int Nk =orb.getNk() * kmesh_times | 1;
			Numerical_Orbital_Lm orb_psi_T, orb_kmesh;
			orb_psi_T.set_orbital_info
			(
				orb.getLabel(),
				orb.getType(),
				orb.getL(),
				orb.getChi(),
				orb.getNr(),
				orb.getRab(),
				orb.getRadial(),
				Numerical_Orbital_Lm::Psi_Type::Psi,
				orb.getPsi(),
				Nk,
				orb.getDk(),
				orb.getDruniform(),
				false,
				true, GlobalV::CAL_FORCE
			);
			orb_kmesh.set_orbital_info
			(
				orb_psi_T.getLabel(),
				orb_psi_T.getType(),
				orb_psi_T.getL(),
				orb_psi_T.getChi(),
				orb_psi_T.getNr(),
				orb_psi_T.getRab(),
				orb_psi_T.getRadial(),
				Numerical_Orbital_Lm::Psi_Type::Psif,
				orb_psi_T.getPsif(),
				Nk,
				orb_psi_T.getDk(),
				orb_psi_T.getDruniform(),
				false,
				true, GlobalV::CAL_FORCE
			);
			pr_orb( file+"-orb",orb );
			pr_orb( file+"-orb_kmesh",orb_kmesh );
		};
		if(false)
		{
			for( int T=0; T!=GlobalC::ORB.get_ntype(); ++T )
			{
				for( int L=0; L<=GlobalC::ORB.Phi[T].getLmax(); ++L )
				{
					for( int N=0; N!=GlobalC::ORB.Phi[T].getNchi(L); ++N )
					{
						pr_orb_all( "orb_"+ModuleBase::GlobalFunc::TO_STRING(T)+"_"+ModuleBase::GlobalFunc::TO_STRING(L)+"_"+ModuleBase::GlobalFunc::TO_STRING(N), GlobalC::ORB.Phi[T].PhiLN(L,N) );
					}
				}
			}
		}
		else
		{
			for( int T=0; T!=GlobalC::ORB.get_ntype(); ++T )
			{
				for( int L=0; L<=GlobalC::ORB.Phi[T].getLmax(); ++L )
				{
					for( int N=0; N!=GlobalC::ORB.Phi[T].getNchi(L); ++N )
					{
						pr_orb_all_kmesh( "orb_"+ModuleBase::GlobalFunc::TO_STRING(T)+"_"+ModuleBase::GlobalFunc::TO_STRING(L)+"_"+ModuleBase::GlobalFunc::TO_STRING(N), GlobalC::ORB.Phi[T].PhiLN(L,N), 5 );
					}
				}
			}
		}
	};

	auto test_exp = [&]()
	{
		std::cout<<"GlobalC::kv.kvec_d:"<<std::endl;
		for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
		{
			std::cout<<GlobalC::kv.kvec_d[ik]<<std::endl;
		}
		std::cout<<"GlobalC::kv.kvec_c:"<<std::endl;
		for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
		{
			std::cout<<GlobalC::kv.kvec_c[ik]<<std::endl;
		}

		const ModuleBase::Vector3<int>BvK_period( GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2] );
		std::vector<ModuleBase::Vector3<double>> boxes;
		for( int x=0; x!=BvK_period.x; ++x )
		{
			for( int y=0; y!=BvK_period.y; ++y )
			{
				for( int z=0; z!=BvK_period.z; ++z )
				{
					boxes.push_back({x,y,z});
				}
			}
		}

		std::cout<<"boxes:"<<std::endl;
		for( size_t i=0; i!=boxes.size(); ++i )
		{
			std::cout<<boxes[i]<<std::endl;
		}
		std::cout<<"box * GlobalC::ucell.latvec:"<<std::endl;
		for( size_t i=0; i!=boxes.size(); ++i )
		{
			std::cout<<boxes[i]*GlobalC::ucell.latvec<<std::endl;
		}

		std::cout<<"k R"<<std::endl;
		for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
		{
			for( size_t i=0; i!=boxes.size(); ++i )
			{
				std::cout<<GlobalC::kv.kvec_c[ik] * (boxes[i]*GlobalC::ucell.latvec)<<"\t";
			}
			std::cout<<std::endl;
		}
		std::cout<<"exp( - 2 pi i k R )"<<std::endl;
		for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
		{
			for( size_t i=0; i!=boxes.size(); ++i )
			{
				std::cout<<exp( -ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT* (GlobalC::kv.kvec_c[ik]* (boxes[i]*GlobalC::ucell.latvec)) )<<"\t";
			}
			std::cout<<std::endl;
		}
		std::cout<<"k R"<<std::endl;
		for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
		{
			for( size_t i=0; i!=boxes.size(); ++i )
			{
				std::cout<<GlobalC::kv.kvec_d[ik] * static_cast<ModuleBase::Vector3<double>>(boxes[i])<<"\t";
			}
			std::cout<<std::endl;
		}
		std::cout<<"exp( - 2 pi i k R )"<<std::endl;
		for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
		{
			for( size_t i=0; i!=boxes.size(); ++i )
			{
				std::cout<<exp( -ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT* (GlobalC::kv.kvec_d[ik]* static_cast<ModuleBase::Vector3<double>>(boxes[i])) )<<"\t";
			}
			std::cout<<std::endl;
		}

		std::cout<<"Rcut:"<<std::endl;
		for( size_t T=0; T!=GlobalC::ORB.get_ntype(); ++T )
		{
			std::cout<<GlobalC::ORB.Phi[T].getRcut()<<std::endl;
		}
		std::cout<<"tau:"<<std::endl;
		for( size_t iat=0; iat!=GlobalC::ucell.nat; ++iat )
		{
			std::cout<<GlobalC::ucell.atoms[ GlobalC::ucell.iat2it[iat] ].tau[ GlobalC::ucell.iat2ia[iat] ]<<std::endl;
		}
		std::cout<<"taud:"<<std::endl;
		for( size_t iat=0; iat!=GlobalC::ucell.nat; ++iat )
		{
			std::cout<<GlobalC::ucell.atoms[ GlobalC::ucell.iat2it[iat] ].taud[ GlobalC::ucell.iat2ia[iat] ]<<std::endl;
		}
		std::cout<<"taud * latvec:"<<std::endl;
		for( size_t iat=0; iat!=GlobalC::ucell.nat; ++iat )
		{
			std::cout<<GlobalC::ucell.atoms[ GlobalC::ucell.iat2it[iat] ].taud[ GlobalC::ucell.iat2ia[iat] ] * GlobalC::ucell.latvec<<std::endl;
		}
		std::cout<<"GlobalC::ucell.latvec:"<<std::endl;
		GlobalC::ucell.latvec.print();
	};

	auto test_matrix = []()
	{
		ModuleBase::matrix m(3,3);
		for( size_t ir=0; ir!=m.nr; ++ir )
		{
			for( size_t ic=0; ic!=m.nc; ++ic )
			{
				m(ir,ic) = ir*10+ic;
			}
		}
		ModuleBase::ComplexMatrix cm = ModuleBase::ComplexMatrix(m) * exp( -ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT* 1.0/3.0 );
		m.print(std::cout, 1E-10)<<std::endl;
		cm.print(std::cout, 1E-10, 1E-10)<<std::endl;
	};

	auto test_nrm2 = []()
	{
		std::vector<double> x = {1,2,3};
		std::cout<<BlasConnector::nrm2( x.size(), ModuleBase::GlobalFunc::VECTOR_TO_PTR(x), 1 )<<std::endl;
		std::vector<std::complex<double>> y = { {1.1,2.2}, {3.3,-4.4}, {-5.5,-6.6} };
		std::cout<<BlasConnector::nrm2( y.size(), ModuleBase::GlobalFunc::VECTOR_TO_PTR(y), 1 )<<std::endl;
		std::vector<double> z = {1,2,3,4,5,6};
		std::cout<<BlasConnector::nrm2( 3, ModuleBase::GlobalFunc::VECTOR_TO_PTR(z), 2 )<<std::endl;
	};

	ModuleBase::TITLE("Exx_Lcao","init");

mkdir_test_dir();

std::ofstream ofs_mpi(test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
timeval t_start,t_start_all;
gettimeofday( &t_start_all, NULL);

//	DM.flag_mix = info.separate_loop ? false : true;
//	DM.flag_mix = false;		// Peize Lin test

#ifdef __MPI
	if(GlobalC::exx_global.info.separate_loop)
	{
		Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::No;
		Hexx_para.mixing_beta = 0;
	}
	else
	{
		if("plain"==GlobalC::CHR.mixing_mode)
		{
			Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Plain;
		}
		else if("pulay"==GlobalC::CHR.mixing_mode)
		{
			Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Pulay;
		}
		else
		{
			throw std::invalid_argument("exx mixing error. exx_separate_loop==false, mixing_mode!=plain or pulay");
		}
		Hexx_para.mixing_beta = GlobalC::CHR.mixing_beta;
	}
#endif

gettimeofday( &t_start, NULL);
	this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs( GlobalC::ORB, this->kmesh_times );
ofs_mpi<<"TIME@ Exx_Abfs::Construct_Orbs::change_orbs\t"<<time_during(t_start)<<std::endl;

ofs_mpi<<info.files_abfs<<std::endl;
#ifdef __MPI
	Exx_Abfs::Util::bcast( info.files_abfs, 0, MPI_COMM_WORLD );
#endif
ofs_mpi<<info.files_abfs<<std::endl;

gettimeofday( &t_start, NULL);
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		abfs_same_atom = Exx_Abfs::Construct_Orbs::abfs_same_atom( lcaos, this->kmesh_times, info.pca_threshold );		// Peize Lin test
	if(info.files_abfs.empty())
	{
		this->abfs = abfs_same_atom;
	}
	else
	{
		this->abfs = Exx_Abfs::IO::construct_abfs( abfs_same_atom, GlobalC::ORB, info.files_abfs, this->kmesh_times );
	}
//	this->abfs = Exx_Abfs::Construct_Orbs::orth_orbs( abfs_origin );		// Peize Lin test
ofs_mpi<<"TIME@ Exx_Abfs::Construct_Orbs::abfs\t"<<time_during(t_start)<<std::endl;

	auto print_psi1 = [](const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs)
	{
		for(size_t N=0; N!=orbs[0][0].size(); ++N)
		{
			for(size_t ir=0; ir!=orbs[0][0][N].getNr(); ++ir)
			{
				std::cout<<orbs[0][0][N].getPsi(ir)<<"\t";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
		for(size_t N=0; N!=orbs[0][0].size(); ++N)
		{
			for(size_t ir=0; ir!=orbs[0][0][N].getNr(); ++ir)
			{
				std::cout<<orbs[0][0][N].getPsi_r(ir)<<"\t";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
		for(size_t N=0; N!=orbs[0][0].size(); ++N)
		{
			for(size_t ik=0; ik!=orbs[0][0][N].getNk(); ++ik)
			{
				std::cout<<orbs[0][0][N].getPsi_k(ik)<<"\t";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	};

	auto print_psi2 = [](
		const std::string & file_name,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orbs)
	{
		std::ofstream ofs(file_name);
		for( size_t T=0; T!=orbs.size(); ++T )
		{
			for( size_t L=0; L!=orbs[T].size(); ++L )
			{
				for( size_t N=0; N!=orbs[T][L].size(); ++N )
				{
//					ofs<<T<<"\t"<<L<<"\t"<<N<<std::endl;
					for(size_t ir=0; ir!=orbs[T][L][N].getNr(); ++ir)
					{
						ofs<<orbs[T][L][N].getPsi(ir)<<"\t";
					}
					ofs<<std::endl;
				}
			}
		}
		ofs.close();
	};

//	Conv_Coulomb_Pot::cal_orbs_ccp( abfs, abfs_ccp, info.ccp_rmesh_times, 1 );
//{
//	std::ofstream ofs("exx_lcao"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
//	ofs<<static_cast<std::underlying_type<Exx_Lcao::Hybrid_Type>::type>(exx_lcao.info.hybrid_type)<<std::endl;
//	ofs.close();
//}

gettimeofday( &t_start, NULL);
	switch(info.hybrid_type)
	{
		case Exx_Global::Hybrid_Type::HF:
		case Exx_Global::Hybrid_Type::PBE0:
			abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp( this->abfs, Conv_Coulomb_Pot_K::Ccp_Type::Ccp, {}, info.ccp_rmesh_times );		break;
		case Exx_Global::Hybrid_Type::HSE:
			abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp( this->abfs, Conv_Coulomb_Pot_K::Ccp_Type::Hse, {{"hse_omega",info.hse_omega}}, info.ccp_rmesh_times );	break;
		default:
			throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));	break;
	}
ofs_mpi<<"TIME@ Conv_Coulomb_Pot_K::cal_orbs_ccp\t"<<time_during(t_start)<<std::endl;

	auto print_psik = [](
		const std::string & file_name,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> & orbs,
		const int power )
	{
		for( size_t i=0; i!=orbs.size(); ++i )
		{
			for( size_t j=0; j!=orbs[i].size(); ++j )
			{
				for( size_t k=0; k!=orbs[i][j].size(); ++k )
				{
					const Numerical_Orbital_Lm & orb = orbs[i][j][k];
					std::ofstream ofs(file_name+"_"+ModuleBase::GlobalFunc::TO_STRING(i)+"_"+ModuleBase::GlobalFunc::TO_STRING(j)+"_"+ModuleBase::GlobalFunc::TO_STRING(k));
					for( size_t ik=0; ik!=orb.getNk(); ++ik )
					{
						ofs<<orb.getPsi_k(ik) / pow(orb.getKpoint(ik),power)<<std::endl;
					}
					ofs.close();
				}
			}
		}
	};

	#if TEST_EXX_LCAO==1
		print_psi2(test_dir.matrix+"r_abfs_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),abfs);
		print_psi2(test_dir.matrix+"r_abfs_ccp_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),abfs_ccp);
	#elif TEST_EXX_LCAO==-1
		#error
	#endif

	for( size_t T=0; T!=abfs.size(); ++T )
	{
		Exx_Abfs::Lmax = std::max( Exx_Abfs::Lmax, static_cast<int>(abfs[T].size())-1 );
	}

ofs_mpi<<"Exx_Abfs::Lmax:\t"<<Exx_Abfs::Lmax<<std::endl;

	const ModuleBase::Element_Basis_Index::Range
		&&range_lcaos = Exx_Abfs::Abfs_Index::construct_range( lcaos );
	index_lcaos = ModuleBase::Element_Basis_Index::construct_index( range_lcaos );

	const ModuleBase::Element_Basis_Index::Range
		&&range_abfs = Exx_Abfs::Abfs_Index::construct_range( abfs );
	index_abfs = ModuleBase::Element_Basis_Index::construct_index( range_abfs );

ofs_mpi<<range_lcaos<<std::endl;
ofs_mpi<<range_abfs<<std::endl;

	auto test_mll = [&]()
	{
		Exx_Abfs::Matrix_Orbs11 mll;
		mll.init(2,1,1);
		mll.init_radial(GlobalC::ORB, GlobalC::ORB);
		mll.init_radial_table();
		std::ofstream ofsS("S.dat");
		mll.cal_overlap_matrix(0,0,GlobalC::ucell.atoms[0].tau[0],GlobalC::ucell.atoms[0].tau[0],index_lcaos,index_lcaos).print(ofsS, 1E-10)<<std::endl<<std::endl;
		mll.cal_overlap_matrix(0,0,GlobalC::ucell.atoms[0].tau[0],GlobalC::ucell.atoms[0].tau[1],index_lcaos,index_lcaos).print(ofsS, 1E-10)<<std::endl<<std::endl;
	};

gettimeofday( &t_start, NULL);
	m_abfs_abfs.init( 2, this->kmesh_times, (1+info.ccp_rmesh_times)/2.0 );
ofs_mpi<<"TIME@ m_abfs_abfs.init\t"<<time_during(t_start)<<std::endl;
gettimeofday( &t_start, NULL);
	m_abfs_abfs.init_radial( abfs_ccp, abfs );
ofs_mpi<<"TIME@ m_abfs_abfs.init_radial\t"<<time_during(t_start)<<std::endl;
//gettimeofday( &t_start, NULL);
//	m_abfs_abfs.init_radial_table();
//ofs_mpi<<"TIME@ m_abfs_abfs.init_radial_table\t"<<time_during(t_start)<<std::endl;

gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init( 1, this->kmesh_times, 1 );
ofs_mpi<<"TIME@ m_abfslcaos_lcaos.init\t"<<time_during(t_start)<<std::endl;
gettimeofday( &t_start, NULL);
	m_abfslcaos_lcaos.init_radial( abfs_ccp, lcaos, lcaos );
ofs_mpi<<"TIME@ m_abfslcaos_lcaos.init_radial\t"<<time_during(t_start)<<std::endl;
//gettimeofday( &t_start, NULL);
//	m_abfslcaos_lcaos.init_radial_table();
//ofs_mpi<<"TIME@ m_abfslcaos_lcaos.init_radial_table\t"<<time_during(t_start)<<std::endl;

	Born_von_Karman_period = ModuleBase::Vector3<int>{GlobalC::kv.nmp[0],GlobalC::kv.nmp[1],GlobalC::kv.nmp[2]};
ofs_mpi<<"TIME@ Exx_Lcao::init\t"<<time_during(t_start_all)<<std::endl;
ofs_mpi.close();

	auto PCA_test = [&]()
	{
		Exx_Abfs::Matrix_Orbs11 m_abfs_abfs;
		m_abfs_abfs.init( 2, this->kmesh_times, 1 );
		m_abfs_abfs.init_radial( abfs, abfs );
		m_abfs_abfs.init_radial_table();

		Exx_Abfs::Matrix_Orbs21 m_abfslcaos_lcaos;
		m_abfslcaos_lcaos.init( 1, this->kmesh_times, 1 );
		m_abfslcaos_lcaos.init_radial( abfs, lcaos, lcaos );
		m_abfslcaos_lcaos.init_radial_table();
		
		pthread_rwlock_t rwlock_Cw;	pthread_rwlock_init(&rwlock_Cw,NULL);
		pthread_rwlock_t rwlock_Vw;	pthread_rwlock_init(&rwlock_Vw,NULL);

		std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> Cws;
		std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> Vws;
		std::shared_ptr<ModuleBase::matrix> C = Abfs::DPcal_C(
			0,
			0,
			{0,0,0},
			m_abfs_abfs,
			m_abfslcaos_lcaos,
			index_abfs,
			index_lcaos,
			0,
			true,
			rwlock_Cw,
			rwlock_Vw,
			Cws,
			Vws);
		C->print(std::cout, 1E-10)<<std::endl;
		
		pthread_rwlock_destroy(&rwlock_Cw);
		pthread_rwlock_destroy(&rwlock_Vw);
	};

	auto overlap_test = [&]()
	{
		const auto lcaos = Exx_Abfs::Construct_Orbs::change_orbs( GlobalC::ORB, 1 );
		Exx_Abfs::Matrix_Orbs11 m_lcaos_lcaos;
		m_lcaos_lcaos.init(1,1,1);
		m_lcaos_lcaos.init_radial(lcaos,lcaos);
		m_lcaos_lcaos.init_radial_table();
		const ModuleBase::matrix m_overlap = m_lcaos_lcaos.cal_overlap_matrix( 0,0, {0,0,0},{0,0,0}, index_lcaos,index_lcaos );

		std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> lcaos_ccp;
		Conv_Coulomb_Pot::cal_orbs_ccp( lcaos, lcaos_ccp );
		Exx_Abfs::Matrix_Orbs11 m_lcaos_ccp;
		m_lcaos_ccp.init(1,1,1);
		m_lcaos_ccp.init_radial(lcaos,lcaos_ccp);
		m_lcaos_ccp.init_radial_table();
		const ModuleBase::matrix m_overlap_coulomb = m_lcaos_ccp.cal_overlap_matrix( 0,0, {0,0,0},{0,0,0}, index_lcaos,index_lcaos );

		std::ofstream ofs("matrix_overlap.dat");
		m_overlap.print(ofs, 1E-10)<<std::endl;
		m_overlap_coulomb.print(ofs, 1E-10)<<std::endl;
		ofs.close();
	};
}

void Exx_Lcao::cal_exx_ions(const Parallel_Orbitals &pv)
{
	ModuleBase::TITLE("Exx_Lcao","cal_exx_ions");
std::ofstream ofs_mpi(test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
timeval t_start, t_start_all;
gettimeofday( &t_start_all, NULL);

	auto cal_atom_centres_core = [](const std::vector<std::pair<size_t,size_t>> &atom_pairs_core) -> std::set<size_t>
	{
		std::set<size_t> atom_centres_core;
		for( const std::pair<size_t,size_t> & atom_pair : atom_pairs_core )
		{
			atom_centres_core.insert(atom_pair.first);
			atom_centres_core.insert(atom_pair.second);
		}
		return atom_centres_core;
	};
	
gettimeofday( &t_start, NULL);
#ifdef __MPI
	if(atom_pairs_core_origin.empty())
		switch(this->info.distribute_type)
		{
			case Exx_Lcao::Distribute_Type::Htime:
				atom_pairs_core_origin = Exx_Abfs::Parallel::Distribute::Htime::distribute( Born_von_Karman_period, info.ccp_rmesh_times );	break;
			case Exx_Lcao::Distribute_Type::Kmeans2:
				atom_pairs_core_origin = Exx_Abfs::Parallel::Distribute::Kmeans::distribute_kmeans2( MPI_COMM_WORLD );	break;
			case Exx_Lcao::Distribute_Type::Kmeans1:
				atom_pairs_core_origin = Exx_Abfs::Parallel::Distribute::Kmeans::distribute_kmeans1( MPI_COMM_WORLD, info.ccp_rmesh_times );	break;
			case Exx_Lcao::Distribute_Type::Order:
				atom_pairs_core_origin = Exx_Abfs::Parallel::Distribute::Order::distribute( info.ccp_rmesh_times );	break;
			default:
				throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));  break;
				//throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(static_cast<std::underlying_type<Exx_Lcao::Distribute_Type>::type>(info.distribute_type))+"\t"+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));	break;
		}
#endif
ofs_mpi<<"atom_pairs_core_origin\t"<<atom_pairs_core_origin.size()<<std::endl;
ofs_mpi<<"TIME@ Htime::distribute\t"<<time_during(t_start)<<std::endl;
//std::ofstream ofs_atom_pair("atom_pair+"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
//for( const auto & i : atom_pairs_core_origin )
//	ofs_atom_pair<<i.first<<"\t"<<i.second<<std::endl;
//ofs_atom_pair.close();
	
	#if TEST_EXX_LCAO==1
		std::ofstream ofs_adjs("adjs.dat");
		test_adjs(ofs_adjs);
		ofs_adjs.close();
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif

gettimeofday( &t_start, NULL);
	init_radial_table_ions( cal_atom_centres_core(atom_pairs_core_origin), atom_pairs_core_origin );
ofs_mpi<<"TIME@ init_radial_table_ions\t"<<time_during(t_start)<<std::endl;

gettimeofday( &t_start, NULL);
	Vs = Abfs::cal_Vs( atom_pairs_core_origin, m_abfs_abfs, index_abfs, info.ccp_rmesh_times, info.v_threshold, Vws );
ofs_mpi<<"TIME@ Abfs::cal_Vs\t"<<time_during(t_start)<<std::endl;
	Abfs::delete_empty_ptrs( Vws );
gettimeofday( &t_start, NULL);
	Vps = Abfs::cal_mps( Born_von_Karman_period, Vs );
ofs_mpi<<"TIME@ Abfs::cal_Vps\t"<<time_during(t_start)<<std::endl;
	atom_pairs_core = Abfs::get_atom_pair(Vps);
ofs_mpi<<"atom_pairs_core\t"<<atom_pairs_core.size()<<std::endl;

	const std::set<size_t> atom_centres_core = cal_atom_centres_core(atom_pairs_core);
ofs_mpi<<"atom_centres_core\t"<<atom_centres_core.size()<<std::endl;

gettimeofday( &t_start, NULL);
	H_atom_pairs_core = Abfs::get_H_pairs_core( atom_pairs_core );
ofs_mpi<<"H_atom_pairs_core\t"<<H_atom_pairs_core.size()<<std::endl;
ofs_mpi<<"TIME@ Exx_Lcao::allocate_Hexx\t"<<time_during(t_start)<<std::endl;
	
gettimeofday( &t_start, NULL);
	Cs = Abfs::cal_Cs( atom_centres_core, m_abfs_abfs,m_abfslcaos_lcaos, index_abfs,index_lcaos, info.c_threshold, Cws,Vws );
ofs_mpi<<"TIME@ Abfs::cal_Cs\t"<<time_during(t_start)<<std::endl;
	Abfs::delete_empty_ptrs( Cws );
gettimeofday( &t_start, NULL);
	Cps = Abfs::cal_mps( Born_von_Karman_period, Cs );
ofs_mpi<<"TIME@ Abfs::cal_Cps\t"<<time_during(t_start)<<std::endl;

gettimeofday( &t_start, NULL);
	schwarz.init( info.schwarz_threshold, info.schwarz_threshold );
	schwarz.cal_max_pair_fock( atom_centres_core, m_abfs_abfs,m_abfslcaos_lcaos, index_abfs,index_lcaos, Born_von_Karman_period, Cws,Vws );
ofs_mpi<<"TIME@ schwarz::init\t"<<time_during(t_start)<<std::endl;
gettimeofday( &t_start, NULL);
	cauchy.init( info.cauchy_threshold, info.cauchy_threshold, Born_von_Karman_period );
	cauchy.cal_norm_C_max( Cps, index_lcaos, index_abfs );
	cauchy.cal_norm_V( Vps );
ofs_mpi<<"TIME@ cauchy::init\t"<<time_during(t_start)<<std::endl;

	#if EXX_DM==2
	DM_para.init( H_atom_pairs_core, info.dm_threshold );
	#elif EXX_DM==3
	DM_para.allreduce.init( MPI_COMM_WORLD, Abfs::get_H_pairs_core_group( atom_pairs_core ), pv);
	#endif
	
	#if EXX_H_COMM==2
	Hexx_para.allreduce2.init(MPI_COMM_WORLD, H_atom_pairs_core, pv);
	#endif
	
ofs_mpi<<"TIME@ Exx_Lcao::cal_exx_ions\t"<<time_during(t_start_all)<<std::endl;

ofs_mpi<<"sizeof_Cs:\t"<<get_sizeof(Cs)<<std::endl;
ofs_mpi<<"sizeof_Vs:\t"<<get_sizeof(Vs)<<std::endl;
ofs_mpi<<"sizeof_Cps:\t"<<get_sizeof(Cps)<<std::endl;
ofs_mpi<<"sizeof_Vps:\t"<<get_sizeof(Vps)<<std::endl;
ofs_mpi<<"sizeof_Cws:\t"<<get_sizeof(Cws)<<std::endl;
ofs_mpi<<"sizeof_Vws:\t"<<get_sizeof(Vws)<<std::endl;

ofs_mpi.close();

	#if TEST_EXX_LCAO==1
		ofs_matrixes(test_dir.matrix+"Cws_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),Cws);
		ofs_matrixes(test_dir.matrix+"Vws_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),Vws);
		ofs_matrixes(test_dir.matrix+"Cs_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),Cs);
		ofs_matrixes(test_dir.matrix+"Cps_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),Cps);
		ofs_matrixes(test_dir.matrix+"Vs_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),Vs);
		ofs_matrixes(test_dir.matrix+"Vps_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),Vps);
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
}

void Exx_Lcao::cal_exx_elec(Local_Orbital_Charge &loc, complex<double>*** wfc_k_grid)
{
	ModuleBase::TITLE("Exx_Lcao","cal_exx_elec");

static int istep=0;
	#if TEST_EXX_LCAO==1
	//	ofs_matrixes("Cws_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_before.dat",Cws);
	//	ofs_matrixes("Vws_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_before.dat",Vws);
	//	ofs_matrixes("Cs_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_before.dat",Cs);
	//	ofs_matrixes("Vs_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_before.dat",Vs);
	//	ofs_matrixes("Vps_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_before.dat",Vps);
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif

//	if( exx_lcao.cal_DM_delta() < exx_lcao.get_DM_threshold() )	break;

std::ofstream ofs_mpi(test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
timeval t_start, t_start_all;
gettimeofday( &t_start_all, NULL);

#if EXX_DM==1
gettimeofday( &t_start, NULL);
	this->DM_para.cal_DM( Born_von_Karman_period, H_atom_pairs_core, info.dm_threshold, loc.DM, loc.DM_R, wfc_k_grid,);
ofs_mpi<<"TIME@ Exx_Lcao::cal_DM\t"<<time_during(t_start)<<std::endl;
#elif EXX_DM==2
gettimeofday( &t_start, NULL);
	if(!GlobalV::GAMMA_ONLY_LOCAL)
		this->DM_para.cal_DM_k( Born_von_Karman_period, H_atom_pairs_core, info.dm_threshold );
ofs_mpi<<"TIME@ Exx_Lcao::cal_DM\t"<<time_during(t_start)<<std::endl;
#elif EXX_DM==3
gettimeofday( &t_start, NULL);
	this->DM_para.cal_DM(info.dm_threshold, loc);
ofs_mpi<<"TIME@ Exx_Lcao::cal_DM\t"<<time_during(t_start)<<std::endl;
#endif

gettimeofday( &t_start, NULL);
#ifdef __MPI
	cauchy.cal_norm_D_max( DM_para.DMr );
#endif
ofs_mpi<<"TIME@ cauchy.cal_norm_D_max\t"<<time_during(t_start)<<std::endl;
#ifdef __MPI
ofs_mpi<<"sizeof_DM\t"<<get_sizeof(DM_para.DMr)<<std::endl;
#endif
gettimeofday( &t_start, NULL);
	// HexxR[is][iat1][iat2][box2]
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> HexxR = cal_Hexx();
ofs_mpi<<"TIME@ Exx_Lcao::cal_Hexx\t"<<time_during(t_start)<<std::endl;

ofs_mpi<<"sizeof_HexxR\t"<<get_sizeof(HexxR)<<std::endl;

gettimeofday( &t_start, NULL);
	this->energy = cal_energy(HexxR);
ofs_mpi<<"TIME@ Exx_Lcao::cal_energy\t"<<time_during(t_start)<<std::endl;

gettimeofday( &t_start, NULL);
#ifdef __MPI
	Hexx_para.Rexx_to_Km2D(*loc.ParaV, HexxR, {GlobalC::pot.init_chg=="file",GlobalC::CHR.out_chg} );
#endif
ofs_mpi<<"TIME@ Hexx_para.Rexx_to_Km2D\t"<<time_during(t_start)<<std::endl;

#ifdef __MPI
ofs_mpi<<"sizeof_Hexx2D\t"<<get_sizeof(Hexx_para.HK_Gamma_m2D)+get_sizeof(Hexx_para.HK_K_m2D)<<std::endl;
#endif

ofs_mpi<<"TIME@ Exx_Lcao::cal_exx_elec\t"<<time_during(t_start_all)<<std::endl;
ofs_mpi.close();

#ifdef __MPI
	auto print_Hexxk = [&]()
	{
		std::ofstream ofs("Hexxk_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
		for(int ik=0; ik!=Hexx_para.HK_K_m2D.size(); ++ik)
		{
			ofs<<"@\t"<<ik<<std::endl;
			Hexx_para.HK_K_m2D[ik].print(ofs, 1E-10, 1E-10)<<std::endl;
		};
		ofs.close();
	};
#endif

	#if TEST_EXX_LCAO==1
		ofs_matrixes("DMk_"+ModuleBase::GlobalFunc::TO_STRING(istep)+".dat",DM.DMk);
		ofs_matrixes("DMr_"+ModuleBase::GlobalFunc::TO_STRING(istep)+".dat",DM.DMr);
		ofs_matrixes("Hexx_"+ModuleBase::GlobalFunc::TO_STRING(istep)+".dat",Hexx);
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif

	auto print_LOCDM = [&]()
	{
		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			std::ofstream ofs("LOC.DM.dat",std::ofstream::app);
			const int it1=0, it2=0;
			for( size_t ia1=0; ia1!=GlobalC::ucell.atoms[it1].na; ++ia1 )
				for( size_t ia2=0; ia2!=GlobalC::ucell.atoms[it2].na; ++ia2 )
					for( size_t is=0; is!=GlobalV::NSPIN; ++is )
					{
						ofs<<"@\t"<<ia1<<"\t"<<ia2<<"\t"<<is<<std::endl;
						for( size_t iw1=0; iw1!=GlobalC::ucell.atoms[it1].nw; ++iw1 )
						{
							for( size_t iw2=0; iw2!=GlobalC::ucell.atoms[it2].nw; ++iw2 )
								ofs<<loc.DM[is][GlobalC::ucell.itiaiw2iwt(it1,ia1,iw1)][GlobalC::ucell.itiaiw2iwt(it2, ia2, iw2)]<<"\t";
							ofs<<std::endl;
						}
						ofs<<std::endl;
					}
			ofs.close();
		}
		else
		{
			throw logic_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
			/*
			static int istep=0;
			for( size_t is=0; is!=GlobalV::NSPIN; ++is )
			{
				std::ofstream ofs("LOC.DM_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(is));
				for(int T1=0; T1<GlobalC::ucell.ntype; T1++)
				{
					for(int I1=0; I1<GlobalC::ucell.atoms[T1].na; I1++)
					{
						const int iat1 = GlobalC::ucell.itia2iat(T1,I1);

						int iv=0;
						for( const auto & atom2 : adjs )
						{
							const size_t iat2 = atom2.first;
							for( const Abfs::Vector3_Order<int> & box2 : atom2.second )
							{
								ofs<<"@\t"<<iat1<<"\t"<<iat2<<"\t"<<box2<<std::endl;
								for( int iw1=0; iw1!=GlobalC::ucell.atoms[T1].nw; ++iw1 )
								{
									for( int iw2=0; iw2!=GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw; ++iw2 )
									{
										ofs<<loc.DM_R[is][GlobalC::LNNR.nlocstartg[iat1]+iv]<<"\t";
										++iv;
									}
									ofs<<std::endl;
								}
								ofs<<std::endl;
							}
						}
					}
				}
				ofs.close();
			}
			*/
		}
	};
	
	auto print_WFC = [&](std::complex<double>*** wfc_k_grid)
	{
		if( GlobalV::GAMMA_ONLY_LOCAL )
		{
			for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
			{
				std::ofstream ofs("GlobalC::LOWF.WFC_GAMMA_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(ik));
				for( size_t ib=0; ib!=GlobalV::NBANDS; ++ib )
				{
					for( size_t iwt=0; iwt!=GlobalV::NLOCAL; ++iwt )
					{
						//---------------------------------------------------------
						// GlobalC::LOWF.WFC_GAMMA has been replaced by wfc_dm_2d.cpp 
						// we need to fix this function in near future.
						// -- mohan add 2021-02-09
						//---------------------------------------------------------
						ModuleBase::WARNING_QUIT("Exx_Abfs::DM::cal_DMk_raw","need to update GlobalC::LOWF.WFC_GAMMA");
						//ofs<<GlobalC::LOWF.WFC_GAMMA[ik][ib][iwt]<<"\t";
					}
					ofs<<std::endl;
				}
				ofs.close();
			}
		}
		else
		{
			for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
			{
				std::ofstream ofs("LOWF.wfc_k_grid_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(ik));
				for( size_t ib=0; ib!=GlobalV::NBANDS; ++ib )
				{
					for( size_t iwt=0; iwt!=GlobalV::NLOCAL; ++iwt )
						ofs<<wfc_k_grid[ik][ib][iwt]<<"\t";
					ofs<<std::endl;
				}
				ofs.close();
			}
		}
	};

#ifdef __MPI
	auto print_Hexx=[&]()		// Peize Lin test 2019-11-14
	{
		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{		
				std::ofstream ofs("Hexx_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(is)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
				this->Hexx_para.HK_Gamma_m2D[is].print(ofs, 1E-10)<<std::endl;
			}
		}
		else
		{
			for(int ik=0; ik<GlobalC::kv.nks; ++ik)
			{
				std::ofstream ofs("Hexx_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(ik)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
				this->Hexx_para.HK_K_m2D[ik].print(ofs, 1E-10, 1E-10)<<std::endl;
			}
		}
	};
#endif

	auto print_wfc=[&](std::vector<ModuleBase::matrix>& wfc_gamma,
        std::vector<ModuleBase::ComplexMatrix>& wfc_k)		// Peize Lin test 2019-11-14
	{
		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{		
				std::ofstream ofs("wfc_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(is)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
				wfc_gamma[is].print(ofs, 1E-10)<<std::endl;
			}
		}
		else
		{
			for(int ik=0; ik<GlobalC::kv.nks; ++ik)
			{
				std::ofstream ofs("wfc_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(ik)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
				wfc_gamma[ik].print(ofs, 1E-10)<<std::endl;
			}
		}
	};

	auto print_ekb=[&]()		// Peize Lin test 2019-11-14
	{
		for(int ik=0; ik<GlobalC::kv.nks; ++ik)
		{
			std::ofstream ofs("ekb_"+ModuleBase::GlobalFunc::TO_STRING(ik)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK), std::ofstream::app);
			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
				ofs<<GlobalC::wf.ekb[ik][ib]<<"\t";
			ofs<<std::endl;
		}
	};	
	
	#if TEST_EXX_LCAO==1
	//	ofs_matrixes("Cws_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_end.dat",Cws);
	//	ofs_matrixes("Vws_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_end.dat",Vws);
	//	ofs_matrixes("Cs_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_end.dat",Cs);
	//	ofs_matrixes("Vs_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_end.dat",Vs);
	//	ofs_matrixes("Vps_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_end.dat",Vps);
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
++istep;

//	std::cout<<"screen.schwarz"<<std::endl;
//	std::cout<<Exx_Abfs::Screen::Schwarz::num_screen<<"\t"
//	    <<Exx_Abfs::Screen::Schwarz::num_cal<<std::endl;
//	std::cout<<"screen.cauchy"<<std::endl;
//	std::cout<<Exx_Abfs::Screen::Cauchy::num_screen1<<"\t"
//	    <<Exx_Abfs::Screen::Cauchy::num_screen2<<"\t"
//		<<Exx_Abfs::Screen::Cauchy::num_screen3<<"\t"
//		<<Exx_Abfs::Screen::Cauchy::num_cal<<std::endl;
}

void Exx_Lcao::cal_exx_elec_nscf(const Parallel_Orbitals &pv)
{
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> HexxR;
#ifdef __MPI
	Hexx_para.Rexx_to_Km2D(pv, HexxR, {GlobalC::pot.init_chg=="file",GlobalC::CHR.out_chg} );
#endif
}

/*
void Exx_Lcao::cal_Hexx_gamma( const std::set<std::pair<size_t,size_t>> &atom_pairs )
{
	for( const std::pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = GlobalC::ucell.iat2it[iat1];
		const size_t it2 = GlobalC::ucell.iat2it[iat2];

		const std::vector<Abfs::Atom_Info> adj1s = Abfs::get_adjs( GlobalC::ucell.atoms[it1].tau[GlobalC::ucell.iat2ia[iat2]] );
		const std::vector<Abfs::Atom_Info> adj2s = Abfs::get_adjs( GlobalC::ucell.atoms[it2].tau[GlobalC::ucell.iat2ia[iat2]] );
		for( const Abfs::Atom_Info & atom3 : adj1s )
		{
			const size_t iat3 = GlobalC::ucell.itia2iat( atom3.T, atom3.I );
			for( const Abfs::Atom_Info & atom4 : adj2s )
			{
				const size_t iat4 = GlobalC::ucell.itia2iat( atom4.T, atom4.I );
				const ModuleBase::matrix matrix_1342 = *Cs[iat1][iat3][atom3.box] * *Vps_gamma[iat1][iat2] * *Cs[iat2][iat4][atom4.box];

				for( size_t iw1=0; iw1!=GlobalC::ucell.atoms[it1].nw; ++iw1 )
					for( size_t iw2=0; iw2!=GlobalC::ucell.atoms[it2].nw; ++iw2 )
						for( size_t iw3=0; iw3!=GlobalC::ucell.atoms[atom3.T].nw; ++iw3 )
							for( size_t iw4=0; iw4!=GlobalC::ucell.atoms[atom4.T].nw; ++iw4 )
							{
								const double matrix_element_1342 = matrix_1342( iw1*GlobalC::ucell.atoms[atom3.T].nw+iw3, iw2*GlobalC::ucell.atoms[atom4.T].nw+iw4 );
								for( size_t is=0; is!=GlobalV::NSPIN; ++is )
								{
									Hexx_gamma[iat1][iat2][is](iw1,iw2) = matrix_element_1342 * DM_Gamma[iat3][iat4][is](iw3,iw4);
									Hexx_gamma[iat1][iat4][is](iw1,iw4) = matrix_element_1342 * DM_Gamma[iat3][iat2][is](iw3,iw2);
									Hexx_gamma[iat3][iat2][is](iw3,iw2) = matrix_element_1342 * DM_Gamma[iat1][iat4][is](iw1,iw4);
									Hexx_gamma[iat3][iat4][is](iw3,iw4) = matrix_element_1342 * DM_Gamma[iat1][iat2][is](iw1,iw2);
								}
							}
			}
		}
	}
}
*/

double Exx_Lcao::cal_energy( 
	const std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &HexxR ) const
{
	ModuleBase::TITLE("Exx_Lcao","cal_energy");
std::ofstream ofs_mpi(test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
timeval t_start;
gettimeofday( &t_start, NULL);

	double energy = 0;
	for( size_t is=0; is!=GlobalV::NSPIN; ++is )
	{
		for( const auto & HA : HexxR[is] )
		{
			const size_t iat1 = HA.first;
			for( const auto & HB : HA.second )
			{
				const size_t iat2 = HB.first;
				for( const auto & HC : HB.second )
				{
					const Abfs::Vector3_Order<int> & box2 = HC.first;
					if( const ModuleBase::matrix*const DMr_ptr = static_cast<const ModuleBase::matrix*const>(ModuleBase::GlobalFunc::MAP_EXIST( DM_para.DMr[is], iat1, iat2, box2 )) )
					{
						const ModuleBase::matrix & H = HC.second;
						assert(DMr_ptr->nr == H.nr);
						assert(DMr_ptr->nc == H.nc);
						energy += BlasConnector::dot( H.nr*H.nc, DMr_ptr->c,1, H.c,1 );
					}
				}
			}
		}
	}
	const std::map<int,double> SPIN_multiple = {{1,2}, {2,1}, {4,1}};							// ???
	energy *= SPIN_multiple.at(GlobalV::NSPIN);			// ?
	energy /= 2;					// /2 for Ry
	
ofs_mpi<<"TIME@ Exx_Lcao::cal_energy_cal\t"<<time_during(t_start)<<std::endl;
	Parallel_Reduce::reduce_double_all(energy);
ofs_mpi<<"E_exx\t"<<energy<<std::endl;
ofs_mpi.close();
	return energy;
}

void Exx_Lcao::add_Hexx( const size_t ik, const double alpha, LCAO_Matrix &lm) const
{
	ModuleBase::TITLE("Exx_Lcao","add_Hexx");
	
	if( GlobalV::GAMMA_ONLY_LOCAL )
	{
		const ModuleBase::matrix & H = Hexx_para.HK_Gamma_m2D[ik];
		for( size_t i=0; i<H.nr*H.nc; ++i )
		{
			lm.Hloc[i] += alpha * H.c[i];
		}
	}
	else
	{
		const ModuleBase::ComplexMatrix & H = Hexx_para.HK_K_m2D[ik];
		for( size_t i=0; i<H.nr*H.nc; ++i )
		{
			lm.Hloc2[i] += alpha * H.c[i];
		}
	}
}


void Exx_Lcao::init_radial_table_ions( const std::set<size_t> &atom_centres_core, const std::vector<std::pair<size_t,size_t>> &atom_pairs_core )
{
	ModuleBase::TITLE("Exx_Lcao::init_radial_table_ions");
	
std::ofstream ofs_mpi(test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
timeval t_start;

	std::map<size_t,std::map<size_t,std::set<double>>> radial_R;
	
	auto print_atom = [&](std::ostream &os)
	{
		os<<"atom_centres_core:"<<"\t";
		for(const auto &a : atom_centres_core)
		{
			os<<a<<"\t";
		}
		os<<std::endl;
		os<<"atom_pairs_core:"<<"\t";
		for(const auto &a : atom_pairs_core)
		{
			os<<"("<<a.first<<","<<a.second<<")"<<"\t";
		}
		os<<std::endl;
	};
	auto print_radial_R = [&](std::ostream &os)
	{
		os<<"radial_R"<<std::endl;
		for(const auto &rA : radial_R)
		{
			for(const auto &rB : rA.second)
			{
				os<<rA.first<<"\t"<<rB.first<<"\t"<<"|"<<"\t";
				for(const auto rC : rB.second)
				{
					os<<rC<<"\t";
				}
				os<<std::endl;
			}
		}
	};
	
	for(int it=0; it!=GlobalC::ucell.ntype; ++it)
	{
		radial_R[it][it].insert(0.0);
	}

	#if TEST_EXX_LCAO==1
		std::ofstream ofs_C("delta_R_C.dat");
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif

gettimeofday( &t_start, NULL);
	for( const int iat1 : atom_centres_core )
	{
		const size_t it1 = GlobalC::ucell.iat2it[iat1];
		const size_t ia1 = GlobalC::ucell.iat2ia[iat1];
		const ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[it1].tau[ia1];

		const std::map<size_t,std::vector<Abfs::Vector3_Order<int>>> adjs = Abfs::get_adjs(iat1);
		for( const auto & atom2 : adjs )
		{
			const size_t iat2 = atom2.first;
			const size_t it2 = GlobalC::ucell.iat2it[iat2];
			const size_t ia2 = GlobalC::ucell.iat2ia[iat2];
			const ModuleBase::Vector3<double> &tau2 = GlobalC::ucell.atoms[it2].tau[ia2];
			for( const Abfs::Vector3_Order<int> &box2 : atom2.second )
			{
				const double delta_R = (-tau1+tau2+(box2*GlobalC::ucell.latvec)).norm();
//				const double delta_R = (-tau1+tau2).norm();

				#if TEST_EXX_LCAO==1
					ofs_C<<-tau1+tau2+(box2*GlobalC::ucell.latvec)<<"\t"<<delta_R<<std::endl;
				#elif TEST_EXX_LCAO==-1
					#error "TEST_EXX_LCAO"
				#endif

				radial_R[it1][it2].insert( delta_R );
				radial_R[it2][it1].insert( delta_R );
			}
		}
	}
ofs_mpi<<"TIME@ radial_R_C\t"<<time_during(t_start)<<std::endl;

	#if TEST_EXX_LCAO==1
		ofs_C.close();
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif

gettimeofday( &t_start, NULL);
#if TEST_EXX_RADIAL>=1
	m_abfslcaos_lcaos.init_radial_table(radial_R);
#else
	m_abfslcaos_lcaos.init_radial_table();
#endif
ofs_mpi<<"TIME@ m_abfslcaos_lcaos.init_radial_table\t"<<time_during(t_start)<<std::endl;

	#if TEST_EXX_LCAO==1
		std::ofstream ofs_V("delta_R_V.dat");
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
gettimeofday( &t_start, NULL);
	std::vector<Abfs::Vector3_Order<int>> Coulomb_potential_boxes = Abfs::get_Coulomb_potential_boxes(info.ccp_rmesh_times);
	for( const std::pair<size_t,size_t> & atom_pair : atom_pairs_core )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = GlobalC::ucell.iat2it[iat1];
		const size_t ia1 = GlobalC::ucell.iat2ia[iat1];
		const size_t it2 = GlobalC::ucell.iat2it[iat2];
		const size_t ia2 = GlobalC::ucell.iat2ia[iat2];
		const ModuleBase::Vector3<double> &tau1 = GlobalC::ucell.atoms[it1].tau[ia1];
		const ModuleBase::Vector3<double> &tau2 = GlobalC::ucell.atoms[it2].tau[ia2];
		const double Rcut = std::min( GlobalC::ORB.Phi[it1].getRcut()*info.ccp_rmesh_times+GlobalC::ORB.Phi[it2].getRcut(), GlobalC::ORB.Phi[it1].getRcut()+GlobalC::ORB.Phi[it2].getRcut()*info.ccp_rmesh_times );

		for( const ModuleBase::Vector3<int> &box2 : Coulomb_potential_boxes )
		{
			const double delta_R = (-tau1+tau2+(box2*GlobalC::ucell.latvec)).norm();

			#if TEST_EXX_LCAO==1
				ofs_V<<-tau1+tau2+(box2*GlobalC::ucell.latvec)<<"\t"<<delta_R<<std::endl;
			#elif TEST_EXX_LCAO==-1
				#error "TEST_EXX_LCAO"
			#endif

			if( delta_R*GlobalC::ucell.lat0 < Rcut )
			{
				radial_R[it1][it2].insert( delta_R );
				radial_R[it2][it1].insert( delta_R );
			}
		}
	}
ofs_mpi<<"TIME@ radial_R_V\t"<<time_during(t_start)<<std::endl;
	#if TEST_EXX_LCAO==1
		ofs_V.close();
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
	
	#if TEST_EXX_LCAO==1
		std::ofstream ofs_V("delta_R_V.dat");
		for( const double & r : radial_R[0][0] )
			ofs_V<<r<<std::endl;
		ofs_V.close();
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif

gettimeofday( &t_start, NULL);
#if TEST_EXX_RADIAL>=1
	m_abfs_abfs.init_radial_table(radial_R);
#else
	m_abfs_abfs.init_radial_table();
#endif
ofs_mpi<<"TIME@ m_abfs_abfs.init_radial_table\t"<<time_during(t_start)<<std::endl;
ofs_mpi.close();
}


std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> Exx_Lcao::cal_Hexx() const
{
#ifdef __MKL
    const int mkl_threads = mkl_get_max_threads();
	mkl_set_num_threads(std::max(1UL,mkl_threads/atom_pairs_core.size()));
#endif
	
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> HexxR(GlobalV::NSPIN);
#ifdef _OPENMP
	omp_lock_t Hexx_lock;
	omp_init_lock(&Hexx_lock);
#endif
	
#ifdef _OPENMP
	#pragma omp parallel
#endif
	{
		// m_new( i2, i1, i3 ) = m( i1, i2, i3 )
		auto transform = [](
			const ModuleBase::matrix & m,
			const size_t n1, const size_t n2, const size_t n3 ) -> ModuleBase::matrix
		{
			assert( n1*n2*n3 == m.nr*m.nc );
			const auto length = sizeof(double)*n3;
			const auto n13 = n1*n3;
			const double *m_ptr = m.c;
			ModuleBase::matrix m_new( n1*n2, n3, false );
			for( size_t i1=0; i1!=n1; ++i1 )
			{
				double *m_new_ptr = m_new.c+i1*n3;
				for( size_t i2=0; i2!=n2; ++i2, m_ptr+=n3, m_new_ptr+=n13 )
					memcpy( m_new_ptr, m_ptr, length );
			}
			return m_new;
		};

		// C = alpha * A.? * B.?
		auto gemm = [](
			const char transA, 		// A / A.T / A.C
			const char transB, 		// B / B.T / B.C
			const int M,			// A.?.nr = C.nr
			const int N,			// B.?.nc = C.nc
			const int K,			// A.?.nc = B.?.nr
			const double alpha,
			const ModuleBase::matrix & A,
			const ModuleBase::matrix & B) -> ModuleBase::matrix
		{
			ModuleBase::matrix C(M,N,false);
			BlasConnector::gemm(
				transA, transB,
				M, N, K,
				alpha,
				A.c, (transA=='N')?K:M,
				B.c, (transB=='N')?N:K,
				0,
				C.c, N );
			return C;
		};

		// insert m_tmp into m_all
		auto insert_matrixes = []( 
			std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> & m_all,
			std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> & m_tmp)
		{
			for( size_t is=0; is!=GlobalV::NSPIN; ++is )
			{
				for( auto & m_tmp1 : m_tmp[is] )
				{
					const size_t iat1 = m_tmp1.first;
					for( auto & m_tmp2 : m_tmp1.second )
					{
						const size_t iat2 = m_tmp2.first;
						for( auto & m_tmp3 : m_tmp2.second )
						{
							const Abfs::Vector3_Order<int> & box2 = m_tmp3.first;
							if( ModuleBase::matrix*const m_all_ptr = static_cast<ModuleBase::matrix*const>(ModuleBase::GlobalFunc::MAP_EXIST( m_all[is], iat1, iat2, box2 )) )
							{
								*m_all_ptr += m_tmp3.second;
							}
							else
							{
								m_all[is][iat1][iat2][box2] = std::move(m_tmp3.second);
							}
						}
					}
				}
				m_tmp[is].clear();
			}
		};

		auto vector_empty = []( const std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> & v ) -> bool
		{
			for( const auto &i : v )
			{
				if(!i.empty())	return false;
			}
			return true;
		};
		
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> HexxR_tmp(GlobalV::NSPIN);
	
#ifdef _OPENMP
		#pragma omp for
#endif
		for(size_t i_atom_pair=0; i_atom_pair<atom_pairs_core.size(); ++i_atom_pair)
		{		
			const size_t iat1 = atom_pairs_core[i_atom_pair].first;
			const size_t iat2 = atom_pairs_core[i_atom_pair].second;
			
			//ofs_thread<<iat1<<"\t"<<iat2<<std::endl;

			auto add_Hexx_one = [&HexxR_tmp, iat1, iat2, this](
				const ModuleBase::matrix &Hexx_one, const size_t is, const size_t iat_I, const size_t iat_J, const Abfs::Vector3_Order<int> &box_IJ )
			{
				if(iat1!=iat2)
				{
					const ModuleBase::matrix Hexx_one_T = transpose(Hexx_one);
					if( ModuleBase::matrix * const HexxR_ptr = static_cast<ModuleBase::matrix*const>(ModuleBase::GlobalFunc::MAP_EXIST( HexxR_tmp[is], iat_J, iat_I, Abfs::Vector3_Order<int>(-box_IJ)%this->Born_von_Karman_period )) )
						*HexxR_ptr += Hexx_one_T;
					else
						HexxR_tmp[is][iat_J][iat_I][Abfs::Vector3_Order<int>(-box_IJ)%this->Born_von_Karman_period] = std::move(Hexx_one_T);
				}	
				if( ModuleBase::matrix * const HexxR_ptr = static_cast<ModuleBase::matrix*const>(ModuleBase::GlobalFunc::MAP_EXIST( HexxR_tmp[is], iat_I, iat_J, Abfs::Vector3_Order<int>(box_IJ)%this->Born_von_Karman_period )) )
					*HexxR_ptr += Hexx_one;
				else
					HexxR_tmp[is][iat_I][iat_J][Abfs::Vector3_Order<int>(box_IJ)%this->Born_von_Karman_period] = std::move(Hexx_one);
			};
			
			for( const auto & Cp1s : Cps.at(iat1) )
			{
				const size_t iat3 = Cp1s.first;

				for( const auto & Cp2s : Cps.at(iat2) )
				{
					const size_t iat4 = Cp2s.first;

					const size_t it1 = GlobalC::ucell.iat2it[iat1];
					const size_t it2 = GlobalC::ucell.iat2it[iat2];
					const size_t it3 = GlobalC::ucell.iat2it[iat3];
					const size_t it4 = GlobalC::ucell.iat2it[iat4];

					for( const auto & Cp1 : Cp1s.second )
					{
						const Abfs::Vector3_Order<int> & box3 = Cp1.first;
						const auto & Cp13 = Cp1.second;

						for( const auto & Cp2 : Cp2s.second )
						{
							const Abfs::Vector3_Order<int> & box4 = Cp2.first;
							const auto & Cp24 = Cp2.second;
							
							if( schwarz.screen( iat1,iat2,iat3,iat4, box3,box4 ) )	continue;

							for( const auto &Vp : Vps.at(iat1).at(iat2) )
							{
								const Abfs::Vector3_Order<int> & box2 = Vp.first;

								const ModuleBase::matrix & C_13 = *Cp13,								// iw1*iw3, \mu1
											 & V    = *Vp.second,							// \mu1, \mu2
											 & C_24 = *Cp24;								// iw2*iw4, \mu2

								const Exx_Abfs::Screen::Cauchy::Info_Step info_step = cauchy.input_info( iat1, iat2, iat3, iat4, box2, box3, box4 );
								int cauchy_postcal;
								if(!( cauchy_postcal = cauchy.postcalA( info_step ) ))	continue;

	//								assert( V.nc==C_24.nc );
								const ModuleBase::matrix VC_24 = gemm(									// \mu1, iw2*iw4
									'N', 'T',
									index_abfs[it1].count_size,
									index_lcaos[it2].count_size * index_lcaos[it4].count_size,
									index_abfs[it2].count_size,
									1, V, C_24 );
								const ModuleBase::matrix VC_24_T = transform( VC_24,					// iw2*\mu1, iw4
									index_abfs[it1].count_size,
									index_lcaos[it2].count_size,
									index_lcaos[it4].count_size );
								if(!( cauchy_postcal = cauchy.postcalB( info_step, VC_24_T, index_lcaos[it2].count_size, index_abfs[it1].count_size, index_lcaos[it4].count_size, cauchy_postcal ) ))	continue;
													
								const ModuleBase::matrix C_13_T = transform( C_13,						// iw3*iw1, \mu1
									index_lcaos[it1].count_size,
									index_lcaos[it3].count_size,
									index_abfs[it1].count_size );
									 
								for( size_t is=0; is!=GlobalV::NSPIN; ++is )
								{
									switch( cauchy_postcal )		// Attention: case and go on calculating
									{
										case 4:
										{
											const ModuleBase::matrix * const DM_ptr = static_cast<const ModuleBase::matrix*const>(ModuleBase::GlobalFunc::MAP_EXIST( DM_para.DMr[is], iat3, iat4, Abfs::Vector3_Order<int>(box2-box3+box4)%Born_von_Karman_period ));
											if( DM_ptr )
											{
												const ModuleBase::matrix DVC_32 = gemm(								// iw3, \mu1*iw2
													'N', 'T',
													index_lcaos[it3].count_size,
													index_abfs[it1].count_size * index_lcaos[it2].count_size,
													index_lcaos[it4].count_size,
													1, 
													*DM_ptr, 
													VC_24 );
												if( cauchy.postcalC( info_step, DVC_32, index_lcaos[it3].count_size, index_abfs[it1].count_size, index_lcaos[it2].count_size, 1 ) )
												{
													const ModuleBase::matrix Hexx_12 = gemm(							// iw1, iw2
														'N', 'N',
														index_lcaos[it1].count_size,
														index_lcaos[it2].count_size,
														index_lcaos[it3].count_size * index_abfs[it1].count_size,
														-2, C_13, DVC_32 );
													if( cauchy.postcalD( Hexx_12 ) )
														add_Hexx_one(Hexx_12, is, iat1, iat2, box2);
												}
											}
										}	// end case 4
										case 3:
										{
											if( const ModuleBase::matrix * const DM_ptr = static_cast<const ModuleBase::matrix*const>(ModuleBase::GlobalFunc::MAP_EXIST( DM_para.DMr[is], iat1, iat4, Abfs::Vector3_Order<int>(box2+box4)%Born_von_Karman_period )) )
											{
												const ModuleBase::matrix DVC_12 = gemm(								// iw1, \mu1*iw2
													'N', 'T',
													index_lcaos[it1].count_size,
													index_abfs[it1].count_size * index_lcaos[it2].count_size,
													index_lcaos[it4].count_size,
													1, 
													*DM_ptr, 
													VC_24 );
												if( cauchy.postcalC( info_step, DVC_12, index_lcaos[it1].count_size, index_abfs[it1].count_size, index_lcaos[it2].count_size, 3 ) )
												{
													const ModuleBase::matrix Hexx_32 = gemm(							// iw3, iw2
														'N', 'N',
														index_lcaos[it3].count_size,
														index_lcaos[it2].count_size,
														index_lcaos[it1].count_size * index_abfs[it1].count_size,
														-2, C_13_T, DVC_12 );
													if( cauchy.postcalD( Hexx_32 ) )
														add_Hexx_one(Hexx_32, is, iat3, iat2, box2-box3);
												}
											}
										}	// end case 3
										case 2:
										{
											if( const ModuleBase::matrix * const DM_ptr = static_cast<const ModuleBase::matrix*const>(ModuleBase::GlobalFunc::MAP_EXIST( DM_para.DMr[is], iat3, iat2, Abfs::Vector3_Order<int>(box2-box3)%Born_von_Karman_period )) )
											{
												const ModuleBase::matrix DVC_34 = gemm(								// iw3, \mu1*iw4
													'N', 'N',
													index_lcaos[it3].count_size,
													index_abfs[it1].count_size * index_lcaos[it4].count_size,
													index_lcaos[it2].count_size,
													1, 
													*DM_ptr, 
													VC_24_T );
												if( cauchy.postcalC( info_step, DVC_34, index_lcaos[it3].count_size, index_abfs[it1].count_size, index_lcaos[it4].count_size, 1 ) )
												{
													const ModuleBase::matrix Hexx_14 = gemm(							// iw1, iw4
														'N', 'N',
														index_lcaos[it1].count_size,
														index_lcaos[it4].count_size,
														index_lcaos[it3].count_size * index_abfs[it1].count_size,
														-2, C_13, DVC_34 );
													if( cauchy.postcalD( Hexx_14 ) )
														add_Hexx_one(Hexx_14, is, iat1, iat4, box2+box4);
												}
											}
										}	// end case 2
										case 1:
										{
											if( const ModuleBase::matrix * const DM_ptr = static_cast<const ModuleBase::matrix*const>(ModuleBase::GlobalFunc::MAP_EXIST( DM_para.DMr[is], iat1, iat2, Abfs::Vector3_Order<int>(box2)%Born_von_Karman_period )) )
											{
												const ModuleBase::matrix DVC_14 = gemm(								// iw1, \mu1*iw4
													'N', 'N',
													index_lcaos[it1].count_size,
													index_abfs[it1].count_size * index_lcaos[it4].count_size,
													index_lcaos[it2].count_size,
													1, 
													*DM_ptr, 
													VC_24_T );
												if( cauchy.postcalC( info_step, DVC_14, index_lcaos[it1].count_size, index_abfs[it1].count_size, index_lcaos[it4].count_size, 3 ) )
												{
													const ModuleBase::matrix Hexx_34 = gemm(							// iw3, iw4
														'N', 'N',
														index_lcaos[it3].count_size,
														index_lcaos[it4].count_size,
														index_lcaos[it1].count_size * index_abfs[it1].count_size,
														-2, C_13_T, DVC_14 );
													if( cauchy.postcalD( Hexx_34 ) )
														add_Hexx_one(Hexx_34, is, iat3, iat4, box2-box3+box4);
												}
											}
										}	// end case 1
										case 0: ;
									}	// end switch cauchy_postcal
								}	// end for is
							}	// end for box2
						}	// end for box4
					}	// end for box3
					
#ifdef _OPENMP
					if( !vector_empty(HexxR_tmp) && omp_test_lock(&Hexx_lock) )
					{
						insert_matrixes(HexxR,HexxR_tmp);
						omp_unset_lock(&Hexx_lock);
					}
#else
					if( !vector_empty(HexxR_tmp) )
					{
						insert_matrixes(HexxR,HexxR_tmp);
					}
#endif
				}	// end for iat4
			}	// end for iat3
		}	// end omp for i_atom_pair

#ifdef _OPENMP
		if(!vector_empty(HexxR_tmp))
		{
			omp_set_lock(&Hexx_lock);
			insert_matrixes(HexxR,HexxR_tmp);
			omp_unset_lock(&Hexx_lock);
		}
#else
		if(!vector_empty(HexxR_tmp))
		{
			insert_matrixes(HexxR,HexxR_tmp);
		}
#endif
	} // end omp parallel

#ifdef _OPENMP
	omp_destroy_lock(&Hexx_lock);
#endif
	
#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
#endif
	
	return HexxR;
}
#endif //__MPI