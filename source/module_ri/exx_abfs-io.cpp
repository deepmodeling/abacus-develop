#include <fstream>
#include <stdexcept>

#include "exx_abfs-io.h"
#include "exx_abfs-jle.h"
#include "exx_abfs-abfs_index.h"
#include "../module_hamilt_pw/hamilt_pwdft/global.h"
#include "../module_basis/module_ao/ORB_read.h"
#include "../module_base/global_function.h"
#include "../module_base/math_integral.h" // mohan add 2021-04-03


std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> Exx_Abfs::IO::construct_abfs(
	const LCAO_Orbitals &orbs,
	const std::vector<std::string> &files_abfs,
	const double kmesh_times )
{
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs( files_abfs.size() );
	for( size_t T=0; T!=files_abfs.size(); ++T )
		abfs[T] = construct_abfs_T( 
			files_abfs[T],
			T,
			static_cast<int>(orbs.get_kmesh() * kmesh_times) | 1,			// Nk must be odd
//			orbs.get_dk() / kmesh_times,
			orbs.get_dk(),								// Peize Lin change 2017-04-16
			orbs.get_dr_uniform() );		
	
	return abfs;
}

std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> Exx_Abfs::IO::construct_abfs( 
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> & abfs_pre,
	const LCAO_Orbitals &orbs,
	const std::vector<std::string> &files_abfs,
	const double kmesh_times )
{
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> 
		&&abfs = construct_abfs( orbs, files_abfs, kmesh_times );
			
	assert( abfs.size() == abfs_pre.size() );
	for( size_t T=0; T!=abfs.size(); ++T )
	{
		if( abfs[T].size() < abfs_pre[T].size() )
			abfs[T].resize( abfs_pre[T].size() );
		for( size_t L=0; L!=abfs_pre[T].size(); ++L )
		{
			abfs[T][L].insert( abfs[T][L].begin(), abfs_pre[T][L].begin(), abfs_pre[T][L].end() );
		}
	}	
	
	return abfs;
}

std::vector<std::vector<Numerical_Orbital_Lm>> Exx_Abfs::IO::construct_abfs_T( 
	const std::string & file_name,
	const int &T,
	const int &nk,
	const double &dk,
	const double &dr_uniform)
{
	std::string label;
	size_t L_size;
	std::map<size_t,size_t> N_size;
	size_t meshr;
	double dr;
	std::map<size_t,std::map<size_t,std::vector<double>>> psis;
	
	/*----------------------
	  1.read abfs
	----------------------*/
	std::string word;
	
	std::ifstream ifs( file_name.c_str() );
	if(!ifs)
		throw std::runtime_error(" Can't find the abfs ORBITAL file.");	
	
	while( ifs.good() )
	{
		ifs >> word;
		
		if( "Element"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, label );
		}
		else if ( "Lmax"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, L_size );
			
			if( L_size>=9 )
			{
				std::stringstream ss;
				ss<<"Lmax>=9 error in "<<__FILE__<<" line "<<__LINE__;
				throw std::invalid_argument(ss.str());
			}
		}
		else if ( "Sorbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[0] );
		}
		else if ( "Porbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[1] );
		}
		else if ( "Dorbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[2] );
		}
		else if ( "Forbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[3] );
		}
		else if ( "Gorbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[4] );
		}
		else if ( "Horbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[5] );
		}
		else if ( "Iorbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[6] );
		}
		else if ( "Jorbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[7] );
		}
		else if ( "Korbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[8] );
		}
		else if ( "END"==word )
		{
			break;
		}		
	}
	
	ModuleBase::CHECK_NAME(ifs, "Mesh");
	ifs >> meshr;
	
	ModuleBase::CHECK_NAME(ifs, "dr");
	ifs >> dr;

	while(ifs.good())
	{
		ifs >> word;
		if(word=="Type")
		{
			std::string s_L, s_N;
			ifs >> s_L >> s_N;
			
			size_t T,L,N;
			ifs >> T >> L >> N;
			
			psis[L][N].resize(meshr);
			for(int ir=0; ir!=meshr; ir++)
			{
				ifs >> psis[L][N][ir];
			}		
		}

	}
	ifs.close();

	 
	/*----------------------
	  2.check L,N orbital
	----------------------*/
	for( size_t L=0; L<=L_size; ++L )
		if( N_size.find(L) == N_size.end() )
		{
			std::stringstream ss;
			ss<<"Can't find N of L="<<L<<" in "<<file_name;
			throw std::domain_error(ss.str());
		}
	for( size_t L=0; L<=L_size; ++L )
		for( size_t N=0; N!=N_size[L]; ++N )
			if( psis.find(L)==psis.end() || psis[L].find(N)==psis[L].end() )
			{
				std::stringstream ss;
				ss<<"Can't find abf of L="<<L<<" T="<<T<<" in "<<file_name;
				throw std::domain_error(ss.str());
			}

			
	/*----------------------
	  3.rab, radial
	----------------------*/		
	if(meshr%2==0)	++meshr;
	
	std::vector<double> rab(meshr);
	std::vector<double> radial(meshr);
	for( int ir=0; ir!=meshr; ++ir )
	{
		rab[ir] = dr;
		radial[ir] = ir*dr; //mohan 2010-04-19
	}

	
	/*----------------------
	  4.normalize psi
	----------------------*/
	for( size_t L=0; L<=L_size; ++L )
	{	
		for( size_t N=0; N!=N_size[L]; ++N )
		{
			std::vector<double> psir(meshr);
			std::vector<double> inner(meshr);
			psis[L][N].resize(meshr);
			for( int ir=0; ir!=meshr; ++ir )
			{
				psir[ir] = psis[L][N][ir] * radial[ir];
				inner[ir] = psir[ir] * psir[ir];
			}
			double unit = 0.0;	
			ModuleBase::Integral::Simpson_Integral(meshr, ModuleBase::GlobalFunc::VECTOR_TO_PTR(inner), ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab), unit);
			for( int ir=0; ir!=meshr; ++ir )
			{
				psis[L][N][ir] /= sqrt(unit);
			}
		}
	}

	
	/*----------------------
	  5.construct abfs
	----------------------*/	
	std::vector<std::vector<Numerical_Orbital_Lm>> abfs_T;

	abfs_T.resize(L_size+1);
	for( size_t L=0; L<=L_size; ++L )
	{
		abfs_T[L].resize(N_size[L]);
		for( size_t N=0; N!=N_size[L]; ++N )
		{
			abfs_T[L][N].set_orbital_info(
				label,
				T, //type
				L, //angular momentum L
				N, // number of orbitals of this L
				meshr, // number of radial mesh
				ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab),
				ModuleBase::GlobalFunc::VECTOR_TO_PTR(radial),// radial mesh value(a.u.)
				Numerical_Orbital_Lm::Psi_Type::Psi,
				ModuleBase::GlobalFunc::VECTOR_TO_PTR(psis[L][N]), // radial wave function
				nk,
				dk,
				dr_uniform,
				false,
				true, GlobalV::CAL_FORCE);		
		}
	}
	
	return abfs_T;
}
