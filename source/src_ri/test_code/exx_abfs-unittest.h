#ifndef EXX_ABFS_UNITTEST_H
#define EXX_ABFS_UNITTEST_H

#include<iostream>
#include<sstream>
using namespace std;

#include "src_ri/exx_abfs.h"

#include "matrix-test.h"
#include "exx_abfs-unittest.h"

static void of_abfs_cpp( const std::string &file_name, const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb ) 
{
	for( int i=0; i!=orb.size(); ++i )
		for( int j=0; j!=orb[i].size(); ++j )
			for( int k=0; k!=orb[i][j].size(); ++k )
			{
				std::stringstream ss;	ss<<file_name<<"_"<<i<<"_"<<j<<"_"<<k;
				std::ofstream ofs(ss.str().c_str());
				for( const auto & c : orb[i][j][k].psi_uniform )
					ofs<<c<<std::endl;
				ofs.close();
			}
}

static void cout_matrix_4( const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &matrix )
{
	for( const auto & m1 : matrix )
	{
		const size_t TA = m1.first;
		{
			for( const auto & m2 : m1.second )
			{
				const size_t IA = m2.first; 
				for( const auto & m3 : m2.second )
				{
					const size_t TB = m3.first;
					for( const auto & m4 : m3.second )
					{
						const size_t IB = m4.first;
						const auto & m = m4.second;
						
						std::cout<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<std::endl;
						m.print(std::cout, 1E-10)<<std::endl;
					}
				}
			}
		}
	}
}

static void ofs_ms( 
	const std::string& file_name, 
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &ms )
{
	std::ofstream ofs(file_name.c_str());
	ofs<<"0000"<<std::endl;
	ms.at(0).at(0).at(0).at(0).print(ofs, 1E-10)<<std::endl;
	ofs<<"0001"<<std::endl;
	ms.at(0).at(0).at(0).at(1).print(ofs, 1E-10)<<std::endl;
	ofs.close();
}

static void ofs_ms( 
	const std::string& file_name, 
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &ms )
{
	std::ofstream ofs(file_name.c_str());
	ofs<<"0000"<<std::endl;
	for( const auto m : ms.at(0).at(0).at(0).at(0) )
		m.print(ofs, 1E-10)<<std::endl;
	ofs<<"0001"<<std::endl;
	for( const auto m : ms.at(0).at(0).at(0).at(1) )
		m.print(ofs, 1E-10)<<std::endl;
	ofs.close();
}

static void ofs_ms( 
	const std::string& file_name, 
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<std::vector<ModuleBase::matrix>>>>>> &ms )
{
	std::ofstream ofs(file_name.c_str());
	ofs<<"0000"<<std::endl;
	for( const auto m1 : ms.at(0).at(0).at(0).at(0) )
		for( const auto m2 : m1 )
			m2.print(ofs, 1E-10)<<std::endl;
	ofs<<"0001"<<std::endl;
	for( const auto m1 : ms.at(0).at(0).at(0).at(1) )
		for( const auto m2 : m1 )
			m2.print(ofs, 1E-10)<<std::endl;
	ofs.close();	
}

#endif
