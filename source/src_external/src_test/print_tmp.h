#ifndef PRINT_TMP_H
#define PRINT_TMP_H

#include<iostream>
#include<fstream>
using namespace std;

template<typename T>
static void print_tmp(const T &S)
{
	for( const auto &s1 : S)
	{
		const auto i1=s1.first;
		for( const auto &s2 : s1.second )
		{
			const auto i2=s2.first;
			for( const auto &s3 : s2.second )
			{
				const auto i3=s3.first;
				std::cout<<i1<<"\t"<<i2<<"\t"<<i3<<std::endl;
			}
		}
	}
}

void print_matrixess(const std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &ms, std::ofstream &ofs)
{
	for( const auto & c1 : ms )
		for( const auto & c2 : c1.second )
			for( const auto & c3 : c2.second )
			{
				ofs<<c1.first<<"\t"<<c2.first<<"\t"<<c3.first<<std::endl;
				c3.second.lock()->print(ofs, 1E-10)<<std::endl;
			}
}

void print_matrixess(const std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> &ms, std::ofstream &ofs)
{
	for( const auto & c1 : ms )
		for( const auto & c2 : c1.second )
			for( const auto & c3 : c2.second )
			{
				ofs<<c1.first<<"\t"<<c2.first<<"\t"<<c3.first<<std::endl;
				c3.second->print(ofs, 1E-10)<<std::endl;
			}
}

#endif