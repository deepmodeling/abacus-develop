#ifndef EXX_LCAO_TEST_H
#define EXX_LCAO_TEST_H

#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "complexmatrix-test.h"
#include "matrix-test.h"
#include <map>
#include <string>
#include<memory>
#include<fstream>
#include<iomanip>
#include<sys/time.h>
using namespace std;

inline double time_during ( const timeval &t_start, const timeval &t_end )
{ 
	return (double)(t_end.tv_sec-t_start.tv_sec) + (double)(t_end.tv_usec-t_start.tv_usec)/1000000.0;
};
inline double time_during ( const timeval &t_start )
{ 
	timeval t_end;
	gettimeofday( &t_end, NULL );
	return time_during(t_start,t_end);
};
inline double time_cut(timeval &t)
{
	const double t_delta = time_during(t);
	gettimeofday(&t, NULL );
	return t_delta;
}



inline size_t get_sizeof( const ModuleBase::matrix & m ){	return sizeof(double)*m.nr*m.nc; }
inline size_t get_sizeof( const ModuleBase::ComplexMatrix & m ){	return sizeof(double)*m.nr*m.nc; }
inline size_t get_sizeof( const std::shared_ptr<ModuleBase::matrix> & m ){	return sizeof(ModuleBase::matrix)+sizeof(double)*m->nr*m->nc; }
inline size_t get_sizeof( const std::shared_ptr<ModuleBase::ComplexMatrix> & m ){	return sizeof(ModuleBase::ComplexMatrix)+sizeof(double)*m->nr*m->nc; }
inline size_t get_sizeof( const std::weak_ptr<ModuleBase::matrix> & m ){	return sizeof(ModuleBase::matrix)+sizeof(double)*m.lock()->nr*m.lock()->nc; }
inline size_t get_sizeof( const std::weak_ptr<ModuleBase::ComplexMatrix> & m ){	return sizeof(ModuleBase::ComplexMatrix)+sizeof(double)*m.lock()->nr*m.lock()->nc; }
template<typename T> static size_t get_sizeof( const std::vector<T> & v )
{
	size_t length = sizeof(T)*v.size();
	for( const T & i : v )	length+=get_sizeof(i);
	return length;
}
template<typename T1,typename T2> static size_t get_sizeof( const std::map<T1,T2> & m )
{
	size_t length = (sizeof(T1)+sizeof(T2))*m.size();
	for( const auto & i : m )	length+=get_sizeof(i.second);
	return length;
}

static void ofs_matrixes( const std::string & file_name, const ModuleBase::matrix & ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	if(flag_print_content)
	{
		ms.print(ofs, 1E-10)<<std::endl;
	}
	ofs.close();
}
static void ofs_matrixes( const std::string & file_name, const ModuleBase::ComplexMatrix & ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	if(flag_print_content)
		ms.print(ofs, 1E-10, 1E-10)<<std::endl;	
	ofs.close();
}

template<typename T>
static void ofs_matrixes( const std::string & file_name, const std::vector<T> & ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	for( size_t i=0; i<ms.size(); ++i )
	{
			ofs<<fixed;
			const auto precision_old = ofs.precision(15);
		ofs<<"@\t"<<i<<std::endl;
			ofs.precision(precision_old);
			ofs.unsetf(std::ostream::floatfield);
		ofs_matrixes( file_name, ms[i], flag_print_content );		
	}
	ofs.close();
}

template<typename T1,typename T2>
static void ofs_matrixes( const std::string & file_name, const std::map<T1,T2> &ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	for( const auto &m : ms )
	{
			ofs<<fixed;
			const auto precision_old = ofs.precision(15);
		ofs<<"@\t"<<m.first<<std::endl;
			ofs.precision(precision_old);
			ofs.unsetf(std::ostream::floatfield);
		ofs_matrixes( file_name, m.second, flag_print_content );		
	}
	ofs.close();
}

template<typename T>
static void ofs_matrixes( const std::string & file_name, const std::weak_ptr<T> & ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	if(flag_print_content)
		ofs<<*ms.lock()<<std::endl;
	ofs.close();
}
template<typename T>
static void ofs_matrixes( const std::string & file_name, const std::shared_ptr<T> & ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	if(flag_print_content)
		ofs<<*ms<<std::endl;
	ofs.close();
}


/*
template<typename T1, typename T2, typename T3, typename M>
static void ofs_matrixes( const std::string & file_name, const std::map<T1,std::map<T2,std::map<T3,std::weak_ptr<M>>>> & ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	for( const auto &m1 : ms )
	{
		const auto i1 = m1.first;
		for( const auto &m2 : m1.second )
		{
			const auto i2 = m2.first;
			for( const auto &m3 : m2.second )
			{
				const auto i3 = m3.first;
					ofs<<fixed;
					const auto precision_old = ofs.precision(15);
				ofs<<"@\t"<<i1<<"\t"<<i2<<"\t"<<i3<<std::endl;
					ofs.precision(precision_old);
					ofs.unsetf(std::ostream::floatfield);
				if(flag_print_content)
					ofs<<*m3.second.lock()<<std::endl;
			}
		}
	}
	ofs.close();
}

template<typename T1, typename T2, typename T3, typename M>
static void ofs_matrixes( const std::string & file_name, const std::map<T1,std::map<T2,std::map<T3,std::shared_ptr<M>>>> & ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	for( const auto &m1 : ms )
	{
		const auto i1 = m1.first;
		for( const auto &m2 : m1.second )
		{
			const auto i2 = m2.first;
			for( const auto &m3 : m2.second )
			{
				const auto i3 = m3.first;
					ofs<<fixed;
					const auto precision_old = ofs.precision(15);
				ofs<<"@\t"<<i1<<"\t"<<i2<<"\t"<<i3<<std::endl;
					ofs.precision(precision_old);
					ofs.unsetf(std::ostream::floatfield);
				if(flag_print_content)
					ofs<<*m3.second<<std::endl;
			}
		}
	}
	ofs.close();
}

template<typename T1, typename T2, typename M>
static void ofs_matrixes( const std::string & file_name, const std::map<T1,std::map<T2,std::vector<M>>> & ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	for( const auto &m1 : ms )
	{
		const auto i1 = m1.first;
		for( const auto &m2 : m1.second )
		{
			const auto i2 = m2.first;
			for( size_t i3=0; i3!=m2.second.size(); ++i3 )
			{
					ofs<<fixed;
					const auto precision_old = ofs.precision(15);
				ofs<<"@\t"<<i1<<"\t"<<i2<<"\t"<<i3<<std::endl;
					ofs.precision(precision_old);
					ofs.unsetf(std::ostream::floatfield);
				if(flag_print_content)
					ofs<<m2.second[i3]<<std::endl;
			}
		}
	}
	ofs.close();
}

template<typename T1, typename T2, typename T3, typename M>
static void ofs_matrixes( const std::string & file_name, const std::map<T1,std::map<T2,std::map<T3,std::vector<M>>>> & ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	for( const auto &m1 : ms )
	{
		const auto i1 = m1.first;
		for( const auto &m2 : m1.second )
		{
			const auto i2 = m2.first;
			for( const auto &m3 : m2.second )
			{
				const auto i3 = m3.first;
				for( size_t i4=0; i4!=m3.second.size(); ++i4 )
				{
						ofs<<fixed;
						const auto precision_old = ofs.precision(15);
					ofs<<"@\t"<<i1<<"\t"<<i2<<"\t"<<i3<<"\t"<<i4<<std::endl;
						ofs.precision(precision_old);
						ofs.unsetf(std::ostream::floatfield);
					if(flag_print_content)
						ofs<<m3.second[i4]<<std::endl;
				}
			}
		}
	}
	ofs.close();
}

template<typename T1, typename T2, typename T3, typename M>
static void ofs_matrixes( const std::string & file_name, const std::map<T1,std::map<T2,std::map<T3,M>>> & ms, const bool flag_print_content=true )
{
	std::ofstream ofs(file_name,std::ofstream::app);
	for( const auto &m1 : ms )
	{
		const auto i1 = m1.first;
		for( const auto &m2 : m1.second )
		{
			const auto i2 = m2.first;
			for( const auto &m3 : m2.second )
			{
				const auto i3 = m3.first;
					ofs<<fixed;
					const auto precision_old = ofs.precision(15);
				ofs<<"@\t"<<i1<<"\t"<<i2<<"\t"<<i3<<std::endl;
					ofs.precision(precision_old);
					ofs.unsetf(std::ostream::floatfield);
				if(flag_print_content)
					ofs<<m3.second<<std::endl;
			}
		}
	}
	ofs.close();
}
*/

#endif