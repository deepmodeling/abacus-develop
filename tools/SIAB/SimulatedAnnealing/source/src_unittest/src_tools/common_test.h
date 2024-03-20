#ifndef COMMON_TEST
#define COMMON_TEST

#include<iostream>
#include<fstream>
#include<iomanip>
using std::cout;
using std::endl;

#include<sys/time.h>

extern double sum_t;

// gettimeofday( &time_now, NULL);
inline double get_time(const timeval &begin, const timeval &end)
{
	return (double)(end.tv_sec-begin.tv_sec) + (double)(end.tv_usec-begin.tv_usec)/1000000.0;
}

inline void cout_matrix( const ComplexMatrix &matrix)
{
	cout<< std::setprecision(10)<< std::setiosflags(ios::fixed);
	const size_t nr(matrix.nr), nc(matrix.nc);
	for(size_t i(0); i<nr; ++i)
	{
		for( size_t j(0); j<nc; ++j)
		{
			cout<<matrix(i,j)<<"\t";
		}
		cout<<endl;
	}
	cout<<endl;	
	cout<< std::resetiosflags(ios::fixed);	
	return;
}

inline void ofs_matrix( const char * file_name, const ComplexMatrix &matrix)
{
	ofstream ofs(file_name);
	ofs<< std::setprecision(10)<< std::setiosflags(ios::fixed);
	const size_t nr(matrix.nr), nc(matrix.nc);
	for(size_t i(0); i<nr; ++i)
	{
		for( size_t j(0); j<nc; ++j)
		{
			ofs<<matrix(i,j)<<"\t";
		}
		ofs<<endl;
	}
	ofs<<endl;	
	ofs<< std::resetiosflags(ios::fixed);
	ofs.close();
	return;
}

#endif