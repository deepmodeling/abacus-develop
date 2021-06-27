#ifndef MATRIX_INL_H
#define MATRIX_INL_H

#include "../module_base/matrix.h"
#include <cstring>

// test
#include <iostream>
using namespace std;

inline matrix::matrix( const int nrows, const int ncols, const bool flag_zero )
	:nr(nrows), nc(ncols), c(nullptr)
{
cout<<__FILE__<<__LINE__<<endl;
	if( nr && nc )
	{
//		auto handler_old = set_new_handler(matrixAlloc);
		c = new double[nr*nc];
//		set_new_handler(handler_old);
		if(flag_zero)	this->zero_out();
	}
}

inline matrix::matrix( const matrix &m_in )
{
cout<<__FILE__<<__LINE__<<endl;
	create( m_in.nr, m_in.nc, false );
cout<<__FILE__<<__LINE__<<endl;
	memcpy( c, m_in.c, nr*nc*sizeof(double) );
cout<<__FILE__<<__LINE__<<endl;
}

// Peize Lin add 2016-08-05
inline matrix::matrix( matrix && m_in )
{
cout<<__FILE__<<__LINE__<<endl;
	nr = m_in.nr; nc = m_in.nc;
	c = m_in.c;
	m_in.nr = m_in.nc = 0;
	m_in.c = nullptr;
}

// Peize Lin change 2018-07-02
inline matrix& matrix::operator=( const matrix & m_in )
{
cout<<__FILE__<<__LINE__<<endl;
	this->create( m_in.nr, m_in.nc, false );
cout<<__FILE__<<__LINE__<<endl;
	memcpy( c, m_in.c, nr*nc*sizeof(double) );
cout<<__FILE__<<__LINE__<<endl;
	return *this;
}

// Peize Lin add 2016-08-05
inline matrix& matrix::operator=( matrix && m_in )
{
cout<<__FILE__<<__LINE__<<endl;
	nr = m_in.nr;		nc = m_in.nc;
	if(c)	delete[] c;
	c = m_in.c;
	m_in.nr = m_in.nc = 0;
	m_in.c = nullptr;
cout<<__FILE__<<__LINE__<<endl;
	return *this;
}

inline double & matrix::operator()(const int ir,const int ic)
{
	assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
	return c[ir*nc+ic];
}

inline const double & matrix::operator()(const int ir,const int ic) const
{
	assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
	return c[ir*nc+ic];
}

#endif