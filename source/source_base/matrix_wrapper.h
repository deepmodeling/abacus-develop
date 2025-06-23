//==========================================================
// AUTHOR : Peize Lin
// DATE : 2018-07-31
//==========================================================

#ifndef MATRIX_WRAPPER_TIANHE2

#ifndef MATRIX_WRAPPER_H
#define MATRIX_WRAPPER_H

#include "matrix.h"
#include <cstring>

// !!! Attention: c is very dangerous, may be changed by other class,
//	e.g. changed the value, size!=nr*nc, deleted
namespace ModuleBase
{

class Matrix_Wrapper
{
public:
	int nr;
	int nc;
	double *c;
	bool flag_delete_c;

	Matrix_Wrapper(): nr(0), nc(0), c(nullptr), flag_delete_c(false){}
	Matrix_Wrapper( const matrix &m ): nr(m.nr), nc(m.nc), c(m.c), flag_delete_c(false){}
	inline void create( const int nr_in, const int nc_in, const bool flag_zero );
	inline matrix to_matrix();
	~Matrix_Wrapper(){ if(flag_delete_c) delete[]c; }
	Matrix_Wrapper( const Matrix_Wrapper &m )=delete;
	Matrix_Wrapper( Matrix_Wrapper &m )=delete;
	inline Matrix_Wrapper( Matrix_Wrapper &&m );
	inline Matrix_Wrapper&operator=( const Matrix_Wrapper&m );
	inline Matrix_Wrapper&operator=( Matrix_Wrapper&&m );
};


inline void Matrix_Wrapper::create( const int nr_in, const int nc_in, const bool flag_zero )
{
	nr = nr_in;	nc = nc_in;
	if(flag_delete_c)
		delete[] c;
	c = new double[nr*nc];
	flag_delete_c = true;
	if(flag_zero)
		memset( c, 0, sizeof(double)*nr*nc );
}

inline matrix Matrix_Wrapper::to_matrix()
{
	assert( flag_delete_c==true );
	flag_delete_c = false;
	matrix m;
	m.nr = nr;	m.nc = nc;
	m.c = c;
	return m;
}

inline Matrix_Wrapper::Matrix_Wrapper( Matrix_Wrapper &&m )
	:nr(m.nr),
	 nc(m.nc),
	 c(m.c),
	 flag_delete_c(m.flag_delete_c)
{
	m.nr = m.nc = 0;
	m.c = nullptr;
	m.flag_delete_c = false;
}

inline Matrix_Wrapper& Matrix_Wrapper::operator=( const Matrix_Wrapper&m )
{
	assert( !m.flag_delete_c );
	nr = m.nr;	nc = m.nc;	c = m.c;		flag_delete_c = m.flag_delete_c;
	return *this;
}

inline Matrix_Wrapper& Matrix_Wrapper::operator=( Matrix_Wrapper&&m )
{
	nr = m.nr;	nc = m.nc;	c = m.c;		flag_delete_c = m.flag_delete_c;
	m.nr = 0;	m.nc = 0;	m.c = nullptr;	m.flag_delete_c = false;
	return *this;
}

}

#endif

#else
#include "matrix_wrapper_tianhe2.h"
#endif
