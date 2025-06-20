//==========================================================
// AUTHOR : Peize Lin
// DATE : 2018-07-31
//==========================================================

#ifndef MATRIX_WRAPPER_H
#define MATRIX_WRAPPER_H

#include "matrix.h"
#include <cstring>
namespace ModuleBase
{

// !!! Attention: c is very dangerous, may point to somewhere deleted

class Matrix_Wrapper
{
public:
	double *c;
	int nr;
	int nc;
	
	Matrix_Wrapper(): nr(0), nc(0), c(nullptr){}
	Matrix_Wrapper( const matrix &m ): nr(m.nr), nc(m.nc), c(m.c){}
	inline void create( const int nr_in, const int nc_in, const bool flag_zero );
	inline matrix to_matrix();
	Matrix_Wrapper( const Matrix_Wrapper &m ): nr(m.nr), nc(m.nc), c(m.c){}
	Matrix_Wrapper( Matrix_Wrapper &m ): nr(m.nr), nc(m.nc), c(m.c){}
	Matrix_Wrapper( Matrix_Wrapper &&m ): nr(m.nr), nc(m.nc), c(m.c){}
	inline Matrix_Wrapper&operator=( const Matrix_Wrapper&m ){ nr=m.nr; nc=m.nc; c=m.c; return *this; };
	inline Matrix_Wrapper&operator=( Matrix_Wrapper&m ){ nr=m.nr; nc=m.nc; c=m.c; return *this; };
	inline Matrix_Wrapper&operator=( Matrix_Wrapper&&m ){ nr=m.nr; nc=m.nc; c=m.c; return *this; };
};


inline void Matrix_Wrapper::create( const int nr_in, const int nc_in, const bool flag_zero )
{
	nr = nr_in;	nc = nc_in;
	c = new double[nr*nc];
	if(flag_zero)
		memset( c, 0, sizeof(double)*nr*nc );
}

inline matrix Matrix_Wrapper::to_matrix()
{
	matrix m;
	m.nr = nr;	m.nc = nc;
	m.c = c;
	nr = nc = 0;
	c = nullptr;
	return m;
}

}
#endif
