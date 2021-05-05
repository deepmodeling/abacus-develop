/* matrix.cpp file */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <limits>

using namespace std;
#include "matrix.h"

#ifdef __NORMAL
#else
#include "lapack_connector.h"
#endif

//*********************************************************
// The init() function is the main initialization routine.
// Sets up sizes and allocates memory for matrix class.
// All constructors call init()
// ********************************************************

//int matrix::mCount = 0;

void matrixAlloc()
{
// mohan add 2021-04-25
#ifdef __NORMAL
	cout << "Allocation error for Matrix" << endl;
	exit(0);
#else
	WARNING_QUIT("matrix","Allocation error for Matrix");
#endif
}

matrix::matrix( const int nrows, const int ncols, const bool flag_zero )
	:nr(nrows),
	 nc(ncols),
	 c(nullptr)
{
	if( nr && nc )
	{
		auto handler_old = set_new_handler(matrixAlloc);
		c = new double[nr*nc];
		set_new_handler(handler_old);
		if(flag_zero)	this->zero_out();
	}
}

matrix::matrix( const matrix &m_in )
	:nr(m_in.nr),
	 nc(m_in.nc),
	 c(nullptr)
{
	if( nr && nc )
	{
		auto handler_old = set_new_handler(matrixAlloc);
		c = new double[nr*nc];
		set_new_handler(handler_old);
		memcpy( c, m_in.c, nr*nc*sizeof(double) );
	}
}

// Peize Lin add 2016-08-05
matrix::matrix( matrix && m_in )
	:nr(m_in.nr),
	 nc(m_in.nc)
{
	c = m_in.c;
	m_in.nr = m_in.nc = 0;
	m_in.c = nullptr;
}

// Peize Lin change 2018-07-02
matrix& matrix::operator=( const matrix & m_in )
{
	this->create( m_in.nr, m_in.nc, false );
	memcpy( c, m_in.c, nr*nc*sizeof(double) );
	return *this;
}

// Peize Lin add 2016-08-05
matrix& matrix::operator=( matrix && m_in )
{
	nr = m_in.nr;		nc = m_in.nc;
	if(c)	delete[] c;
	c = m_in.c;
	m_in.nr = m_in.nc = 0;
	m_in.c = nullptr;
	return *this;
}

/*
double & matrix::operator()(const int ir,const int ic)
{
	assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
	return c[ir*nc+ic];
}

const double & matrix::operator()(const int ir,const int ic) const
{
	assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
	return c[ir*nc+ic];
}
*/

//*************
//
// destructor
//
//*************
matrix::~matrix()
{
	if(c)					// Peize Lin add 2016-08-05
	{
		delete [] c;
		c = nullptr;
	}
}

//******************************
// reallocate memory for matrix
//******************************
// Peize Lin change 2018-07-29
void matrix::create( const int nrow, const int ncol, const bool flag_zero )
{
	if( nrow && ncol )
	{
		if(c)
		{
			const int size=nrow*ncol;
			if( size!=nr*nc )
			{
				delete[] c;
				auto handler_old = set_new_handler(matrixAlloc);			
				c = new double[size];
				set_new_handler(handler_old);
			}
		}
		else
		{
			auto handler_old = set_new_handler(matrixAlloc);
			c = new double[nrow * ncol];
			set_new_handler(handler_old);
		}			
			
		nr = nrow;
		nc = ncol;
		if(flag_zero)	zero_out();				// Peize Lin change 2018-03-12
	}
	else
	{
		if(c)	delete[] c;
		c = nullptr;
		nr = nrow;
		nc = ncol;
	}
}

/* Adding matrices, as a friend */
matrix operator+(const matrix &m1, const matrix &m2)
{
	assert(m1.nr == m2.nr);
	assert(m2.nc == m2.nc);

	matrix tm(m1);
	const int size = m1.nr*m1.nc;
	for (int i = 0; i < size; i++) 
		tm.c[i] += m2.c[i];
	return tm;
}

/* Subtracting matrices, as a friend */
matrix operator-(const matrix &m1, const matrix &m2)
{
	assert(m1.nr == m2.nr);
	assert(m2.nc == m2.nc);

	matrix tm(m1);
	const int size = m1.nr*m1.nc;
	for(int i = 0; i < size; i++) 
		tm.c[i] -= m2.c[i];
	return tm;
}

//***************************************
//
// Multiplying matrices, as a friend 
//
// *************************************
matrix operator*(const matrix &m1, const matrix &m2)
{
	// fixed bug 2010-01-26
	assert(m1.nc == m2.nr);
	
    // allocate the result and zero it out
    matrix mprod( m1.nr, m2.nc, false );

#ifdef __NORMAL
	mprod.zero_out();
    // do the multiply and return
    for (int i = 0;i < m1.nr;i++)
	{
        for (int j = 0;j < m2.nc;j++)
		{
            for (int k = 0;k < m1.nc;k++)
			{
                mprod(i, j) += m1(i, k) * m2(k, j);
			}
		}
	}
#else
	// Peize Lin accelerate 2017-10-27
	LapackConnector::gemm(
		'N', 'N', 
		m1.nr, m2.nc, m1.nc,
		1, m1.c, m1.nc, m2.c, m2.nc, 
		0, mprod.c, mprod.nc);
#endif

	return mprod;
}

/* Scale a matrix */
matrix operator*(const double &s, const matrix &m)
{
	matrix sm(m);
	const int size=m.nr*m.nc;
	for (int i = 0; i < size; i++) 
		sm.c[i] *= s;
	return sm;
}

/* matrix * double */
matrix operator*(const matrix &m,const double &s)
{
	matrix sm(m);
	const int size=m.nr*m.nc;
	for (int i = 0; i < size; i++)
		sm.c[i] *= s;
	return sm;
}

/* Scale a matrix in place */
void matrix::operator*=(const double &s)
{
	const int size=nc*nr;
	for (int i = 0; i < size; i++) 
		c[i] *= s;
}

/* Accumulate to a matrix in place */
void matrix::operator+=(const matrix & m)
{
	assert( nr==m.nr );
	assert( nc==m.nc );
	const int size=nc*nr;
	const double * const c_in = m.c;
	for( int i = 0; i < size; ++i ) 
		c[i] += c_in[i];
}


/* decumulate to a matrix in place */
void matrix::operator-=(const matrix & m)
{
	assert( nr==m.nr );
	assert( nc==m.nc );
	const int size=nc*nr;
	const double * const c_in = m.c;
	for( int i = 0; i < size; ++i ) 
		c[i] -= c_in[i];
}

/* zero out the matrix */
void matrix::zero_out(void)
{
	const int size = nr*nc;
	for(int i = 0; i < size; i++)
		c[i] = 0.0;
}


matrix transpose(const matrix &m)
{
	matrix tm( m.nc, m.nr, false );
	for (int i = 0;i < m.nr;i++)
		for (int j = 0;j < m.nc;j++)
			tm(j, i) = m(i, j);
	return tm;
}

double matrix::trace_on(void) const
{
    assert(nr == nc);
    int inch = nc + 1;
    int size = nr * nc;
    double tr = 0.0;
    for (int i = 0; i < size; i += inch)
    {
        tr += c[i];
    }
    return tr;
}

void matrix::get_extreme_eigen_values(double &ev_lower, double &ev_upper)const
{
    double *a = new double[nr];
    double *b = new double[nr];
    for (int i = 0; i < nr; ++i)
    {
        double sum = 0.0;
        for(int j = 0; j < nc; ++j)
        {
            sum += fabs(c[i * nc + j]);
        }
        sum -= fabs(c[i * nc + i]);
        a[i] = c[i * nc + i] - sum;
        b[i] = c[i * nc + i] + sum;
    }

    ev_lower = a[0];
    ev_upper = b[0];

    for (int i = 1; i < nr; ++i)
    {
        if (a[i] < ev_lower) ev_lower = a[i];
        if (b[i] > ev_upper) ev_upper = b[i];
    }
    delete[] a;
    delete[] b;
}

// Peize Lin add 2017-05-27
void matrix::reshape( const double nr_new, const double nc_new )
{
	assert( nr*nc == nr_new*nc_new );
	nr=nr_new;
	nc=nc_new;
}

double trace_on(const matrix &A, const matrix &B)
{
    double tr = 0.0;
    for (int i = 0; i < A.nr; ++i)
        for (int k = 0; k < A.nc; ++k)
            tr += A(i,k) * B(k, i);
    return tr;
}

double mdot(const matrix &A, const matrix &B)
{
    assert (A.nr == B.nr);
    assert (A.nc == B.nc);
    const int size = A.nr * A.nc;

    double sum = 0.0;
    for (int i = 0; i < size; ++i)
        sum += A.c[i] * B.c[i];
    return sum;
}

// Peize Lin add 2016-09-08
std::ostream & operator<<( std::ostream & os, const matrix & m )
{
	for( int ir=0; ir!=m.nr; ++ir )
	{
		for( int ic=0; ic!=m.nc; ++ic )
		{
			if(abs(m(ir,ic))>1E-10)
				os<<m(ir,ic)<<"\t";
			else
				os<<0<<"\t";
		}
		os<<std::endl;
	}	
	return os;
}

// Peize Lin add 2016-09-08
double matrix::max() const
{
	double value = std::numeric_limits<double>::min();
	const int size = nr * nc;
	for( int i=0; i<size; ++i )
		value = std::max( value, c[i] );
	return value;
}

// Peize Lin add 2016-09-08
double matrix::min() const
{
	double value = std::numeric_limits<double>::max();
	const int size = nr * nc;
	for( int i=0; i<size; ++i )
	{
		value = std::min( value, c[i] );
	}
	return value;
}

// Peize Lin add 2018-07-02
double matrix::absmax() const
{
	double value = 0;
	const int size = nr * nc;
	for( int i=0; i<size; ++i )
	{
		value = std::max( value, std::abs(c[i]) );
	}
	return value;
}

double matrix::norm() const
{
// mohan add 2021-04-25, no tests.
#ifdef  __NORMAL
	double nn = 0.0;
	for(int i=0; i<nr*nc; ++i)
	{
		nn += c[i]*c[i];
	}	
	return sqrt(nn);
#else
	return LapackConnector::nrm2(nr*nc,c,1);
#endif
}
