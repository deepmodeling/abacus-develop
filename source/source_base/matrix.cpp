/* matrix.cpp file */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <cmath>
#include <limits>

#include "matrix.h"

#ifdef __NORMAL
#else
#include "module_external/blas_connector.h"
#endif

//*********************************************************
// The init() function is the main initialization routine.
// Sets up sizes and allocates memory for matrix class.
// All constructors call init()
// ********************************************************

//int matrix::mCount = 0;
namespace ModuleBase
{

/*create a matrix with nrows*ncols size */
matrix::matrix( const int nrows, const int ncols, const bool flag_zero )
    :nr(nrows),
     nc(ncols),
     c(nullptr)
{
    const std::size_t size = matrix_size(nr, nc);
    if( size )
    {
        c = new double[size];
        if(flag_zero)	this->zero_out();
    }
}

matrix::matrix( const matrix &m_in )
    :nr(m_in.nr),
     nc(m_in.nc),
     c(nullptr)
{
    const std::size_t size = matrix_size(nr, nc);
    if( size )
    {
        c = new double[size];
        memcpy( c, m_in.c, size*sizeof(double) );
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
    // `create` leaves `c` null for an empty matrix, and memcpy's arguments are declared
    // non-null even for a zero count -- so assigning an empty matrix is undefined behaviour.
    const std::size_t size = matrix_size(nr, nc);
    if( size )
    {
        memcpy( c, m_in.c, size*sizeof(double) );
    }
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
    const std::size_t new_size = matrix_size(nrow, ncol);
    if( new_size == 0 )
    {
        if(c)	delete[] c;
        c = nullptr;
        nr = nrow;
        nc = ncol;
        return;
    }

    const std::size_t old_size = ( c != nullptr ) ? matrix_size(nr, nc) : 0;
    if( c == nullptr || new_size != old_size )
    {
        // Allocate before releasing the old buffer so a bad_alloc leaves the
        // matrix unchanged and the caller can still report the failure.
        double* new_data = new double[new_size];
        delete[] c;
        c = new_data;
    }

    nr = nrow;
    nc = ncol;
    if(flag_zero)	zero_out();				// Peize Lin change 2018-03-12
}

/* Unary minus*/
matrix operator-(const matrix& m1)
{
    matrix tm(m1);
    const std::size_t size = matrix_size(m1.nr, m1.nc);
    for (std::size_t i = 0; i < size; i++)
        tm.c[i] = -tm.c[i];
    return tm;
}

/* Adding matrices, as a friend */
matrix operator+(const matrix &m1, const matrix &m2)
{
    assert(m1.nr == m2.nr);
    assert(m1.nc == m2.nc);

    matrix tm(m1);
    const std::size_t size = matrix_size(m1.nr, m1.nc);
    for (std::size_t i = 0; i < size; i++) 
    {
        tm.c[i] += m2.c[i];
    }
    return tm;
}

/* Subtracting matrices, as a friend */
matrix operator-(const matrix &m1, const matrix &m2)
{
    assert(m1.nr == m2.nr);
    assert(m1.nc == m2.nc);

    matrix tm(m1);
    const std::size_t size = matrix_size(m1.nr, m1.nc);
    for(std::size_t i = 0; i < size; i++) 
    {
        tm.c[i] -= m2.c[i];
    }
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
    BlasConnector::gemm(
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
    const std::size_t size = matrix_size(m.nr, m.nc);
    for (std::size_t i = 0; i < size; i++) 
    {
        sm.c[i] *= s;
    }
    return sm;
}

/* matrix * double */
matrix operator*(const matrix &m,const double &s)
{
    matrix sm(m);
    const std::size_t size = matrix_size(m.nr, m.nc);
    for (std::size_t i = 0; i < size; i++)
    {
        sm.c[i] *= s;
    }
    return sm;
}

/* Scale a matrix in place */
void matrix::operator*=(const double &s)
{
    const std::size_t size = matrix_size(nr, nc);
    for (std::size_t i = 0; i < size; i++) 
        c[i] *= s;
}

/* Accumulate to a matrix in place */
void matrix::operator+=(const matrix & m)
{
    if(nr == 0 || nc == 0)
    {
        return;
    }
    assert( nr==m.nr );
    assert( nc==m.nc );
    const std::size_t size = matrix_size(nr, nc);
    const double * const c_in = m.c;
    for( std::size_t i = 0; i < size; ++i ) 
    {
        c[i] += c_in[i];
    }
}


/* decumulate to a matrix in place */
void matrix::operator-=(const matrix & m)
{
    if(nr == 0 || nc == 0)
    {
        return;
    }
    assert( nr==m.nr );
    assert( nc==m.nc );
    const std::size_t size = matrix_size(nr, nc);
    const double * const c_in = m.c;
    for( std::size_t i = 0; i < size; ++i ) 
    {
        c[i] -= c_in[i];
    }
}

/* zero out the matrix */
void matrix::zero_out(void)
{
    const std::size_t size = matrix_size(nr, nc);
    if(size == 0)
    {
        return;
    }
    for(std::size_t i = 0; i < size; i++)
        c[i] = 0.0;
}

/* fill out the matrix with a constant double*/
void matrix::fill_out(const double x)
{
    const std::size_t size = matrix_size(nr, nc);
    if(size == 0)
    {
        return;
    }
    for(std::size_t i = 0; i < size; i++)
        c[i] = x;
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
    const std::size_t size = matrix_size(nr, nc);
    const std::size_t inch = static_cast<std::size_t>(nc) + 1;
    double tr = 0.0;
    for (std::size_t i = 0; i < size; i += inch)
    {
        tr += c[i];
    }
    return tr;
}

// Peize Lin add 2017-05-27
void matrix::reshape( const int nr_new, const int nc_new, const bool flag_zero )
{
    const std::size_t old_size = matrix_size(nr, nc);
    const std::size_t new_size = matrix_size(nr_new, nc_new);
    assert( old_size == new_size );
    nr=nr_new;
    nc=nc_new;

    if(flag_zero) this-> zero_out();
}

double trace_on(const matrix &A, const matrix &B)
{
    assert(A.nr == B.nc);
    assert(A.nc == B.nr);

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
    const std::size_t size = matrix_size(A.nr, A.nc);

    double sum = 0.0;
    for (std::size_t i = 0; i < size; ++i)
        sum += A.c[i] * B.c[i];
    return sum;
}

/*
// Peize Lin add 2016-09-08
std::ostream & operator<<( std::ostream & os, const matrix & m )
{
    for( int ir=0; ir!=m.nr; ++ir )
    {
        for( int ic=0; ic!=m.nc; ++ic )
        {
            if(std::abs(m(ir,ic))>1E-10)
                os<<m(ir,ic)<<"\t";
            else
                os<<0<<"\t";
        }
        os<<std::endl;
    }	
    return os;
}
*/

// Peize Lin add 2021.09.08
/*output the element in a matrix whose absolute value is larger than the input threshold*/
std::ostream & matrix::print( std::ostream & os, const double threshold ) const
{
    for( int ir=0; ir!=this->nr; ++ir )
    {
        for( int ic=0; ic!=this->nc; ++ic )
        {
            if(std::abs((*this)(ir,ic))>threshold)
                os<<(*this)(ir,ic)<<"\t";
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
    const std::size_t size = matrix_size(nr, nc);
    for( std::size_t i=0; i<size; ++i )
        value = std::max( value, c[i] );
    return value;
}

// Peize Lin add 2016-09-08
double matrix::min() const
{
    double value = std::numeric_limits<double>::max();
    const std::size_t size = matrix_size(nr, nc);
    for( std::size_t i=0; i<size; ++i )
    {
        value = std::min( value, c[i] );
    }
    return value;
}

// Peize Lin add 2018-07-02
double matrix::absmax() const
{
    double value = 0;
    const std::size_t size = matrix_size(nr, nc);
    for( std::size_t i=0; i<size; ++i )
    {
        value = std::max( value, std::abs(c[i]) );
    }
    return value;
}

double matrix::norm() const
{
// mohan add 2021-04-25, no tests.
    const std::size_t size = matrix_size(nr, nc);
#ifdef  __NORMAL
    double nn = 0.0;
    for(std::size_t i=0; i<size; ++i)
    {
        nn += c[i]*c[i];
    }	
    return sqrt(nn);
#else
    // BLAS takes an int length.  Chunk the vector so a matrix with more than
    // INT_MAX elements is never passed as a negative or wrapped length.
    double nn = 0.0;
    const std::size_t chunk_limit = static_cast<std::size_t>(1) << 30;
    for (std::size_t offset = 0; offset < size; offset += chunk_limit)
    {
        const std::size_t remaining = size - offset;
        const std::size_t chunk_size = remaining > chunk_limit ? chunk_limit : remaining;
        const int chunk = static_cast<int>(chunk_size);
        const double partial = BlasConnector::nrm2(chunk, c + offset, 1);
        nn += partial * partial;
    }
    return sqrt(nn);
#endif
}

}
