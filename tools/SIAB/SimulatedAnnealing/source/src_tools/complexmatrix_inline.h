//==========================================================
// AUTHOR : Lixin He, Mohan Chen
// LAST UPDATE : 2009-03-23 modify "=" operator
//==========================================================

#ifndef COMPLEXMATRIX_INLINE_H
#define COMPLEXMATRIX_INLINE_H

#include <cassert>
#include <new>
#include <cstdlib>
#include <iostream>
#include <exception>
#include <cstring>
#include "complexmatrix.h"
#include "lapack_connector.h"
#include "complexarray.h"

// constructor with sizes
inline ComplexMatrix::ComplexMatrix(int nrows, int ncols)
{
	this->init(nrows, ncols);
	++mCount;
}

// zero out the ComplexMatrix
inline void ComplexMatrix::zero_out(void)
{
//	for (int i=0; i<size; i++) c[i] = complex<double>(0.0,0.0);
	memset ( c, 0, size*sizeof(complex<double>) );		// Peize Lin update 2015-12-08
}

/*
void need_more_memory()
{
	cout << "\n Sorry to crash... but the running need more momory! Exit." << endl;
	exit(0);
}
*/

inline void ComplexMatrix::init(const int nrows,const int ncols)
{
	this->nr = nrows;
	this->nc = ncols;
	this->size = nr * nc;
//	std::set_new_handler( need_more_memory );
	this->c = new complex<double>[size];
	assert(c != 0);
	this->zero_out();
}//mohan the principle:as simple as possible:modify inline 2007-10-13

// Copy constructor
inline ComplexMatrix::ComplexMatrix(const ComplexMatrix &m1)
{
	this->init(m1.nr, m1.nc);
	++mCount;
	for(int i=0; i<this->size; i++) c[i] = m1.c[i];
}

// deconstructor
inline ComplexMatrix::~ComplexMatrix()
{
	this->freemem();
}

// Free up memory for ComplexMatrix
inline void ComplexMatrix::freemem(void)
{
	delete[] c;
	c = NULL;
}

// reallocate memory for Complex Matrix
inline void ComplexMatrix::create(const int nrow,const int ncol)
{
	// because c has been 'new' in  init function.
	delete[] c;
	this->init(nrow, ncol);
	return;
}

inline void ComplexMatrix::set_as_identity_matrix()
{
/*	for(int i=0; i<nr; i++)
	{
		for(int j=0; j<nc; j++)
		{
			if(i==j) c[nc * i + j] = complex<double>(1.0, 0.0);  
			else c[nc * i + j] = complex<double>(0.0, 0.0); 
		}
	}
*/
	// Peize Lin update 2015-12-08
	zero_out();
	const int nrc_min = min(nr,nc);
	for( int i=0; i<nrc_min; ++i)
	{
		c[nc*i+i] = complex<double>(1.0,0.0);
	}
	return;
}

// Adding matrices, as a friend 
inline ComplexMatrix operator+(const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	ComplexMatrix tm(m1);
	tm+=m2;
	return tm;
}

// Subtracting matrices, as a friend 
inline ComplexMatrix operator-(const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	ComplexMatrix tm(m1);
	tm-=m2;
	return tm;
}

// Multiplying matrices, as a friend
// mprod = m1 * m2
inline ComplexMatrix operator*(const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	assert(m1.nc == m2.nr);	
	ComplexMatrix mprod(m1.nr, m2.nc);

/*							
	complex<double> z;
	for (int i = 0;i < m1.nr;i++)
	{
		for (int j = 0;j < m2.nc;j++)
		{
			z = complex<double>(0,0);
			for (int k = 0;k < m1.nc;k++)
			{
				z += m1(i, k) * m2(k, j);
			}
			mprod(i, j) = z;
		}
	}
*/
	// Peize Lin update 2015-12-07
	LapackConnector::zgemm('N', 'N',
		m1.nr, m2.nc, m1.nc,
		1.0, m1.c, m1.nc, m2.c, m2.nc,
		0.0, mprod.c, mprod.nc);
	return mprod;
}

// Peize Lin add 2015-12-08
// m1^T=m1('T') or m1^H=m1('H')
// mprod=m1*m2('L'eft) or mprod=m2*m1('R'ight)
// m1 stored in 'U'pper or 'L'ower
inline ComplexMatrix multiply_special(const char Symmetry, const char Side, const char Uplo, const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	assert(m1.nr == m1.nc);
	assert(m1.nc == m2.nr);
	ComplexMatrix mprod(m1.nr, m2.nc);
	if('T'==Symmetry)
	{
		LapackConnector::zsymm(Side, Uplo,
			  m1.nr, m2.nc,
			  1.0, m1.c, m1.nc, m2.c, m2.nc,
			  0.0, mprod.c, mprod.nc);
	}
	else if('H'==Symmetry)
	{
		LapackConnector::zhemm(Side, Uplo,
			  m1.nr, m2.nc,
			  1.0, m1.c, m1.nc, m2.c, m2.nc,
			  0.0, mprod.c, mprod.nc);
	}
	else
	{
		throw invalid_argument("Symmetry must be 'T' or 'H'");
	}
	return mprod;	
}

// Scale a ComplexMatrix
inline ComplexMatrix operator*(const complex<double> &c,const ComplexMatrix &m)
{
	ComplexMatrix sm(m);
	for (int i=0 ;i<m.size; i++) sm.c[i] *= c;
	return sm;
}

// ComplexMatrix scalar
inline ComplexMatrix operator*(const ComplexMatrix &m,const complex<double> &c)
{
	ComplexMatrix sm(m);
	for (int i = 0;i < m.size;i++) sm.c[i] *= c;
	return sm;
}

inline ComplexMatrix operator*(const double &r,const ComplexMatrix &m)
{
	ComplexMatrix sm(m);
	for(int i=0; i<m.size; i++) sm.c[i]*= r;
	return sm;
}

inline ComplexMatrix operator*(const ComplexMatrix &m,const double &r)
{
	ComplexMatrix sm(m);
	for (int i=0; i<m.size; i++) sm.c[i] *= r;
	return sm;
}

inline ComplexMatrix &ComplexMatrix::add(const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	for(int i=0; i<size; ++i) this->c[i] = m1.c[i] + m2.c[i];
	return *this;
}

inline ComplexMatrix &ComplexMatrix::minus(const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	for(int i=0; i<size; ++i) this->c[i] = m1.c[i] - m2.c[i];
	return *this;
}
			   
inline ComplexMatrix &ComplexMatrix::multiply(const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	LapackConnector::zgemm('N', 'N',
		m1.nr, m2.nc, m1.nc,
		1.0, m1.c, m1.nc, m2.c, m2.nc,
		0.0, this->c, this->nc);
	return *this;
}
		
// Peize Lin add 2015-12-08
// m1^T=m1('T') or m1^H=m1('H')
// mprod=m1*m2('L'eft) or mprod=m2*m1('R'ight)
// m1 stored in 'U'pper or 'L'ower		
inline ComplexMatrix &ComplexMatrix::multiply_special(const char Symmetry, const char Side, const char Uplo, const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	if('T'==Symmetry)
	{
		LapackConnector::zsymm(Side, Uplo,
			  m1.nr, m2.nc,
			  1.0, m1.c, m1.nc, m2.c, m2.nc,
			  0.0, this->c, this->nc);
	}
	else if('H'==Symmetry)
	{
		LapackConnector::zhemm(Side, Uplo,
			  m1.nr, m2.nc,
			  1.0, m1.c, m1.nc, m2.c, m2.nc,
			  0.0, this->c, this->nc);
	}
	else
	{
		throw invalid_argument("Symmetry must be 'T' or 'H'");
	}
	return *this;	
}

inline ComplexMatrix& ComplexMatrix::operator=(const ComplexMatrix &m)
{
	if(m.nr!=this->nr || m.nc!=this->nc)
	{
		cout << "\n row/col number can't match in ComplexMatrix '=' operator\n";
		cout << " this nr = " << this->nr;
		cout << " this nc = " << this->nc;
		cout << " in nr = " << m.nr;
		cout << " in nc = " << m.nc;
		exit(0);
	}
/*	else {
		this->create(m.nr, m.nc);
	} */	// why ?

//	for(int i=0;i<this->size;i++) c[i] = m.c[i];
	memcpy ( c, m.c, this->size * sizeof(complex<double>) );	// Peize Lin update 2015-12-08
	return *this;
}

inline ComplexMatrix& ComplexMatrix::operator*=(const complex<double> &s)
{
	for (int i = 0;i < this->size;i++) c[i] *= s;
	return *this;
}

// Accumulate to a ComplexMatrix in place
inline ComplexMatrix& ComplexMatrix::operator+=(const ComplexMatrix &m)
{
	for(int i=0; i<size; i++) this->c[i] += m.c[i];
	return *this;
}

// decumulate to a ComplexMatrix in place
inline ComplexMatrix& ComplexMatrix::operator-=(const ComplexMatrix &m)
{
	for(int i=0; i<size; i++) this->c[i] -= m.c[i];
	return *this;
}

// Returns trace of ComplexMatrix
inline complex<double> trace(const ComplexMatrix &m)
{
	complex<double> tr=complex<double>(0,0);
	assert(m.nr == m.nc);
	for (int i=0; i<m.nr; i++) tr += m(i, i);
	return tr;
}

// Do mout += s*min
inline void scale_accumulate(const complex<double> &s,
                      const ComplexMatrix &min,
                      ComplexMatrix &mout)
{
	assert(min.nr == mout.nr);
	assert(min.nc == mout.nc);
	for (int j=0; j<min.size; j++)
	{
		mout.c[j] += s * min.c[j];
	}
	return;
}

// Do mout[i] += s*min[i]
inline void scale_accumulate(const int &nmat,
                      const complex<double> &s,
                      ComplexMatrix **min,
                      ComplexMatrix **mout)
{
	assert(nmat>=0);
	for (int i=0; i<nmat; i++)
	{
		scale_accumulate(s, *min[i], *mout[i]);
	}
	return;
}

// Do mout = s1*m1 + s2*m2
inline void scaled_sum(const complex<double> &s1,
                const ComplexMatrix &m1,
                const complex<double> &s2,
                const ComplexMatrix &m2,
                ComplexMatrix &mout)
{
	assert(m1.nr == m2.nr);
	assert(m1.nr == mout.nr);
	assert(m1.nc == m2.nc);
	assert(m1.nc == mout.nc);

	for(int i=0; i<m1.size; i++)
	{
		mout.c[i] = s1 * m1.c[i] + s2 * m2.c[i];
	}
	return;
}

// Does mout[i] = s1*m1[i] + s2*m2[i]
inline void scaled_sum(const int &nmat,
                const complex<double> &s1,
                ComplexMatrix **m1,
                const complex<double> &s2,
                ComplexMatrix **m2,
                ComplexMatrix **mout)
{
	assert(nmat>0);
	for(int i=0; i<nmat; i++)
	{
		scaled_sum(s1, *m1[i], s2, *m2[i], *mout[i]);
	}
	return;
}

inline double abs2_row(const ComplexMatrix &m,const int ir)
{
	double r=0.0;
	complex<double> z;
	for(int ic=0;ic<m.nc;ic++)
	{
		z = m.c[ m.nc*ir + ic];
		r += z.real()*z.real() + z.imag()*z.imag();
	}
	return r;
}

inline double abs2_column(const ComplexMatrix &m,const int ic)
{
	double r=0.0;
	complex<double> z;
	for(int ir=0;ir<m.nr;ir++)
	{
		z = m.c[ m.nc*ir + ic ];
		r += z.real()*z.real() + z.imag()*z.imag();
	}
	return r;
}

// returns absolute square magnitude of sum of all ComplexMatrix elements
inline double abs2(const ComplexMatrix &m)
{
	double r=0.0;
	complex<double> z;
	
	for (int i = 0;i < m.size;i++)
	{
		z = m.c[i];
		r += z.real() * z.real() + z.imag() * z.imag();
	}
	return r;
}

// Same for an array of matrices
inline double abs2(const int nmat, ComplexMatrix **m)
{
	double r = 0.0;
	for (int i = 0;i < nmat;i++)
	{
		r += abs2(*m[i]);
	}
	return r;
}

inline double abs2_triangle_matrix( const char Uplo, ComplexMatrix &m)
{
	double r(0.0);
	complex<double> z;
	
	if('U'==Uplo)
	{
		for (int i=0; i<m.nr; ++i)
		{
			for (int j=i+1; j<m.nc; ++j)
			{
				z = m(i,j);
				r += z.real() * z.real() + z.imag() * z.imag();
			}
		}		
	}
	else if('L'==Uplo)
	{
		for (int i=0; i<m.nr; ++i)
		{
			for (int j=0; j<i; ++j)
			{
				z = m(i,j);
				r += z.real() * z.real() + z.imag() * z.imag();				
			}
		}
	}
	else
	{
		throw invalid_argument("Uplo must be 'U' or 'L'");
	}
}

inline ComplexMatrix transpose(const ComplexMatrix &m, const bool &conjugate)
{
	ComplexMatrix tm(m.nc, m.nr);
	if(conjugate)
	{
		for (int i = 0;i < m.nr;i++)
		{
			for (int j = 0;j < m.nc;j++)
			{
				tm(j, i) = conj ( m(i, j) );
			}
		}
	}
	else
	{
		for (int i = 0;i < m.nr;i++)
		{
			for (int j = 0;j < m.nc;j++)
			{
				tm(j, i) = m(i, j);
			}
		}
	}
	return tm;
}

#endif