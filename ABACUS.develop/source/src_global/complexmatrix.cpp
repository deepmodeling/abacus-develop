//==========================================================
// AUTHOR : Lixin He, Mohan Chen
// LAST UPDATE : 2009-03-23 modify "=" operator
//==========================================================

#include <cassert>
#include <new>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "complexmatrix.h"
#include "lapack_connector.h"

// constructor with sizes
ComplexMatrix::ComplexMatrix(const int nrows, const int ncols, const bool flag_zero)
	:nr(nrows),
	 nc(ncols),
	 size(nrows*ncols),
	 c(nullptr)
{
	if( size )
	{
		c = new complex<double>[size];
		if(flag_zero)	zero_out();
	}
}

// zero out the ComplexMatrix
void ComplexMatrix::zero_out(void)
{
	for (int i=0; i<size; i++) c[i] = complex<double>(0.0,0.0);
}

/*
void need_more_memory()
{
	cout << "\n Sorry to crash... but the running need more momory! Exit." << endl;
	exit(0);
}
*/

// Copy constructor
ComplexMatrix::ComplexMatrix(const ComplexMatrix &m1)
	:nr(m1.nr),
	 nc(m1.nc),
	 size(m1.size),
	 c(nullptr)
{
	if(size)
	{
		c = new complex<double>[size];
		memcpy( c, m1.c, size*sizeof(complex<double>) );
	}
}

// Peize Lin add 2016-08-05
ComplexMatrix::ComplexMatrix( ComplexMatrix && m1 )
	:nr(m1.nr),
	 nc(m1.nc),
	 size(m1.size),
	 c(m1.c)
{
	m1.nr = m1.nc = m1.size = 0;
	m1.c = nullptr;
}

// Peize Lin add 2017-03-29
ComplexMatrix::ComplexMatrix(const matrix &m)
	:nr(m.nr),
	 nc(m.nc),
	 size(m.nr*m.nc),
	 c(nullptr)
{
	if( size )
	{
		c = new complex<double>[size];
		for( int i=0; i<size; ++i)
			c[i] = m.c[i];
	}
}

// deconstructor
ComplexMatrix::~ComplexMatrix()
{
	if(c)
	{
		delete[] c;
		c = nullptr;
	}
}

// reallocate memory for Complex Matrix
void ComplexMatrix::create(const int nr_in, const int nc_in, const bool flag_zero)
{
	if( nr_in && nc_in )
	{
		if(c)
		{
			const int size_in=nr_in*nc_in;
			if( size_in!=nr*nc )
			{
				delete[] c;
				c = new complex<double>[size_in];
			}
		}
		else
		{
			c = new complex<double>[nr_in * nc_in];
		}

		nr = nr_in;
		nc = nc_in;
		size = nr*nc;
		if(flag_zero)	zero_out();
	}
	else
	{
		if(c)	delete[] c;
		c = nullptr;
		nr = nr_in;
		nc = nc_in;
		size = nr*nc;
	}
}

void ComplexMatrix::set_as_identity_matrix()
{
	for(int i=0; i<nr; i++)
	{
		for(int j=0; j<nc; j++)
		{
			if(i==j) c[nc * i + j] = complex<double>(1.0, 0.0);
			else c[nc * i + j] = complex<double>(0.0, 0.0);
		}
	}
	return;
}

// Adding matrices, as a friend
ComplexMatrix operator+(const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	assert(m1.nr == m2.nr);
	assert(m2.nc == m2.nc);
	
	ComplexMatrix tm(m1);
	tm+=m2;
	return tm;
}

// Subtracting matrices, as a friend
ComplexMatrix operator-(const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	assert(m1.nr == m2.nr);
	assert(m2.nc == m2.nc);
	
	ComplexMatrix tm(m1);
	tm-=m2;
	return tm;
}

// Multiplying matrices, as a friend
// mprod = m1 * m2
ComplexMatrix operator*(const ComplexMatrix &m1, const ComplexMatrix &m2)
{
	assert(m1.nc == m2.nr);
	ComplexMatrix mprod(m1.nr, m2.nc);

	complex<double> z;
//	for (int i = 0;i < m1.nr;i++)
//	{
//		for (int j = 0;j < m2.nc;j++)
//		{
//			z = complex<double>(0,0);
//			for (int k = 0;k < m1.nc;k++)
//			{
//				z += m1(i, k) * m2(k, j);
//			}
//			mprod(i, j) = z;
//		}
//	}
	// Peize Lin accelerate 2017-10-27
	LapackConnector::gemm('N', 'N', m1.nr, m2.nc, m1.nc,
		1, m1.c, m1.nc, m2.c, m2.nc,
		0, mprod.c, mprod.nc);
	return mprod;
}

// Scale a ComplexMatrix
ComplexMatrix operator*(const complex<double> &c,const ComplexMatrix &m)
{
	ComplexMatrix sm(m);
	for (int i=0 ;i<m.size; i++) sm.c[i] *= c;
	return sm;
}

// ComplexMatrix scalar
ComplexMatrix operator*(const ComplexMatrix &m,const complex<double> &c)
{
	ComplexMatrix sm(m);
	for (int i = 0;i < m.size;i++) sm.c[i] *= c;
	return sm;
}

ComplexMatrix operator*(const double &r,const ComplexMatrix &m)
{
	ComplexMatrix sm(m);
	for(int i=0; i<m.size; i++) sm.c[i]*= r;
	return sm;
}

ComplexMatrix operator*(const ComplexMatrix &m,const double &r)
{
	ComplexMatrix sm(m);
	for (int i=0; i<m.size; i++) sm.c[i] *= r;
	return sm;
}

ComplexMatrix& ComplexMatrix::operator=(const ComplexMatrix &m)
{
	this->create(m.nr, m.nc, false);
	memcpy( c, m.c, size*sizeof(complex<double>) );
	return *this;
}

// Peize Lin add 2016-08-05
ComplexMatrix& ComplexMatrix::operator=( ComplexMatrix && m )
{
	nr = m.nr;		nc = m.nc;		size = m.size;
	if(c)	delete[] c;		
	c  = m.c;
	m.nr = m.nc = m.size = 0;
	m.c = nullptr;
	return *this;
}

ComplexMatrix& ComplexMatrix::operator*=(const complex<double> &s)
{
	for (int i = 0;i < this->size;i++) c[i] *= s;
	return *this;
}

// Accumulate to a ComplexMatrix in place
ComplexMatrix& ComplexMatrix::operator+=(const ComplexMatrix &m)
{
	for(int i=0; i<size; i++) this->c[i] += m.c[i];
	return *this;
}

// decumulate to a ComplexMatrix in place
ComplexMatrix& ComplexMatrix::operator-=(const ComplexMatrix &m)
{
	for(int i=0; i<size; i++) this->c[i] -= m.c[i];
	return *this;
}

// Peize Lin add 2017-03-29
matrix ComplexMatrix::real() const
{
	matrix m(nr,nc,false);
	for( int i=0; i<this->size; ++i) m.c[i] = c[i].real();
	return m;
}

// Returns trace of ComplexMatrix
complex<double> trace(const ComplexMatrix &m)
{
	complex<double> tr=complex<double>(0,0);
	assert(m.nr == m.nc);
	for (int i=0; i<m.nr; i++) tr += m(i, i);
	return tr;
}

// Do mout += s*min
void scale_accumulate(const complex<double> &s,
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
void scale_accumulate(const int &nmat,
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
void scaled_sum(const complex<double> &s1,
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
void scaled_sum(const int &nmat,
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


double abs2_row(const ComplexMatrix &m,const int ir)
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

double abs2_column(const ComplexMatrix &m,const int ic)
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
double abs2(const ComplexMatrix &m)
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
double abs2(const int nmat, ComplexMatrix **m)
{
	double r = 0.0;
	for (int i = 0;i < nmat;i++)
	{
		r += abs2(*m[i]);
	}
	return r;
}

ComplexMatrix transpose(const ComplexMatrix &m, const bool &conjugate)
{
	ComplexMatrix tm(m.nc, m.nr, false);
	if(conjugate)
		for (int i = 0;i < m.nr;i++)
			for (int j = 0;j < m.nc;j++)
				tm(j, i) = conj ( m(i, j) );
	else
		for (int i = 0;i < m.nr;i++)
			for (int j = 0;j < m.nc;j++)
				tm(j, i) = m(i, j);
	return tm;
}

ComplexMatrix conj(const ComplexMatrix &m)
{
	ComplexMatrix cm( m.nr, m.nc, false );
	for(int i=0; i!=m.size; ++i)
		cm.c[i] = conj(m.c[i]);
	return cm;
}