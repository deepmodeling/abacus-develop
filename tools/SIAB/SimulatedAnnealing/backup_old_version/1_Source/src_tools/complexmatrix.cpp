//==========================================================
// AUTHOR : Lixin He, Mohan Chen
// LAST UPDATE : 2009-03-23 modify "=" operator
//==========================================================

#include <cassert>
#include <new>
#include <cstdlib>
#include <iostream>
#include "complexmatrix.h"

int ComplexMatrix::mCount = 0;

// constructor with sizes
ComplexMatrix::ComplexMatrix(int nrows, int ncols)
{
	this->init(nrows, ncols);
	++mCount;
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

void ComplexMatrix::init(const int nrows,const int ncols)
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
ComplexMatrix::ComplexMatrix(const ComplexMatrix &m1)
{
	this->init(m1.nr, m1.nc);
	++mCount;
	for(int i=0; i<this->size; i++) c[i] = m1.c[i];
}

// deconstructor
ComplexMatrix::~ComplexMatrix()
{
	this->freemem();
}

// Free up memory for ComplexMatrix
void ComplexMatrix::freemem(void)
{
	delete[] c;
	c = NULL;
}

// reallocate memory for Complex Matrix
void ComplexMatrix::create(const int nrow,const int ncol)
{
	// because c has been 'new' in  init function.
	delete[] c;
	this->init(nrow, ncol);
	return;
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
	ComplexMatrix tm(m1);
	tm+=m2;
	return tm;
}

// Subtracting matrices, as a friend 
ComplexMatrix operator-(const ComplexMatrix &m1, const ComplexMatrix &m2)
{
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
	if(m.nr!=this->nr || m.nc!=this->nc)
	{
		cout << "\n row/col number can't match in ComplexMatrix '=' operator\n";
		cout << " this nr = " << this->nr;
		cout << " this nc = " << this->nc;
		cout << " in nr = " << m.nr;
		cout << " in nc = " << m.nc;
		exit(0);
	}
	else {
		this->create(m.nr, m.nc);
	}
	for(int i=0;i<this->size;i++) c[i] = m.c[i];
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
