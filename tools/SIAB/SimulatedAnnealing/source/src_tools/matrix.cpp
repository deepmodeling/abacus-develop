/* matrix.cpp file */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>

#include <cmath>
#include <cstdlib>

using namespace std;
#include "matrix.h"

//*********************************************************
// The init() function is the main initialization routine.
// Sets up sizes and allocates memory for matrix class.
// All constructors call init()
// ********************************************************

int matrix::mCount = 0;

void matrixAlloc()
{
	cout << "\n Allocation error for Matrix " << endl;
	abort();
}

void matrix::init(const int nrows,const int ncols)
{
	nr = nrows;
	nc = ncols;
	c = NULL;
//	hermetian = 0;
	set_new_handler(matrixAlloc);

	if(nr*nc == 0) c = NULL;
	else
	{
		c = new double[nrows*ncols]();
		this->zero_out();
		assert(c!=0);// terminate if memory not allocated
	}

	mCount++;
}

// Free up memory for matrix
void matrix::freemem(void)
{
	delete [] c;
	c = NULL;
}

// constructor with sizes
matrix::matrix(const int nrows,const int ncols)
{
	this->init(nrows, ncols);
}

//******************
//
// Copy constructor
//
//******************
matrix::matrix(const matrix &m1)
{
	init(m1.nr, m1.nc);
	//hermetian = m1.hermetian;
	// Copy over m1 contents
	for (int i = 0; i < nr*nc; i++) c[i] = m1.c[i];
}

//*************
//
// destructor
//
//*************
matrix::~matrix()
{
	this->freemem();
}

//******************************
// reallocate memory for matrix
//******************************
void matrix::create(const int nrow,const int ncol)
{
	delete [] c;
	nr = nrow;
	nc = ncol;
	c = new double[nrow * ncol]();
	return;
}

/* Assignment:  nonstandard in that it returns void.  To make it standard,
 * replace void -> matrix and uncomment the return *this; */

void matrix::operator=(const matrix & m1)
{
	for (int i = 0; i < nr*nc; i++) c[i] = m1.c[i];
}

/* Adding matrices, as a friend */
matrix operator+(const matrix &m1, const matrix &m2)
{
	assert(m1.nr == m2.nr);
	assert(m2.nc == m2.nc);

	matrix tm(m1);
	for (int i = 0; i < m1.nr*m1.nc; i++) tm.c[i] += m2.c[i];
	return tm;
}

/* Subtracting matrices, as a friend */
matrix operator-(const matrix &m1, const matrix &m2)
{
	assert(m1.nr == m2.nr);
	assert(m2.nc == m2.nc);

	matrix tm(m1);
	for(int i = 0; i < m1.nr*m1.nc; i++) tm.c[i] -= m2.c[i];
	return tm;
}

//***************************************
//
// Multiplying matrices, as a friend 
//
// *************************************
matrix operator*(const matrix &m1, const matrix &m2)
{
	assert(m1.nr == m2.nr);
	assert(m2.nc == m2.nc);
	
	// allocate the result and zero it out
	matrix mprod(m2.nr, m1.nc);

	mprod.zero_out();

	// do the multiply and return
	for (int i = 0;i < m2.nr;i++)
	{
		for (int j = 0;j < m1.nc;j++)
		{
			for (int k = 0;k < m2.nc;k++)
			{
				mprod(i, j) += m2(i, k) * m1(k, j);
			}
		}
	}

	return mprod;
}

/* Scale a matrix */
matrix operator*(const double &s, const matrix &m)
{
	matrix sm(m);
	for (int i = 0; i < m.nr*m.nc; i++) sm.c[i] *= s;
	return sm;
}

/* matrix * double */
matrix operator*(const matrix &m,const double &s)
{
	matrix sm(m);
	for (int i = 0; i < m.nr*m.nc; i++)sm.c[i] *= s;
	return sm;
}

/* Scale a matrix in place */
void matrix::operator*=(const double &s)
{
	for (int i = 0; i < nc*nr; i++) c[i] *= s;
}

/* Accumulate to a matrix in place */
void matrix::operator+=(const matrix & m)
{
	if (nr != m.nr || nc != m.nc)
	{
		cout << " void matrix::operator+=(const matrix &m) has size mismatch\n";
		abort();
	}
	for (int i = 0; i < nc*nr; i++) c[i] += m.c[i];
}


/* decumulate to a matrix in place */
void matrix::operator-=(const matrix & m)
{
	if (nr != m.nr || nc != m.nc)
	{
		cout << "void matrix::operator-=(const matrix &m) has size mismatch\n";
		abort();
	}
	for (int i = 0; i < nc*nr; i++) c[i] -= m.c[i];
}

/* zero out the matrix */
void matrix::zero_out(void)
{
	for(int i = 0; i < nr*nc; i++)
	{
		c[i] = 0.0;
	}
	return;
}


matrix transpose(matrix m)
{
	matrix tm(m.nc, m.nr);

	for (int i = 0;i < m.nr;i++)
	{
		for (int j = 0;j < m.nc;j++)
		{
			tm(j, i) = m(i, j);
		}
	}
	return tm;
}
