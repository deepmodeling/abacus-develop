//==========================================================
// Author : Lixin He, Mohan Chen
// Last Update : 2009-3-8
//==========================================================
#ifndef COMPLEXMATRIX_H
#define COMPLEXMATRIX_H

#include <complex>
using namespace std;

class ComplexMatrix
{

public:

	complex<double> *c;
	int nr;
	int nc;
	int size;

	//==============
	// Constructors
	//==============
	ComplexMatrix(const int nrows = 1,const int ncols = 1);

	//===================
	// Copy constructor
	//==================
	ComplexMatrix(const ComplexMatrix &m1);
	~ComplexMatrix();

	//============
	// Operators
	//============
	complex<double> &operator()(const int i,const int j)
	{
		return c[nc * i + j];//mohan modify in-line 2007-10-1
	}
	const complex<double> &operator()(const int i,const int j)const
	{
		return c[nc * i + j];//mohan modify in-line 2007-10-13
	}

	friend ComplexMatrix operator+(const ComplexMatrix &m1,  const ComplexMatrix &m2);
	friend ComplexMatrix operator-(const ComplexMatrix &m1,  const ComplexMatrix &m2);
	friend ComplexMatrix operator*(const ComplexMatrix &m1,  const ComplexMatrix &m2);
	friend ComplexMatrix operator*(const complex<double> &s, const ComplexMatrix &m);
	friend ComplexMatrix operator*(const ComplexMatrix &m,   const complex<double> &s);
	friend ComplexMatrix operator*(const double &s,          const ComplexMatrix &m);
	friend ComplexMatrix operator*(const ComplexMatrix &m,   const double &s);
	ComplexMatrix& operator=(const ComplexMatrix &m);
	ComplexMatrix& operator*=(const complex<double> &s);
	ComplexMatrix& operator+=(const ComplexMatrix &m);
	ComplexMatrix& operator-=(const ComplexMatrix &m);

	//==================
	// member function:
	//==================
	void create(const int nrow,const int ncol);
	void zero_out(void);
	static int& getMCount(void){return mCount;}
	void set_as_identity_matrix(void);
private:

	static int mCount;
	void freemem(void);//mohan add 2007-11-20
	void init(const int nrows,const int ncols);
};

complex<double> trace(const ComplexMatrix &m);
// mohan add 2008-7-1
double abs2_row(const ComplexMatrix &m,const int ir);
// mohan add 2008-7-1
double abs2_column(const ComplexMatrix &m,const int ic);
double abs2(const ComplexMatrix &m);
double abs2(const int nmat, ComplexMatrix **m);

ComplexMatrix transpose(const ComplexMatrix &m, const bool &conjugate);

void scale_accumulate(
		const complex<double> &s, 
		const ComplexMatrix &min, 
		ComplexMatrix &mout);

void scale_accumulate(
		const int &nmat, 
		const complex<double> &s, 
		ComplexMatrix **min, 
		ComplexMatrix **mout);

void scaled_sum(
		const complex<double> &s1, 
		const ComplexMatrix &m1, 
		const complex<double> &s2, 
		const ComplexMatrix &m2, 
		ComplexMatrix &mout);

void scaled_sum(
		const int &nmat,
		const complex<double> &s1, 
		ComplexMatrix **m1, 
		const complex<double> &s2,
		ComplexMatrix **m2, 
		ComplexMatrix **mout);
#endif
