/*
    DFT++ is a density functional package developed by the research group
    of Professor Tomas Arias

    Copyright 1996-2003 Sohrab Ismail-Beigi
*/

#ifndef COMPLEX_ARRAY_H
#define COMPLEX_ARRAY_H

#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

class ComplexArray
{
public:
	int ndata;  // number of data elements;
//	const static int datasize;	// = sizeof(complex);
	complex<double> *ptr; // data array
	int bound1, bound2, bound3, bound4;

	ComplexArray(int bnd1 = 1);
	ComplexArray(int bnd1, int bnd2, int bnd3);
	ComplexArray(int bnd1, int bnd2, int bnd3, int bnd4);

	~ComplexArray();

	void init(int);
	void freemem();

	void create(int bnd1, int bnd2, int bnd3);
	void create(int bnd1, int bnd2, int bnd3, int bnd4);

	// get and release scratch space
//  void get_temp(int);
//  void release_temp();

	void operator=(const ComplexArray &cd);
	inline void operator=(complex <double> c);
//  inline std::complex < double>  &operator()(int i, int j, int k)const
//{return d[(i * bound2 + j) * bound3 +k];}
//  inline void operator=(complex < double> c);

	ComplexArray operator+(const ComplexArray &cd);
	void operator+=(const ComplexArray &cd);
	ComplexArray operator-(const ComplexArray &cd);
	void operator-=(const ComplexArray &cd);
	ComplexArray operator*(const double r);
	ComplexArray operator*(const std::complex < double> c);
	void operator*=(const double r);
	void operator*=(const std::complex < double> c);

	void operator*=(const ComplexArray &in);

	// subscript operator
	complex < double> &operator()(int, int, int);
	complex < double> &operator()(int, int, int, int);
//  complex < double> &operator()(int, int, int, int, int);

	const complex < double> &operator()(int, int, int)const;
	const complex < double> &operator()(int, int, int, int)const;
//  const complex < double> &operator()(int, int, int, int, int)const;

	void zero_out(void);
	void negate(void);
	void randomize(void);//uniform distribution

//  void write(char *fname);
//  void write(FILE *fp);
//  void writea(char *fname);
//  void read(char *fname);
	void print();
};

ComplexArray operator*(double r, const ComplexArray &cd);
ComplexArray operator*(std::complex < double> c, const ComplexArray &cd);

double abs2(const ComplexArray &cd);

void
add_scale_abs2(const std::complex < double> &c, const ComplexArray & in,
               ComplexArray &out);

std::complex < double> dot(const ComplexArray &cd1, const ComplexArray &cd2);
void scale_accumulate(double r, const ComplexArray &cd1, ComplexArray &cd2);
void scale_accumulate(std::complex < double> c, const ComplexArray &cd1, ComplexArray &cd2);
void scaled_sum(double r1, const ComplexArray &cd1,
                double r2, const ComplexArray &cd2,
                ComplexArray &cd3);

void scaled_sum(double r1, ComplexArray &cd1,
                double r2, ComplexArray &cd2,
                ComplexArray &cd3);

void point_mult(ComplexArray &a1, ComplexArray &in2, ComplexArray &out);

// set elements of u as zero which u is 1_d complex array
template <class T>
void zeros(complex <T> *u, int n)
{
	if (n == 0 || u == 0)
	{
		cout << "\n error in zeros(),n or u = 0";
		return;
	}

	for (int i = 0;i < n;i++)
	{
		u[i] = complex <T> (0.0, 0.0);
	}
}

#endif // COMPLEX_ARRAY_H
