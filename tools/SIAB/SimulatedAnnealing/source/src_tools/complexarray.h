/*******************************************
 * ESCP:Electro-Structure Calculate Package.
 ********************************************/
#ifndef COMPLEXARRAY_H
#define COMPLEXARRAY_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <complex>

using namespace std;

class ComplexArray
{
public:
	complex<double> *ptr;

	ComplexArray(const int d1 = 1 ,const int d2 = 1,const int d3 = 1);
	ComplexArray(const int d1, const int d2,const int d3,const int d4);
	~ComplexArray();

	void create(const int d1,const int d2,const int d3);
	void create(const int d1,const int d2,const int d3,const int d4);

	const ComplexArray &operator=(const ComplexArray &right);
	const ComplexArray &operator=(const complex<double> &right);

	complex<double> &operator()(const int d1,const int d2,const int d3);
	complex<double> &operator()(const int d1,const int d2,const int d3,const int d4);

	const complex<double> &operator()(const int d1,const int d2,const int d3)const;
	const complex<double> &operator()(const int d1,const int d2,const int d3,const int d4)const;

	void zero_out(void);

	const int getSize() const
	{ return size;}

	const int getDim() const
	{ return dim;}

	const int getBound1() const
	{ return bound1;}

	const int getBound2() const
	{ return bound2;}

	const int getBound3() const
	{ return bound3;}

	const int getBound4() const
	{ return bound4;}

	static const int getArrayCount(void)
	{ return arrayCount;}

private:
	int size;
	int dim;
	int bound1, bound2, bound3, bound4;
	static int arrayCount;

	void freemem();
};

#endif	// ComplexArray class
