/*******************************************
 * ESCP:Electro-Structure Calculate Package.
 ********************************************/

#ifndef REALARRAY_H
#define REALARRAY_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

using namespace std;

class realArray
{
public:
	double *ptr;

	realArray(const int d1 = 1 ,const int d2 = 1,const int d3 = 1);
	realArray(const int d1, const int d2,const int d3,const int d4);
	~realArray();

	void create(const int d1,const int d2,const int d3);
	void create(const int d1,const int d2,const int d3,const int d4);

	const realArray &operator=(const realArray &right);
	const realArray &operator=(const double &right);

	double &operator()(const int d1,const int d2,const int d3);
	double &operator()(const int d1,const int d2,const int d3,const int d4);

	const double &operator()(const int d1,const int d2,const int d3)const;
	const double &operator()(const int d1,const int d2,const int d3,const int d4)const;

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

//**************************************************
// set elements of a as zeros which a is 1_d array.
//**************************************************
template<class T>
void zeros(T *u,const int n)
{
	assert(n>0);
	for (int i = 0;i < n;i++) u[i] = 0;
}

#endif	// realArray class
