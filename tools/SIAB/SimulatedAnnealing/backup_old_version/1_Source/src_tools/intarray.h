/*******************************************
 * ESCP:Electro-Structure Calculate Package.
 ********************************************/

#ifndef INTARRAY_H
#define INTARRAY_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

using namespace std;

class IntArray
{
public:
	int *ptr;

	// Constructors for different dimesnions
	IntArray(const int d1 = 1, const int d2 = 1);
	IntArray(const int d1, const int d2,const int d3);
	IntArray(const int d1, const int d2,const int d3,const int d4);
	IntArray(const int d1, const int d2,const int d3,const int d4,const int d5);
	IntArray(const int d1, const int d2,const int d3,const int d4,const int d5,const int d6);

	~IntArray();

	void create(const int d1, const int d2);
	void create(const int d1, const int d2, const int d3);
	void create(const int d1, const int d2, const int d3, const int d4);
	void create(const int d1, const int d2, const int d3, const int d4, const int d5);
	void create(const int d1, const int d2, const int d3, const int d4, const int d5, const int d6);

	const IntArray &operator=(const IntArray &right);
	const IntArray &operator=(const int &right);

	int &operator()(const int d1, const int d2);
	int &operator()(const int d1, const int d2, const int d3);
	int &operator()(const int d1, const int d2, const int d3,const int d4);
	int &operator()(const int d1, const int d2, const int d3, const int d4, const int d5);
	int &operator()(const int d1, const int d2, const int d3, const int d4, const int d5, const int d6);

	const int &operator()(const int d1,const int d2)const;
	const int &operator()(const int d1,const int d2,const int d3)const;
	const int &operator()(const int d1,const int d2,const int d3,const int d4)const;
	const int &operator()(const int d1,const int d2,const int d3,const int d4, const int d5)const;
	const int &operator()(const int d1,const int d2,const int d3,const int d4, const int d5, const int d6)const;

	void zero_out(void);

	const int getSize() const{ return size;}
	const int getDim() const{ return dim;}
	const int getBound1() const{ return bound1;}
	const int getBound2() const{ return bound2;}
	const int getBound3() const{ return bound3;}
	const int getBound4() const { return bound4;}
	const int getBound5() const { return bound5;}
	const int getBound6() const { return bound6;}

	static const int getArrayCount(void)
	{ return arrayCount;}

private:
	int size;
	int dim;
	int bound1, bound2, bound3, bound4, bound5, bound6;
	static int arrayCount;
	void freemem();
};

#endif	// IntArray class
