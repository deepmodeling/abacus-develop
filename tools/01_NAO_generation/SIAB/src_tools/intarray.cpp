/*******************************************
 * ESCP:Electro-Structure Calculate Package.
 ********************************************/
#include <cstdlib>
#include "intarray.h"

int IntArray::arrayCount = 0;

void IntArrayAlloc()
{
	cout << "\n Allocation error for IntArray " << endl;
	exit(0);
}

IntArray::IntArray(const int d1,const int d2)
{
	dim = 2;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = bound4 = 0;
	size = bound1 * bound2;
	ptr = new int[size];
	assert( ptr != 0);
	++arrayCount;
}

IntArray::IntArray(const int d1,const int d2,const int d3)
{
	dim = 3;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = 0;
	set_new_handler(IntArrayAlloc);
	size = bound1 * bound2 * bound3 ;	//* sizeof(float);
	ptr = new int[size]();assert(ptr != 0);
	zero_out();
	++arrayCount;
}

IntArray::IntArray(const int d1,const int d2,const int d3,const int d4)
{
	dim = 4;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = (d4 <= 0) ? 1 : d4;
	set_new_handler(IntArrayAlloc);
	size = bound1 * bound2 * bound3 * bound4 ;	//* sizeof(float);
	ptr = new int[size]();assert(ptr != 0);
	zero_out();
	++arrayCount;
}

IntArray::IntArray(const int d1,const int d2,const int d3,
		const int d4,const int d5)
{
	dim = 5;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = (d4 <= 0) ? 1 : d4;
	bound5 = (d5 <= 0) ? 1 : d5;
	set_new_handler(IntArrayAlloc);
	size = bound1 * bound2 * bound3 * bound4 * bound5;
	ptr = new int[size]();assert(ptr != 0);
	zero_out();
	++arrayCount;
}

IntArray::IntArray(const int d1,const int d2,const int d3,
		const int d4,const int d5,const int d6)
{
	dim = 6;
	bound1 = (d1 <= 0) ? 1 : d1;
    bound2 = (d2 <= 0) ? 1 : d2;
    bound3 = (d3 <= 0) ? 1 : d3;
    bound4 = (d4 <= 0) ? 1 : d4;
    bound5 = (d5 <= 0) ? 1 : d5;
	bound6 = (d6 <= 0) ? 1 : d6;
    set_new_handler(IntArrayAlloc);
    size = bound1 * bound2 * bound3 * bound4 * bound5 * bound6;
	ptr = new int[size]();assert(ptr != 0);
	zero_out();
	++arrayCount;
}

//********************************
// Destructor for class IntArray
//********************************
IntArray ::~IntArray()
{
    freemem();
}

void IntArray::freemem()
{
	delete [] ptr;
	ptr = NULL;
}

void IntArray::create(const int d1,const int d2,const int d3,const int d4,const int d5,const int d6)
{
	size = d1 * d2 * d3 * d4 * d5 * d6;assert(size>0);
	dim = 6;
	bound1 = d1;bound2 = d2;bound3 = d3;bound4 = d4;bound5 = d5;bound6 = d6;
	delete[] ptr; ptr = new int[size];
	assert(ptr != 0);zero_out();
}

void IntArray::create(const int d1,const int d2,const int d3,const int d4,const int d5)
{
	size = d1 * d2 * d3 * d4 * d5;assert(size>0);
	dim = 5;
	bound1 = d1;bound2 = d2;bound3 = d3;bound4 = d4;bound5 = d5;
	delete[] ptr; ptr = new int[size];
	assert(ptr != 0);zero_out();
}

void IntArray::create(const int d1,const int d2,const int d3,const int d4)
{
	size = d1 * d2 * d3 * d4;assert(size>0);
	dim = 4;
	bound1 = d1;bound2 = d2;bound3 = d3;bound4 = d4;
	delete[] ptr; ptr = new int[size];
	assert(ptr != 0);zero_out();
}

void IntArray::create(const int d1,const int d2,const int d3)
{
	size = d1 * d2 * d3;assert(size>0);
	dim = 3;
	bound1 = d1;bound2 = d2;bound3 = d3;bound4 = 1;
	delete [] ptr;ptr = new int[size];
	assert(ptr != 0);zero_out();
}

void IntArray::create(const int d1, const int d2)
{
	size = d1 * d2;assert(size>0);
	dim = 2;
	bound1 = d1;bound2 = d2;bound3 = bound4 = 1;
	delete[] ptr;ptr = new int[size];
	assert(ptr !=0 );zero_out();
}

const IntArray &IntArray::operator=(const IntArray &right)
{
	for (int i = 0;i < size;i++) ptr[i] = right.ptr[i];
	return *this;// enables x = y = z;
}

const IntArray &IntArray::operator=(const int &value)
{
	for (int i = 0;i < size;i++) ptr[i] = value;
	return *this;// enables x = y = z;
}

//********************************************************
// overloaded subscript operator for const Int Array
// const reference return create an cvakue
//********************************************************
const int &IntArray::operator()
(const int ind1,const int ind2)const
{return ptr[ ind1 * bound2 + ind2 ];}

const int &IntArray::operator()
(const int ind1,const int ind2,const int ind3)const
{return ptr[ (ind1 * bound2 + ind2) * bound3 + ind3 ];}

const int &IntArray::operator()
(const int ind1,const int ind2,const int ind3,const int ind4)const
{return ptr[ ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4 ];}

const int &IntArray::operator()
(const int ind1,const int ind2,const int ind3,const int ind4,const int ind5)const
{return ptr[ (((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4) * bound5 + ind5 ];}

const int &IntArray::operator()
(const int ind1,const int ind2,const int ind3,const int ind4,const int ind5,const int ind6)const
{return ptr[ ((((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4) * bound5 + ind5) * bound6 + ind6 ];}

//********************************************************
// overloaded subscript operator for non-const Int Array
// const reference return creates an lvakue
//********************************************************
int &IntArray::operator()(const int ind1,const int ind2)
{return ptr[ind1 * bound2 + ind2];}

int &IntArray::operator()(const int ind1,const int ind2,const int ind3)
{return ptr[ (ind1 * bound2 + ind2) * bound3 + ind3 ];}

int &IntArray::operator()(const int ind1,const int ind2,const int ind3,const int ind4)
{return ptr[ ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4 ];}

int &IntArray::operator()
(const int ind1,const int ind2,const int ind3,const int ind4,const int ind5)
{return ptr[ (((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4) * bound5 + ind5 ];}

int &IntArray::operator()
(const int ind1,const int ind2,const int ind3,const int ind4,const int ind5,const int ind6)
{return ptr[ ((((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4) * bound5 + ind5) * bound6 + ind6 ];}

//****************************
// zeroes out the whole array
//****************************
void IntArray::zero_out(void)
{
	if (size <= 0) return;
	for (int i = 0;i < size; i++) ptr[i] = 0;
	return;
}
