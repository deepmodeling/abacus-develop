/*******************************************
 * ESCP:Electro-Structure Calculate Package.
 ********************************************/

#include <cstdlib>
#include "realarray.h"

int realArray::arrayCount = 0;

void realArrayAlloc()
{
	cout << "\n Allocation error for realArray " << endl;
	exit(0);
}

realArray::realArray(const int d1,const int d2,const int d3)
{
	dim = 3;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = 0;

	set_new_handler(realArrayAlloc);

	size = bound1 * bound2 * bound3 ;	//* sizeof(float);

	ptr = new double[size]();
	assert(ptr != 0);
	for(int i=0; i<size; i++) ptr[i] = 0.0;

	++arrayCount;
}

realArray::realArray(const int d1,const int d2,const int d3,const int d4)
{
	dim = 4;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = (d4 <= 0) ? 1 : d4;

	set_new_handler(realArrayAlloc);

	size = bound1 * bound2 * bound3 * bound4 ;	//* sizeof(float);

	ptr = new double[size]();
	assert(ptr != 0);
	for(int i=0; i<size; i++) ptr[i] = 0.0;

	++arrayCount;
}

//********************************
//
// Destructor for class realArray
//
//********************************
realArray ::~realArray()
{
    freemem();
}

void realArray::freemem()
{
	delete [] ptr;
	ptr = NULL;
}

void realArray::create(const int d1,const int d2,const int d3,const int d4)
{
	size = d1 * d2 * d3 * d4;
	assert(size>0);

	dim = 4;

	bound1 = d1;
	bound2 = d2;
	bound3 = d3;
	bound4 = d4;

	delete [] ptr;
	ptr = new double[size]();
	assert(ptr != 0);
	for(int i=0; i<size; i++) ptr[i] = 0.0;
}

void realArray::create(const int d1,const int d2,const int d3)
{
	size = d1 * d2 * d3;
	assert(size>0);

	dim = 3;

	bound1 = d1;
	bound2 = d2;
	bound3 = d3;
	bound4 = 1;

	delete [] ptr;
	ptr = new double[size]();
	assert(ptr != 0);
	for(int i=0; i<size; i++) ptr[i] = 0.0;
}

const realArray &realArray::operator=(const realArray &right)
{
	for (int i = 0;i < size;i++) ptr[i] = right.ptr[i];
	return *this;// enables x = y = z;
}

const realArray &realArray::operator=(const double &right)
{
	for (int i = 0;i < size;i++) ptr[i] = right;
	return *this;// enables x = y = z;
}

//********************************************************
//
// overloaded subscript operator for const real Array
// const reference return create an cvakue
//
//********************************************************
const double &realArray::operator()
(const int ind1,const int ind2,const int ind3)const
{
	const int ind = (ind1 * bound2 + ind2) * bound3 + ind3 ;
	return ptr[ind];
}


//********************************************************
//
// overloaded subscript operator for const real Array
// const reference return creates an cvakue
//
//********************************************************
const double &realArray::operator()
(const int ind1,const int ind2,const int ind3,const int ind4)const
{
	const int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;
	return ptr[ind];
}

//********************************************************
//
// overloaded subscript operator for non-const real Array
// const reference return creates an lvakue
//
//********************************************************
double &realArray::operator()(const int ind1,const int ind2,const int ind3)
{
	const int ind = (ind1 * bound2 + ind2) * bound3 + ind3;
	return ptr[ind]; // reference return
}

//********************************************************
//
// overloaded subscript operator for non-const real Array
// const reference return creates an lvakue
//
//********************************************************
double &realArray::operator()(const int ind1,const int ind2,const int ind3,const int ind4)
{
	const int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;
	return ptr[ind];// reference return
}

//****************************
//
// zeroes out the whole array
//
//****************************
void realArray::zero_out(void)
{
	if (size <= 0) return;
	for (int i = 0;i < size; i++) ptr[i] = 0;
	return;
}
