/*******************************************
 * ESCP:Electro-Structure Calculate Package.
 ********************************************/

#include <cstdlib>
#include "complexarray.h"

int ComplexArray::arrayCount = 0;

void ComplexArrayAlloc()
{
	cout << "\n Allocation error for ComplexArray " << endl;
	exit(0);
}

ComplexArray::ComplexArray(const int d1,const int d2,const int d3)
{
	dim = 3;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = 0;

	set_new_handler(ComplexArrayAlloc);

	size = bound1 * bound2 * bound3 ;	//* sizeof(float);

	ptr = new complex<double>[size]();
	assert(ptr != 0);

	++arrayCount;
}

ComplexArray::ComplexArray(const int d1,const int d2,const int d3,const int d4)
{
	dim = 4;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = (d4 <= 0) ? 1 : d4;

	set_new_handler(ComplexArrayAlloc);

	size = bound1 * bound2 * bound3 * bound4 ;	//* sizeof(float);

	ptr = new complex<double>[size]();
	assert(ptr != 0);

	++arrayCount;
}

//********************************
//
// Destructor for class ComplexArray
//
//********************************
ComplexArray ::~ComplexArray()
{
    freemem();
}

void ComplexArray::freemem()
{
	delete [] ptr;
	ptr = NULL;
}

void ComplexArray::create(const int d1,const int d2,const int d3,const int d4)
{
	size = d1 * d2 * d3 * d4;
	assert(size>0);

	dim = 4;

	bound1 = d1;
	bound2 = d2;
	bound3 = d3;
	bound4 = d4;

	delete [] ptr;
	ptr = new complex<double>[size]();
	assert(ptr != 0);
}

void ComplexArray::create(const int d1,const int d2,const int d3)
{
	size = d1 * d2 * d3;
	assert(size>0);

	dim = 3;

	bound1 = d1;
	bound2 = d2;
	bound3 = d3;
	bound4 = 1;

	delete [] ptr;
	ptr = new complex<double>[size]();
	assert(ptr != 0);
}

const ComplexArray &ComplexArray::operator=(const ComplexArray &right)
{
	for (int i = 0;i < size;i++) ptr[i] = right.ptr[i];
	return *this;// enables x = y = z;
}

const ComplexArray &ComplexArray::operator=(const complex<double> &right)
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
const complex<double> &ComplexArray::operator()
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
const complex<double> &ComplexArray::operator()
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
complex<double> &ComplexArray::operator()(const int ind1,const int ind2,const int ind3)
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
complex<double> &ComplexArray::operator()(const int ind1,const int ind2,const int ind3,const int ind4)
{
	const int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;
	return ptr[ind];// reference return
}

//****************************
//
// zeroes out the whole array
//
//****************************
void ComplexArray::zero_out(void)
{
	if (size <= 0) return;
	for (int i = 0;i < size; i++) ptr[i] = 0;
	return;
}
