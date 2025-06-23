/*******************************************
 * ESCP:Electro-Structure Calculate Package.
 ********************************************/

#include "realarray.h"

#include <cstdlib>

namespace ModuleBase
{

int realArray::arrayCount = 0;

void realArrayAlloc()
{
	std::cout << "\n Allocation error for realArray " << std::endl;
	exit(0);
}

realArray::realArray(const int d1,const int d2,const int d3)
{
	dim = 3;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = 0;

	size = bound1 * bound2 * bound3 ;	//* sizeof(float);

	auto handler_old = std::set_new_handler(realArrayAlloc);
	ptr = new double[size];
	std::set_new_handler(handler_old);
	zero_out();
	assert(ptr != 0);

	++arrayCount;
}

realArray::realArray(const int d1,const int d2,const int d3,const int d4)
{
	dim = 4;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = (d4 <= 0) ? 1 : d4;

	size = bound1 * bound2 * bound3 * bound4 ;	//* sizeof(float);

	auto handler_old = std::set_new_handler(realArrayAlloc);
	ptr = new double[size];
	std::set_new_handler(handler_old);
	zero_out();

	++arrayCount;
}

realArray::realArray(const realArray &cd)
{
	this->size = cd.getSize();
	this->ptr = new double[size];
	for (int i = 0; i < size; i++)
		this->ptr[i] = cd.ptr[i];
	this->dim = cd.dim;
	this->bound1 = cd.bound1;
	this->bound2 = cd.bound2;
	this->bound3 = cd.bound3;
	this->bound4 = cd.bound4;

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
	ptr = new double[size];

	zero_out(); // mohan modify 2009-09-17
	
	assert(ptr != 0);
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
    ptr = new double[size];
    zero_out();
    assert(ptr != 0);
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
// overloaded subscript operator for const real Array
// const reference return create an cvakue
//********************************************************
// Peize Lin add assert 2016-08-22
const double &realArray::operator()
(const int ind1,const int ind2,const int ind3)const
{
	assert(ind1>=0);	assert(ind1<bound1);
	assert(ind2>=0);	assert(ind2<bound2);
	assert(ind3>=0);	assert(ind3<bound3);	
	return ptr[ (ind1 * bound2 + ind2) * bound3 + ind3 ];
}

const double &realArray::operator()
(const int ind1,const int ind2,const int ind3,const int ind4)const
{
	assert(ind1>=0);	assert(ind1<bound1);
	assert(ind2>=0);	assert(ind2<bound2);
	assert(ind3>=0);	assert(ind3<bound3);
	assert(ind4>=0);	assert(ind4<bound4);	
	return ptr[ ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4 ];
}

//********************************************************
// overloaded subscript operator for non-const real Array
// const reference return creates an lvakue
//********************************************************
double &realArray::operator()(const int ind1,const int ind2,const int ind3)
{
	assert(ind1>=0);	assert(ind1<bound1);
	assert(ind2>=0);	assert(ind2<bound2);
	assert(ind3>=0);	assert(ind3<bound3);
	return ptr[ (ind1 * bound2 + ind2) * bound3 + ind3 ];
}

double &realArray::operator()(const int ind1,const int ind2,const int ind3,const int ind4)
{
	assert(ind1>=0);	assert(ind1<bound1);
	assert(ind2>=0);	assert(ind2<bound2);
	assert(ind3>=0);	assert(ind3<bound3);
	assert(ind4>=0);	assert(ind4<bound4);
	return ptr[ ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4 ];
}

//****************************
// zeroes out the whole array
//****************************
void realArray::zero_out(void)
{
	if (size <= 0) return;
	for (int i = 0;i < size; i++) ptr[i] = 0;
	return;
}

}
