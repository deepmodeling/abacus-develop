/*******************************************
 * ESCP:Electro-Structure Calculate Package.
 ********************************************/
#include <cstdlib>
#include "intarray.h"

namespace ModuleBase
{
void IntArrayAlloc()
{
	std::cout << "\n Allocation error for IntArray " << std::endl;
	exit(0);
}

IntArray::IntArray(const int d1,const int d2)
{
	dim = 2;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = bound4 = bound5 = bound6 = 0;
	size = bound1 * bound2;
	ptr = new int[size];zero_out();
	assert( ptr != nullptr);
}

IntArray::IntArray(const int d1,const int d2,const int d3)
{
	dim = 3;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = bound5 = bound6 = 0;
	//set_new_handler(IntArrayAlloc);
	size = bound1 * bound2 * bound3 ;	//* sizeof(float);
	ptr = new int[size];zero_out();
	assert(ptr != nullptr);
}

IntArray::IntArray(const int d1,const int d2,const int d3,const int d4)
{
	dim = 4;
	bound1 = (d1 <= 0) ? 1 : d1;
	bound2 = (d2 <= 0) ? 1 : d2;
	bound3 = (d3 <= 0) ? 1 : d3;
	bound4 = (d4 <= 0) ? 1 : d4;
	bound5 = bound6 = 0;
	//set_new_handler(IntArrayAlloc);
	size = bound1 * bound2 * bound3 * bound4 ;	//* sizeof(float);
	ptr = new int[size];zero_out();
	assert(ptr != nullptr);
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
	//set_new_handler(IntArrayAlloc);
	size = bound1 * bound2 * bound3 * bound4 * bound5;
	ptr = new int[size];zero_out();
	assert(ptr != nullptr);
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
    //set_new_handler(IntArrayAlloc);
    size = bound1 * bound2 * bound3 * bound4 * bound5 * bound6;
	ptr = new int[size];zero_out();
	assert(ptr != nullptr);
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
	if(ptr!= nullptr)
	{
		delete [] ptr;
		ptr = nullptr;
	}
}

void IntArray::create(const int d1,const int d2,const int d3,const int d4,const int d5,const int d6)
{
	size = d1 * d2 * d3 * d4 * d5 * d6;assert(size>0);
	dim = 6;
	bound1 = d1;bound2 = d2;bound3 = d3;bound4 = d4;bound5 = d5;bound6 = d6;
	delete[] ptr; ptr = new int[size];
	assert(ptr != nullptr);zero_out();
}

void IntArray::create(const int d1,const int d2,const int d3,const int d4,const int d5)
{
	size = d1 * d2 * d3 * d4 * d5;assert(size>0);
	dim = 5;
	bound1 = d1;bound2 = d2;bound3 = d3;bound4 = d4;bound5 = d5;
	delete[] ptr; ptr = new int[size];
	assert(ptr != nullptr);zero_out();
}

void IntArray::create(const int d1,const int d2,const int d3,const int d4)
{
	size = d1 * d2 * d3 * d4;assert(size>0);
	dim = 4;
	bound1 = d1;bound2 = d2;bound3 = d3;bound4 = d4;
	delete[] ptr; ptr = new int[size];
	assert(ptr != nullptr);zero_out();
}

void IntArray::create(const int d1,const int d2,const int d3)
{
	size = d1 * d2 * d3;assert(size>0);
	dim = 3;
	bound1 = d1;bound2 = d2;bound3 = d3;bound4 = 1;
	delete [] ptr;ptr = new int[size];
	assert(ptr != nullptr);zero_out();
}

void IntArray::create(const int d1, const int d2)
{
	size = d1 * d2;assert(size>0);
	dim = 2;
	bound1 = d1;bound2 = d2;bound3 = bound4 = 1;
	delete[] ptr;ptr = new int[size];
	assert(ptr != nullptr );zero_out();
}

//****************************
// zeroes out the whole array
//****************************
void IntArray::zero_out()
{
	if (size <= 0) 
	{
		return;
	}
	for (int i = 0;i < size; i++)
	{
		ptr[i] = 0;
	}
	return;
}

}