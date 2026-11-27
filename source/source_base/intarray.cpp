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
    try 
    {
        ptr = new int[size];
        zero_out();
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error for IntArray: " << e.what() << std::endl;
        ptr = nullptr;
        size = 0;
        throw;
    }
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
    size = bound1 * bound2 * bound3 ;    //* sizeof(float);
    try 
    {
        ptr = new int[size];
        zero_out();
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error for IntArray: " << e.what() << std::endl;
        ptr = nullptr;
        size = 0;
        throw;
    }
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
    size = bound1 * bound2 * bound3 * bound4 ;    //* sizeof(float);
    try 
    {
        ptr = new int[size];
        zero_out();
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error for IntArray: " << e.what() << std::endl;
        ptr = nullptr;
        size = 0;
        throw;
    }
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
    try 
    {
        ptr = new int[size];
        zero_out();
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error for IntArray: " << e.what() << std::endl;
        ptr = nullptr;
        size = 0;
        throw;
    }
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
    try 
    {
        ptr = new int[size];
        zero_out();
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error for IntArray: " << e.what() << std::endl;
        ptr = nullptr;
        size = 0;
        throw;
    }
    assert(ptr != nullptr);
}

// Copy constructor
IntArray::IntArray(const IntArray& other)
{
    size = other.size;
    dim = other.dim;
    bound1 = other.bound1;
    bound2 = other.bound2;
    bound3 = other.bound3;
    bound4 = other.bound4;
    bound5 = other.bound5;
    bound6 = other.bound6;
    try 
    {
        ptr = new int[size];
        for (int i = 0; i < size; i++) 
        {
            ptr[i] = other.ptr[i];
        }
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error in IntArray copy constructor: " << e.what() << std::endl;
        ptr = nullptr;
        size = 0;
        throw;
    }
    assert(ptr != nullptr);
}

// Move constructor
IntArray::IntArray(IntArray&& other) noexcept
    : size(other.size),
      dim(other.dim),
      bound1(other.bound1),
      bound2(other.bound2),
      bound3(other.bound3),
      bound4(other.bound4),
      bound5(other.bound5),
      bound6(other.bound6),
      ptr(other.ptr)
{
    other.ptr = nullptr;
    other.size = 0;
    other.dim = 0;
    other.bound1 = other.bound2 = other.bound3 = other.bound4 = other.bound5 = other.bound6 = 0;
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
    if (ptr != nullptr)
    {
        delete[] ptr;
        ptr = nullptr;
    }
}

// Move assignment operator
IntArray& IntArray::operator=(IntArray&& other) noexcept
{
    if (this != &other)
    {
        freemem();
        size = other.size;
        dim = other.dim;
        bound1 = other.bound1;
        bound2 = other.bound2;
        bound3 = other.bound3;
        bound4 = other.bound4;
        bound5 = other.bound5;
        bound6 = other.bound6;
        ptr = other.ptr;
        other.ptr = nullptr;
        other.size = 0;
        other.dim = 0;
        other.bound1 = other.bound2 = other.bound3 = other.bound4 = other.bound5 = other.bound6 = 0;
    }
    return *this;
}

void IntArray::create(const int d1,const int d2,const int d3,const int d4,const int d5,const int d6)
{
    size = d1 * d2 * d3 * d4 * d5 * d6;assert(size>0);
    dim = 6;
    bound1 = d1;bound2 = d2;bound3 = d3;bound4 = d4;bound5 = d5;bound6 = d6;
    int* new_ptr = nullptr;
    try 
    {
        new_ptr = new int[size];
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error in IntArray::create: " << e.what() << std::endl;
        assert(new_ptr != nullptr);
        return;
    }
    delete[] ptr;
    ptr = new_ptr;
    zero_out();
}

void IntArray::create(const int d1,const int d2,const int d3,const int d4,const int d5)
{
    size = d1 * d2 * d3 * d4 * d5;assert(size>0);
    dim = 5;
    bound1 = d1;bound2 = d2;bound3 = d3;bound4 = d4;bound5 = d5;
    int* new_ptr = nullptr;
    try 
    {
        new_ptr = new int[size];
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error in IntArray::create: " << e.what() << std::endl;
        assert(new_ptr != nullptr);
        return;
    }
    delete[] ptr;
    ptr = new_ptr;
    zero_out();
}

void IntArray::create(const int d1,const int d2,const int d3,const int d4)
{
    size = d1 * d2 * d3 * d4;assert(size>0);
    dim = 4;
    bound1 = d1;bound2 = d2;bound3 = d3;bound4 = d4;
    int* new_ptr = nullptr;
    try 
    {
        new_ptr = new int[size];
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error in IntArray::create: " << e.what() << std::endl;
        assert(new_ptr != nullptr);
        return;
    }
    delete[] ptr;
    ptr = new_ptr;
    zero_out();
}

void IntArray::create(const int d1,const int d2,const int d3)
{
    size = d1 * d2 * d3;assert(size>0);
    dim = 3;
    bound1 = d1;bound2 = d2;bound3 = d3;bound4 = 1;
    int* new_ptr = nullptr;
    try 
    {
        new_ptr = new int[size];
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error in IntArray::create: " << e.what() << std::endl;
        assert(new_ptr != nullptr);
        return;
    }
    delete [] ptr;
    ptr = new_ptr;
    zero_out();
}

void IntArray::create(const int d1, const int d2)
{
    size = d1 * d2;assert(size>0);
    dim = 2;
    bound1 = d1;bound2 = d2;bound3 = bound4 = 1;
    int* new_ptr = nullptr;
    try 
    {
        new_ptr = new int[size];
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << "Allocation error in IntArray::create: " << e.what() << std::endl;
        assert(new_ptr != nullptr);
        return;
    }
    delete[] ptr;
    ptr = new_ptr;
    zero_out();
}

//****************************
// zeroes out the whole array
//****************************
void IntArray::zero_out()
{
    if (size <= 0 || ptr == nullptr)
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
