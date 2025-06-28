/*******************************************
 * ESCP:Electro-Structure Calculate Package.
 ********************************************/

#ifndef INTARRAY_H
#define INTARRAY_H

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace ModuleBase
{
/**
 * @brief Integer array
 *
 */
class IntArray
{
  public:
    int *ptr;

    /**
     * @brief Construct a new Int Array object
     *
     * @param d1 The first dimension size
     * @param d2 The second dimension size
     */
    IntArray(const int d1 = 1, const int d2 = 1);
    IntArray(const int d1, const int d2, const int d3);
    IntArray(const int d1, const int d2, const int d3, const int d4);
    IntArray(const int d1, const int d2, const int d3, const int d4, const int d5);
    IntArray(const int d1, const int d2, const int d3, const int d4, const int d5, const int d6);

    ~IntArray();

    /**
     * @brief Create integer arrays
     *
     * @param[in] d1
     * @param[in] d2
     */
    void create(const int d1, const int d2);
    void create(const int d1, const int d2, const int d3);
    void create(const int d1, const int d2, const int d3, const int d4);
    void create(const int d1, const int d2, const int d3, const int d4, const int d5);
    void create(const int d1, const int d2, const int d3, const int d4, const int d5, const int d6);

    /**
     * @brief copy assignment
     *
     * @param right
     * @return const IntArray&
     */
    IntArray &operator=(const IntArray &other)
    {
        if(this != &other)
        {
            delete[] ptr;
            size = other.size;
            dim = other.dim;
            bound1 = other.bound1;
            bound2 = other.bound2;
            bound3 = other.bound3;
            bound4 = other.bound4;
            bound5 = other.bound5;
            bound6 = other.bound6;
            ptr = new int[size];
            for (int i = 0;i < size;i++)
            { ptr[i] = other.ptr[i]; }
        }
        return *this;
    }

    /**
     * @brief Equal all elements of an IntArray to an
     * integer
     *
     * @param right
     * @return const IntArray&
     */
    const IntArray &operator=(const int &right)
    {
        for (int i = 0;i < size;i++) ptr[i] = right;
        return *this;// enables x = y = z;
    }

    /**
     * @brief Access elements by using operator "()"
     *
     * @param d1
     * @param d2
     * @return int&
     */
    int &operator()(const int d1, const int d2)
    {
        assert( d1 < bound1 );
        assert( d2 < bound2 );
        return ptr[ d1 * bound2 + d2 ];
    }
    int &operator()(const int d1, const int d2, const int d3)
    {
        assert( d1 < bound1 );
        assert( d2 < bound2 );
        assert( d3 < bound3 );
        return ptr[ (d1 * bound2 + d2) * bound3 + d3 ];
    }
    int &operator()(const int d1, const int d2, const int d3, const int d4)
    {
        assert( d1 < bound1 );
        assert( d2 < bound2 );
        assert( d3 < bound3 );
        assert( d4 < bound4 );
        return ptr[ ((d1 * bound2 + d2) * bound3 + d3) * bound4 + d4 ];
    }
    int &operator()(const int d1, const int d2, const int d3, const int d4, const int d5)
    {
        assert( d1 < bound1 );
        assert( d2 < bound2 );
        assert( d3 < bound3 );
        assert( d4 < bound4 );
        assert( d5 < bound5 );
        return ptr[ (((d1 * bound2 + d2) * bound3 + d3) * bound4 + d4) * bound5 + d5 ];
    }
    int &operator()(const int d1, const int d2, const int d3, const int d4, const int d5, const int d6)
    {
        assert( d1 < bound1 );
        assert( d2 < bound2 );
        assert( d3 < bound3 );
        assert( d4 < bound4 );
        assert( d5 < bound5 );
        assert( d6 < bound6 );
        return ptr[ ((((d1 * bound2 + d2) * bound3 + d3) * bound4 + d4) * bound5 + d5) * bound6 + d6 ];
    }

    /**
     * @brief Access elements by using "()" through pointer
     * without changing its elements
     *
     * @param d1
     * @param d2
     * @return const int&
     */
    const int &operator()(const int d1, const int d2) const
    {
        assert( d1 < bound1 );
        assert( d2 < bound2 );
        return ptr[ d1 * bound2 + d2 ];
    }
    const int &operator()(const int d1, const int d2, const int d3) const
    {
        assert( d1 < bound1 );
        assert( d2 < bound2 );
        assert( d3 < bound3 );
        return ptr[ (d1 * bound2 + d2) * bound3 + d3 ];
    }
    const int &operator()(const int d1, const int d2, const int d3, const int d4) const
    {
        assert( d1 < bound1 );
        assert( d2 < bound2 );
        assert( d3 < bound3 );
        assert( d4 < bound4 );
        return ptr[ ((d1 * bound2 + d2) * bound3 + d3) * bound4 + d4 ];
    }
    const int &operator()(const int d1, const int d2, const int d3, const int d4, const int d5) const
    {
        assert( d1 < bound1 );
        assert( d2 < bound2 );
        assert( d3 < bound3 );
        assert( d4 < bound4 );
        assert( d5 < bound5 );
        return ptr[ (((d1 * bound2 + d2) * bound3 + d3) * bound4 + d4) * bound5 + d5 ];
    }
    const int &operator()(const int d1, const int d2, const int d3, const int d4, const int d5, const int d6) const
    {
        assert( d1 < bound1 );
        assert( d2 < bound2 );
        assert( d3 < bound3 );
        assert( d4 < bound4 );
        assert( d5 < bound5 );
        assert( d6 < bound6 );
        return ptr[ ((((d1 * bound2 + d2) * bound3 + d3) * bound4 + d4) * bound5 + d5) * bound6 + d6 ];
    }

    /**
     * @brief Set all elements of an IntArray to zero
     *
     */
    void zero_out(void);

    int getSize() const
    {
        return size;
    }
    int getDim() const
    {
        return dim;
    }
    int getBound1() const
    {
        return bound1;
    }
    int getBound2() const
    {
        return bound2;
    }
    int getBound3() const
    {
        return bound3;
    }
    int getBound4() const
    {
        return bound4;
    }
    int getBound5() const
    {
        return bound5;
    }
    int getBound6() const
    {
        return bound6;
    }

  private:
    int size=0;
    int dim=0;
    int bound1=0;
    int bound2=0; 
    int bound3=0;
    int bound4=0; 
    int bound5=0; 
    int bound6=0;
    void freemem();
};
} // namespace ModuleBase

#endif // IntArray class
