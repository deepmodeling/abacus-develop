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
     * @brief Equal an IntArray to another one
     *
     * @param right
     * @return const IntArray&
     */
    const IntArray &operator=(const IntArray &right);

    /**
     * @brief Equal all elements of an IntArray to an
     * integer
     *
     * @param right
     * @return const IntArray&
     */
    const IntArray &operator=(const int &right);

    /**
     * @brief Access elements by using operator "()"
     *
     * @param d1
     * @param d2
     * @return int&
     */
    int &operator()(const int d1, const int d2);
    int &operator()(const int d1, const int d2, const int d3);
    int &operator()(const int d1, const int d2, const int d3, const int d4);
    int &operator()(const int d1, const int d2, const int d3, const int d4, const int d5);
    int &operator()(const int d1, const int d2, const int d3, const int d4, const int d5, const int d6);

    /**
     * @brief Access elements by using "()" through pointer
     * without changing its elements
     *
     * @param d1
     * @param d2
     * @return const int&
     */
    const int &operator()(const int d1, const int d2) const;
    const int &operator()(const int d1, const int d2, const int d3) const;
    const int &operator()(const int d1, const int d2, const int d3, const int d4) const;
    const int &operator()(const int d1, const int d2, const int d3, const int d4, const int d5) const;
    const int &operator()(const int d1, const int d2, const int d3, const int d4, const int d5, const int d6) const;

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

    /**
     * @brief Get the Array Count object
     *
     * @return int
     */
    static int getArrayCount(void)
    {
        return arrayCount;
    }

  private:
    int size;
    int dim;
    int bound1, bound2, bound3, bound4, bound5, bound6;
    static int arrayCount;
    void freemem();
};
} // namespace ModuleBase

#endif // IntArray class
