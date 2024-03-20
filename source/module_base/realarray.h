/*******************************************
 * ESCP:Electro-Structure Calculate Package.
 ********************************************/

#ifndef REALARRAY_H
#define REALARRAY_H

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>


namespace ModuleBase
{
/**
 * @brief double float array
 *
 */
class realArray
{
  public:
    double *ptr;

    realArray(const int d1 = 1, const int d2 = 1, const int d3 = 1);
    realArray(const int d1, const int d2, const int d3, const int d4);
    ~realArray();

    /**
     * @brief create 3 dimensional real array
     *
     * @param[in] d1 The first dimension size
     * @param[in] d2 The second dimension size
     * @param[in] d3 The third dimension size
     */
    void create(const int d1, const int d2, const int d3);
    void create(const int d1, const int d2, const int d3, const int d4);

    realArray(const realArray &cd);

    /**
     * @brief Equal a realArray to another one
     *
     * @param right
     * @return const realArray&
     */
    const realArray &operator=(const realArray &right);
    /**
     * @brief Set all value of an array to a double float number
     *
     * @param right
     * @return const realArray&
     */
    const realArray &operator=(const double &right);

    /**
     * @brief Access elements by using operator "()"
     *
     * @param d1
     * @param d2
     * @param d3
     * @return double&
     */
    double &operator()(const int d1, const int d2, const int d3);
    double &operator()(const int d1, const int d2, const int d3, const int d4);

    /**
     * @brief Access elements by using "()" through pointer
     * without changing its elements
     *
     * @param d1
     * @param d2
     * @param d3
     * @return const double&
     */
    const double &operator()(const int d1, const int d2, const int d3) const;
    const double &operator()(const int d1, const int d2, const int d3, const int d4) const;

    /**
     * @brief Set all elements of an IntArray to zero
     *
     */
    void zero_out(void);

    /**
     * @brief Get the Size object
     *
     * @return int
     */
    int getSize() const
    {
        return size;
    }

    /**
     * @brief Get the Dim object
     * i.e. the dimension of a real array
     *
     * @return int
     */
    int getDim() const
    {
        return dim;
    }

    /**
     * @brief Get the Bound1 object
     * i.e. the first dimension size
     *
     * @return int
     */
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
    int bound1, bound2, bound3, bound4;
    static int arrayCount;

    void freemem();
};

//**************************************************
// set elements of a as zeros which a is 1_d array.
//**************************************************
template <class T> void zeros(T *u, const int n)
{
    assert(n > 0);
    for (int i = 0; i < n; i++)
        u[i] = 0;
}

} // namespace ModuleBase

#endif // realArray class
