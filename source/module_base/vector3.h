#ifndef VECTOR3_H
#define VECTOR3_H

#ifdef _MCD_CHECK
#include "mcd.h"
#endif

#include <cmath>
#include <iomanip>
#include <iostream>
#include <array>
namespace ModuleBase
{

/**
 * @brief 3 elements vector
 *
 * @tparam T
 */
template <class T> class Vector3
{
  public:
    T x;
    T y;
    T z;

    /**
     * @brief Construct a new Vector 3 object
     *
     * @param x1
     * @param y1
     * @param z1
     */
    Vector3(const T &x1 = 0, const T &y1 = 0, const T &z1 = 0) : x(x1), y(y1), z(z1){};
    Vector3(const Vector3<T> &v) : x(v.x), y(v.y), z(v.z){}; // Peize Lin add 2018-07-16
	explicit Vector3(const std::array<T,3> &v) :x(v[0]), y(v[1]), z(v[2]){}

    /**
     * @brief set a 3d vector
     *
     * @param x1
     * @param y1
     * @param z1
     */
    void set(const T &x1, const T &y1, const T &z1)
    {
        x = x1;
        y = y1;
        z = z1;
    }

    /**
     * @brief Overload operator "=" for Vector3
     *
     * @param u
     * @return Vector3<T>&
     */
    Vector3<T> &operator=(const Vector3<T> &u)
    {
        x = u.x;
        y = u.y;
        z = u.z;
        return *this;
    }

    Vector3<T> &operator=(const T &u)
    {
        x = u;
        y = u;
        z = u;
        return *this;
    }

    /**
     * @brief Overload operator "+=" for Vector3
     *
     * @param u
     * @return Vector3<T>&
     */
    Vector3<T> &operator+=(const Vector3<T> &u)
    {
        x += u.x;
        y += u.y;
        z += u.z;
        return *this;
    }

    /**
     * @brief Overload operator "-=" for Vector3
     *
     * @param u
     * @return Vector3<T>&
     */
    Vector3<T> &operator-=(const Vector3<T> &u)
    {
        x -= u.x;
        y -= u.y;
        z -= u.z;
        return *this;
    }

    /**
     * @brief Overload operator "*=" for (Vector3)*scalar
     *
     * @param s
     * @return Vector3<T>&
     */
    Vector3<T> &operator*=(const T &s)
    {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    /**
     * @brief Overload operator "/=" for (Vector3)/scalar
     *
     * @param s
     * @return Vector3<T>&
     */
    Vector3<T> &operator/=(const T &s)
    {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }

    /**
     * @brief Overload operator "-" to get (-Vector3)
     *
     * @return Vector3<T>
     */
    Vector3<T> operator-() const
    {
        return Vector3<T>(-x, -y, -z);
    } // Peize Lin add 2017-01-10

    /**
     * @brief Over load "[]" for accessing elements with pointers
     *
     * @param index
     * @return T
     */
    T operator[](int index) const
    {
        //return (&x)[index]; // this is undefind behavior and breaks with icpx
        T const* ptr[3] = {&x, &y, &z};
        return *ptr[index];
    }

    /**
     * @brief Overload operator "[]" for accesing elements
     *
     * @param index
     * @return T&
     */
    T &operator[](int index)
    {
        //return (&x)[index]; // this is undefind behavior and breaks with icpx
        T* ptr[3] = {&x, &y, &z};
        return *ptr[index];
    }

    /**
     * @brief Get the square of nomr of a Vector3
     *
     * @return T
     */
    T norm2(void) const
    {
        return x * x + y * y + z * z;
    }

    /**
     * @brief Get the norm of a Vector3
     *
     * @return T
     */
    T norm(void) const
    {
        return sqrt(norm2());
    }

    /**
     * @brief Normalize a Vector3
     *
     * @return Vector3<T>&
     */
    Vector3<T> &normalize(void)
    {
        const T m = norm();
        x /= m;
        y /= m;
        z /= m;
        return *this;
    } // Peize Lin update return 2019-09-08

    /**
     * @brief Get (-Vector3)
     *
     * @return Vector3<T>&
     */
    Vector3<T> &reverse(void)
    {
        x = -x;
        y = -y;
        z = -z;
        return *this;
    } // Peize Lin update return 2019-09-08

    /**
     * @brief Print a Vector3 on standard output
     * with formats
     *
     */
    void print(void) const; // mohan add 2009-11-29
};

/**
 * @brief Overload "+" for two Vector3
 *
 * @param[in] u
 * @param[in] v
 * @return Vector3<T>
 */
template <class T> inline Vector3<T> operator+(const Vector3<T> &u, const Vector3<T> &v)
{
    return Vector3<T>(u.x + v.x, u.y + v.y, u.z + v.z);
}

/**
 * @brief Overload "-" for two Vector3
 *
 * @param[in] u
 * @param[in] v
 * @return Vector3<T>
 */
template <class T> inline Vector3<T> operator-(const Vector3<T> &u, const Vector3<T> &v)
{
    return Vector3<T>(u.x - v.x, u.y - v.y, u.z - v.z);
}

/**
 * @brief Overload "*" to calculate the dot product
 * of two Vector3
 *
 * @param u
 * @param v
 * @return template <class T>
 */
template <class T> inline T operator*(const Vector3<T> &u, const Vector3<T> &v)
{
    return (u.x * v.x + u.y * v.y + u.z * v.z);
}

/**
 * @brief Overload "*" to calculate (Vector3)*scalar
 *
 * @param[in] s
 * @param[in] u
 * @return Vector3<T>
 */
template <class T> inline Vector3<T> operator*(const T &s, const Vector3<T> &u)
{
    return Vector3<T>(u.x * s, u.y * s, u.z * s);
}

/**
 * @brief Overload "*" to calculate scalar*(Vector3)
 *
 * @param u
 * @param s
 * @return Vector3<T>
 */
template <class T> inline Vector3<T> operator*(const Vector3<T> &u, const T &s)
{
    return Vector3<T>(u.x * s, u.y * s, u.z * s);
} // mohan add 2009-5-10

/**
 * @brief Overload "/" to calculate Vector3/scalar
 *
 * @tparam T
 * @param u
 * @param s
 * @return Vector3<T>
 */
template <class T> inline Vector3<T> operator/(const Vector3<T> &u, const T &s)
{
    return Vector3<T>(u.x / s, u.y / s, u.z / s);
}

/**
 * @brief Overload "/" to calculate scalar/Vector3
 *
 * @tparam T
 * @param s
 * @param u
 * @return Vector3<T>
 */
template <class T> inline Vector3<T> operator/(const T &s, const Vector3<T> &u)
{
    return Vector3<T>(s/u.x, s/u.y, s/u.z);
}

/**
 * @brief Dot productor of two Vector3
 *
 * @param u
 * @param v
 * @return T
 * @note u.v=(ux*vx)+(uy*vy)+(uz*vz)
 */
template <class T> inline T dot(const Vector3<T> &u, const Vector3<T> &v)
{
    return (u.x * v.x + u.y * v.y + u.z * v.z);
}

/**
 * @brief Overload "^" for cross product of two Vector3
 *
 * @param u
 * @param v
 * @return template <class T>
 * @note
 * | i  j  k  |
 * | ux uy uz |
 * | vx vy vz |
 * u.v=(uy*vz-uz*vy)i+(-ux*vz+uz*vx)j+(ux*vy-uy*vx)k
 */
template <class T> inline Vector3<T> operator^(const Vector3<T> &u, const Vector3<T> &v)
{
    return Vector3<T>(u.y * v.z - u.z * v.y, -u.x * v.z + u.z * v.x, u.x * v.y - u.y * v.x);
}

/**
 * @brief Cross product of two Vector3
 *
 * @param u
 * @param v
 * @return template <class T>
 * @note
 * | i  j  k  |
 * | ux uy uz |
 * | vx vy vz |
 * u.v=(uy*vz-uz*vy)i+(-ux*vz+uz*vx)j+(ux*vy-uy*vx)k
 */
template <class T> inline Vector3<T> cross(const Vector3<T> &u, const Vector3<T> &v)
{
    return Vector3<T>(u.y * v.z - u.z * v.y, -u.x * v.z + u.z * v.x, u.x * v.y - u.y * v.x);
}
// s = u.(v x w)
// template <class T> T TripleScalarProduct(Vector3<T> u, Vector3<T> v, Vector3<T> w)
//{
//	return T((u.x * (v.y * w.z - v.z * w.y)) +
//	         (u.y * (-v.x * w.z + v.z * w.x)) +
//	         (u.z * (v.x * w.y - v.y * w.x)));
// }

// whether m1 != m2
template <class T> inline bool operator!=(const Vector3<T> &u, const Vector3<T> &v)
{
    return !(u == v);
}
// whether u == v
template <class T> inline bool operator==(const Vector3<T> &u, const Vector3<T> &v)
{
    if (u.x == v.x && u.y == v.y && u.z == v.z)
        return true;
    return false;
}

/**
 * @brief Print a Vector3 on standard output with formats
 *
 */
template <class T> void Vector3<T>::print(void) const
{
    std::cout.precision(5);
    std::cout << "(" << std::setw(10) << x << "," << std::setw(10) << y << "," << std::setw(10) << z << ")"
              << std::endl;
    return;
}

/**
 * @brief Overload "<<" tor print out a
 * Vector3 on standard output
 *
 * @tparam T
 * @param[in] os
 * @param[in] u
 * @return std::ostream&
 */
template <class T> static std::ostream &operator<<(std::ostream &os, const Vector3<T> &u)
{
    os << "(" << std::setw(10) << u.x << "," << std::setw(10) << u.y << "," << std::setw(10) << u.z << ")";
    return os;
}

} // namespace ModuleBase

#endif
