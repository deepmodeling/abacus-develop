//==========================================================
//
//==========================================================
#ifndef VECTOR3_H
#define VECTOR3_H

#ifdef _MCD_CHECK
#include "./src_parallel/mcd.h"
#endif

#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

template <class T>
class Vector3
{
public:
	T x;
	T y;
	T z;
	Vector3(const T &x1 = 0,const T &y1 = 0,const T &z1 = 0);

	Vector3<T>& operator=(const Vector3<T> &u);
	Vector3<T>& operator+=(const Vector3<T> &u);
	Vector3<T>& operator-=(const Vector3<T> &u);
	Vector3<T>& operator*=(const Vector3<T> &u);
	Vector3<T>& operator*=(const T &s);
	Vector3<T>& operator/=(const Vector3<T> &u);
	Vector3<T>& operator/=(const T &s);

	void set(const T &x1, const T &y1,const T &z1);

	inline T norm(void) const
	{ return sqrt(x*x + y*y + z*z); }

	T norm2(void)const;
	void normalize(void);
	void reverse(void);

	// mohan add 2009-11-29
	void print(void)const ;
};

template <class T>
Vector3<T> operator+(const Vector3<T> &u,const Vector3<T> &v);
template <class T>
Vector3<T> operator-(const Vector3<T> &u,const Vector3<T> &v);
// | i  j  k  |
// | ux uy uz |
// | vx vy vz |
//u.v=(uy*vz-uz*vy)i+(-ux*vz+uz*vx)j+(ux*vy-uy*vx)k
template <class T>
Vector3<T> operator^(const Vector3<T> &u,const Vector3<T> &v);
//u.v=(ux*vx)+(uy*vy)+(uz*vz)
template <class T>
T operator*(const Vector3<T> &u,const Vector3<T> &v);
template <class T>
Vector3<T> operator*(const T &s,const Vector3<T> &u);
template <class T>
Vector3<T> operator/(const Vector3<T> &u,const T &s);
template <class T>
Vector3<T> cross(const Vector3<T> &u,const Vector3<T> &v);
//u.v=(ux*vx)+(uy*vy)+(uz*vz)
template <class T>
T dot(const Vector3<T> &u,const Vector3<T> &v);
//template <class T>
//T TripleSalarProduct(Vector3<T> u, Vector3<T> v, Vector3<T> w);

//==========================================================
// Define Constructor
//==========================================================
template <class T>
Vector3<T>::Vector3(const T &x1,const T &y1, const T &z1)
:x(x1),y(y1),z(z1)
{}

template <class T>
void Vector3<T>::set(const T &x1,const T &y1,const T &z1)
{
	x = x1; y = y1; z = z1;
}

template <class T>
T Vector3<T>::norm2(void)const
{
	return x*x + y*y + z*z;
}

//|v| = sqrt(x*x+y*y+z*z)=1
template <class T>
void Vector3<T>::normalize(void)
{
	T m = sqrt(x*x + y*y + z*z);
	x /= m;y /= m;z /= m;
}

template <class T>
void Vector3<T>::reverse(void)
{
	x = -x;y = -y;z = -z;
}

template <class T>
Vector3<T>& Vector3<T>::operator=(const Vector3<T> &u)
{
	x = u.x;y = u.y;z = u.z;
	return *this;
}

template <class T>
Vector3<T>& Vector3<T>::operator+=(const Vector3<T> &u)
{
	x += u.x;y += u.y;z += u.z;
	return *this;
}

template <class T>
Vector3<T>& Vector3<T>::operator-=(const Vector3<T> &u)
{
	x -= u.x;y -= u.y;z -= u.z;
	return *this;
}

template <class T>
Vector3<T>& Vector3<T>::operator*=(const T &s)
{
	x *= s;y *= s;z *= s;
	return *this;
}

template <class T>
Vector3<T>& Vector3<T>::operator/=(const T &s)
{
	x /= s;y /= s;z /= s;
	return *this;
}

template <class T>
Vector3<T> operator+(const Vector3<T> &u,const Vector3<T> &v)
{
	return Vector3<T>(u.x + v.x, u.y + v.y, u.z + v.z);
}

template <class T>
Vector3<T> operator-(const Vector3<T> &u,const Vector3<T> &v)
{
	return Vector3<T>(u.x - v.x, u.y - v.y, u.z - v.z);
}

//	u.v=(uy*vz-uz*vy)i+(-ux*vz+uz*vx)j+(ux*vy-uy*vzx)k
template <class T>
Vector3<T> operator^(const Vector3<T> &u,const Vector3<T> &v)
{
	return Vector3<T> (u.y * v.z - u.z * v.y,
	                   -u.x * v.z + u.z * v.x,
	                   u.x * v.y - u.y * v.x);
}

//u.v=(ux*vx)+(uy*vy)+(uz*vz)
template <class T>
T operator*(const Vector3<T> &u,const Vector3<T> &v)
{
	return (u.x * v.x + u.y * v.y + u.z * v.z);
}

template <class T>// mohan add 2009-5-10
Vector3<T> operator*(const Vector3<T> &u, const T &s)
{
	return Vector3<T>(u.x * s, u.y * s, u.z * s);
}

template <class T>
Vector3<T> operator*(const T &s,const Vector3<T> &u)
{
	return Vector3<T>(u.x * s, u.y * s, u.z * s);
}

template <class T>
Vector3<T> operator /(const Vector3<T> &u,const T &s)
{
	return Vector3<T>(u.x / s, u.y / s, u.z / s);
}

//	u.v=(uy*vz-uz*vy)i+(-ux*vz+uz*vx)j+(ux*vy-uy*vzx)k
template <class T>
Vector3<T> cross(const Vector3<T> &u, const Vector3<T> &v)
{
	return Vector3<T> (u.y * v.z - u.z * v.y,
	                   -u.x * v.z + u.z * v.x,
	                   u.x * v.y - u.y * v.x);
}

template <class T>
T dot(const Vector3<T> &u,const Vector3<T> &v)
{
	return (u.x * v.x + u.y * v.y + u.z * v.z);
}

//whether u == v
template <class T>
bool operator ==(const Vector3<T> &u, const Vector3<T> &v)
{
	if(u.x == v.x && u.y == v.y && u.z == v.z)
	{
		return true;
	}
	return false;
}

//whether m1 != m2
template <class T>
bool operator !=(const Vector3<T> &u, const Vector3<T> &v)
{
    return !(u == v);
}

template <class T>
void Vector3<T>::print(void)const
{
	cout.precision(5) ;
	cout << "(" << setw(10) << x << "," << setw(10) << y << ","
	<< setw(10) << z  << ")"  << endl ;
	return ;
}

//s = u.(v x w)
//template <class T>
//T TripleScalarProduct(Vector3<T> u, Vector3<T> v, Vector3<T> w)
//{
//	return T((u.x * (v.y * w.z - v.z * w.y)) +
//	         (u.y * (-v.x * w.z + v.z * w.x)) +
//	         (u.z * (v.x * w.y - v.y * w.x)));
//}

#endif
