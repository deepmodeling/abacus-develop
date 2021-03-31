//==========================================================
// AUTHOR: Lixin He, Mohan Chen
// LAST UPDATE : 2009-03-24 add operator == and !=
//==========================================================
#ifndef MatriX3_H
#define MatriX3_H

#ifdef _MCD_CHECK
#include "../src_parallel/mcd.h"
#endif

#include "src_global/vector3.h"
#include "src_global/matrix.h"

class Matrix3
{
	/* data */
public:
	// element eij:i_column,j_row
	double e11, e12, e13, e21, e22, e23, e31, e32, e33;

	/* Constructors and destructor */
	Matrix3(const double &r11 = 1,const double &r12 = 0,const double &r13 = 0,
	        const double &r21 = 0,const double &r22 = 1,const double &r23 = 0,
	        const double &r31 = 0,const double &r32 = 0,const double &r33 = 1);

	void Reset(void);
	void Identity(void);
	double Det(void) const ;
	Matrix3	Transpose(void) const ;
	Matrix3	Inverse(void) const ;

	Matrix3& operator=(const Matrix3 &m);
	Matrix3& operator+=(const Matrix3 &m);
	Matrix3& operator-=(const Matrix3 &m);
	Matrix3& operator*=(const double &s);
	Matrix3& operator/=(const double &s);

	void print(void)const;
	matrix to_matrix(void)const;
};

Matrix3 operator +(const Matrix3 &m1, const Matrix3 &m2);	//m1+m2
Matrix3 operator -(const Matrix3 &m1, const Matrix3 &m2);	//m1-m2
Matrix3 operator /(const Matrix3 &m,const double &s);		//m/s
Matrix3 operator *(const Matrix3 &m1,const  Matrix3 &m2);	//m1*m2
Matrix3 operator *(const Matrix3 &m,const double &s);		//m*s
Matrix3 operator *(double &s, const Matrix3 &m);		//s*m
template<typename T> Vector3<double> operator *(const Matrix3 &m, const Vector3<T> &u);	//m*u				// Peize Lin change Vector3<T> 2017-01-10
template<typename T> Vector3<double> operator *(const Vector3<T> &u, const Matrix3 &m);	//u*m				// Peize Lin change Vector3<T> 2017-01-10

bool operator ==(const Matrix3 &m1, const Matrix3 &m2); //whether m1 == m2
bool operator !=(const Matrix3 &m1, const Matrix3 &m2); //whethor m1 != m2




//m*u
template<typename T>
Vector3<double> operator *(const Matrix3 &m, const Vector3<T> &u)
{
	return Vector3<double>(m.e11*u.x + m.e12*u.y + m.e13*u.z,
	                       m.e21*u.x + m.e22*u.y + m.e23*u.z,
	                       m.e31*u.x + m.e32*u.y + m.e33*u.z);
}

//u*m
template<typename T>
Vector3<double> operator *(const Vector3<T> &u, const Matrix3 &m)
{
	return Vector3<double>(u.x*m.e11 + u.y*m.e21 + u.z*m.e31,
	                       u.x*m.e12 + u.y*m.e22 + u.z*m.e32,
	                       u.x*m.e13 + u.y*m.e23 + u.z*m.e33);
}

#endif // MATRIX3_H

