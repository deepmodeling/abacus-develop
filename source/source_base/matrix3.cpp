#include "matrix3.h"

namespace ModuleBase
{

Matrix3::Matrix3(const double &r11, const double &r12, const double &r13,
                 const double &r21, const double &r22, const double &r23,
                 const double &r31, const double &r32, const double &r33)
{
	e11 = r11;e12 = r12;e13 = r13;
	e21 = r21;e22 = r22;e23 = r23;
	e31 = r31;e32 = r32;e33 = r33;
}

void Matrix3::Identity(void)
{
	e11 = 1;e12 = 0;e13 = 0;
	e21 = 0;e22 = 1;e23 = 0;
	e31 = 0;e32 = 0;e33 = 1;
}

void Matrix3::Zero(void)
{
	e11 = 0;e12 = 0;e13 = 0;
	e21 = 0;e22 = 0;e23 = 0;
	e31 = 0;e32 = 0;e33 = 0;
}

double Matrix3::Det(void) const 
{
	return	e11*e22*e33 -
	        e11*e32*e23 +
	        e21*e32*e13 -
	        e21*e12*e33 +
	        e31*e12*e23 -
	        e31*e22*e13;
}

Matrix3 Matrix3::Transpose(void) const
{
	return Matrix3(e11, e21, e31, e12, e22, e32, e13, e23, e33);
}

Matrix3 Matrix3::Inverse(void) const
{
	double d = this->Det();

	if(d == 0)d = 1;

	return Matrix3((e22*e33 - e23*e32) / d,
	               -(e12*e33 - e13*e32) / d,
	               (e12*e23 - e13*e22) / d,
	               -(e21*e33 - e23*e31) / d,
	               (e11*e33 - e13*e31) / d,
	               -(e11*e23 - e13*e21) / d,
	               (e21*e32 - e22*e31) / d,
	               -(e11*e32 - e12*e31) / d,
	               (e11*e22 - e12*e21) / d);
}

Matrix3& Matrix3::operator = (const Matrix3 &m)
{
	e11 = m.e11;e12 = m.e12;e13 = m.e13;
	e21 = m.e21;e22 = m.e22;e23 = m.e23;
	e31 = m.e31;e32 = m.e32;e33 = m.e33;
	return *this;
}

Matrix3& Matrix3::operator +=(const Matrix3 &m)
{
	e11 += m.e11;e12 += m.e12;e13 += m.e13;
	e21 += m.e21;e22 += m.e22;e23 += m.e23;
	e31 += m.e31;e32 += m.e32;e33 += m.e33;
	return *this;
}

Matrix3& Matrix3::operator -=(const Matrix3 &m)
{
	e11 -= m.e11;e12 -= m.e12;e13 -= m.e13;
	e21 -= m.e21;e22 -= m.e22;e23 -= m.e23;
	e31 -= m.e31;e32 -= m.e32;e33 -= m.e33;
	return *this;
}

Matrix3& Matrix3::operator *=(const double &s)
{
	e11 *= s;e12 *= s;e13 *= s;
	e21 *= s;e22 *= s;e23 *= s;
	e31 *= s;e32 *= s;e33 *= s;
	return *this;
}

Matrix3& Matrix3::operator /=(const double &s)
{
	e11 /= s;e12 /= s;e13 /= s;
	e21 /= s;e22 /= s;e23 /= s;
	e31 /= s;e32 /= s;e33 /= s;
	return *this;
}

//m1+m2
Matrix3 operator +(const Matrix3 &m1, const Matrix3 &m2)
{
	return Matrix3(m1.e11 + m2.e11,
	               m1.e12 + m2.e12,
	               m1.e13 + m2.e13,
	               m1.e21 + m2.e21,
	               m1.e22 + m2.e22,
	               m1.e23 + m2.e23,
	               m1.e31 + m2.e31,
	               m1.e32 + m2.e32,
	               m1.e33 + m2.e33);
}

//m1-m2
Matrix3 operator -(const Matrix3 &m1, const Matrix3 &m2)
{
	return Matrix3(m1.e11 - m2.e11,				// Zujian Dai fix bug 2019-01-21
	               m1.e12 - m2.e12,
	               m1.e13 - m2.e13,
	               m1.e21 - m2.e21,
	               m1.e22 - m2.e22,
	               m1.e23 - m2.e23,
	               m1.e31 - m2.e31,
	               m1.e32 - m2.e32,
	               m1.e33 - m2.e33);
}

//m/s
Matrix3 operator /(const Matrix3 &m, const double &s)
{
	return Matrix3(m.e11 / s, m.e12 / s, m.e13 / s,
	               m.e21 / s, m.e22 / s, m.e23 / s,
	               m.e31 / s, m.e32 / s, m.e33 / s);
}

//m1*m2
Matrix3 operator *(const Matrix3 &m1, const Matrix3 &m2)
{
	return Matrix3(m1.e11*m2.e11 + m1.e12*m2.e21 + m1.e13*m2.e31,
	               m1.e11*m2.e12 + m1.e12*m2.e22 + m1.e13*m2.e32,
	               m1.e11*m2.e13 + m1.e12*m2.e23 + m1.e13*m2.e33,
	               m1.e21*m2.e11 + m1.e22*m2.e21 + m1.e23*m2.e31,
	               m1.e21*m2.e12 + m1.e22*m2.e22 + m1.e23*m2.e32,
	               m1.e21*m2.e13 + m1.e22*m2.e23 + m1.e23*m2.e33,
	               m1.e31*m2.e11 + m1.e32*m2.e21 + m1.e33*m2.e31,
	               m1.e31*m2.e12 + m1.e32*m2.e22 + m1.e33*m2.e32,
	               m1.e31*m2.e13 + m1.e32*m2.e23 + m1.e33*m2.e33);
}

//m*s
Matrix3 operator *(const Matrix3 &m,const double &s)
{
	return Matrix3(m.e11*s, m.e12*s, m.e13*s,
	               m.e21*s, m.e22*s, m.e23*s,
	               m.e31*s, m.e32*s, m.e33*s);
}

//s*m
Matrix3 operator *(const double &s, const Matrix3 &m)
{
	return Matrix3(m.e11*s, m.e12*s, m.e13*s,
	               m.e21*s, m.e22*s, m.e23*s,
	               m.e31*s, m.e32*s, m.e33*s);
}

// whether m1==m2
bool operator==(const Matrix3 &m1, const Matrix3 &m2)
{
	if(m1.e11 == m2.e11 &&
	   m1.e12 == m2.e12 &&
	   m1.e13 == m2.e13 &&
	   m1.e21 == m2.e21 &&
	   m1.e22 == m2.e22 &&
	   m1.e23 == m2.e23 &&
	   m1.e31 == m2.e31 &&
	   m1.e32 == m2.e32 &&
	   m1.e33 == m2.e33)
	{
		return true;
	}
	return false;
}

//whether m1 != m2
bool operator!=(const Matrix3 &m1, const Matrix3 &m2)
{
    return !(m1 == m2); //!= defined in terms of operator ==
}


void Matrix3::print(void)const
{
	std::cout << e11 << std::setw(15) << e12 << std::setw(15) << e13 << std::endl ;
	std::cout << e21 << std::setw(15) << e22 << std::setw(15) << e23 << std::endl ;
	std::cout << e31 << std::setw(15) << e32 << std::setw(15) << e33 << std::endl ;
	return;
}

ModuleBase::matrix Matrix3::to_matrix(void)const			// Peize Lin add 2021.03.09
{
	ModuleBase::matrix m(3,3);
	m(0,0)=e11;	m(0,1)=e12;	m(0,2)=e13;
	m(1,0)=e21;	m(1,1)=e22;	m(1,2)=e23;
	m(2,0)=e31;	m(2,1)=e32;	m(2,2)=e33;
	return m;
}

}