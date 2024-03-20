#ifndef COMPLEX_ARRAY_H
#define COMPLEX_ARRAY_H

#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace ModuleBase
{

/// @brief A basic type of data for complex array
class ComplexArray
{
public:
	std::complex<double> *ptr=nullptr; // data array
	
	ComplexArray(const int bnd1=0, const int bnd2=1, const int bnd3=1, const int bnd4=1);

	~ComplexArray();

	void freemem();

	void create(const int bnd1=0, const int bnd2=1, const int bnd3=1, const int bnd4=1);

	ComplexArray(const ComplexArray &cd);
	ComplexArray(ComplexArray &&cd);

	/****************************************************
 	* OPERATOR FUNCTIONS
 	***************************************************/
	ComplexArray& operator=(ComplexArray &&cd);
	ComplexArray &operator=(const ComplexArray &cd);
	/// Assignment of scalar:  all entries set to c.
	void operator=(std::complex <double> c);
	/// Add two ComplexArray
	ComplexArray operator+(const ComplexArray &cd);
	/// Accumulate sum of ComplexArray
	void operator+=(const ComplexArray &cd);
	/// Subtract two ComplexArray
	ComplexArray operator-(const ComplexArray &cd);
	/// Accumulate difference of arrays
	void operator-=(const ComplexArray &cd);
	/// Scale a ComplexArray by real r
	ComplexArray operator*(const double r);
	/// Scale a ComplexArray by a std::complex number c
	ComplexArray operator*(const std::complex <double> c);
	/// Scale a ComplexArray by real number in place
	void operator*=(const double r);
	/// Scale a ComplexArray by std::complex c in place
	void operator*=(const std::complex <double> c);
	/// accumulate pointwise multiply
	void operator*=(const ComplexArray &cd);
	/// Judge if two ComplexArray is equal
	bool operator== (const ComplexArray &cd2)const;
	/// Judge if two ComplexArray is not equal
	bool operator!= (const ComplexArray &cd2)const;

	/// overloaded subscript operator for non-const std::complex Array const reference return creates an lvakue
	std::complex <double> &operator()
		(const int ind1=0, const int ind2=0, const int ind3=0, const int ind4=0);
	//  std::complex < double> &operator()(int, int, int, int, int);
	/// overloaded subscript operator for const std::complex Array const reference return creates an cvakue
	const std::complex <double> &operator()
		(const int ind1=0, const int ind2=0, const int ind3=0, const int ind4=0)const;
	//  const std::complex < double> &operator()(int, int, int, int, int)const;

	/****************************************************
 	* MEMBER FUNCTIONS
 	***************************************************/
	/// set all elements to be {0.0,0.0}
	void zero_out(void);

	/// Negates all the entries in the array
	void negate(void);

	/// set all elements to a random number whose real/image is between [-0.5,0.5).
	void randomize(void);
	int getBound1()const{ return bound1; }
	int getBound2()const{ return bound2; }
	int getBound3()const{ return bound3; }
	int getBound4()const{ return bound4; }
	int getSize()const{ return bound1*bound2*bound3*bound4; }

private:
	int bound1, bound2, bound3, bound4;	
	void init(const int size);
};
/// Scale a ComplexArray cd by real r
ComplexArray operator*(const double r, const ComplexArray &cd);
/// Scale a ComplexArray cd by std::complex number c
ComplexArray operator*(const std::complex <double> c, const ComplexArray &cd);

/// Sum of absolute squares of all elements in cd
double abs2(const ComplexArray &cd);

// void add_scale_abs2(const std::complex <double> &c, const ComplexArray & in,
//                ComplexArray &out);

/// Take "dot-product" of two ComplexArray:  sum of cd1(conjugate)[i] * cd2[i]
std::complex <double> dot(const ComplexArray &cd1, const ComplexArray &cd2);

/// Does cd2 += r * cd1
void scale_accumulate(double r, const ComplexArray &cd1, ComplexArray &cd2);

/// Does cd2 += c * cd1
void scale_accumulate(std::complex <double> c, const ComplexArray &cd1, ComplexArray &cd2);

/// Does cd3 = r1*cd1 + r2*cd2
void scaled_sum(double r1, const ComplexArray &cd1,
                double r2, const ComplexArray &cd2,
                ComplexArray &cd3);

/// Does cd3 = c1*cd1 + c2*cd2
void scaled_sum(std::complex <double> c1, const ComplexArray &cd1,
           std::complex <double> c2, const ComplexArray &cd2,
           ComplexArray &cd3);

/// out[i] = a1[i] * in2[i]
void point_mult(ComplexArray &a1, ComplexArray &in2, ComplexArray &out);

/// set elements of u as zero which u is 1_d std::complex array
template <class T>
void zeros(std::complex <T> *u, int n)
{
	if (n == 0 || u == 0)
	{
		std::cout << "\n error in zeros(),n or u = 0";
		return;
	}

	for (int i = 0;i < n;i++)
	{
		u[i] = std::complex <T> (0.0, 0.0);
	}
}
}

#endif // COMPLEX_ARRAY_H
