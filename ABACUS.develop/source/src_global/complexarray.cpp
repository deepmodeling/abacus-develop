/***********************************************************
    DFT++ is a density functional package developed
	by the research group
    of Professor Tomas Arias

    Copyright 1996-2003 Sohrab Ismail-Beigi
************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cstdlib>
#include <cassert>

using namespace std;

#include "complexarray.h"

const complex <double> ZERO(0.0, 0.0), ONE(1.0, 0.0);

void complexArrayxAlloc()
{
	cout << "\n Allocation error for complexArray " << endl;
	abort();
}

ComplexArray::ComplexArray(int bnd1)
{
	ndata = bnd1;
	bound1 = bnd1;
	bound2 = 0;
	init(ndata);
}

ComplexArray::ComplexArray(int bnd1, int bnd2, int bnd3)
{
	ndata = bnd1 * bnd2 * bnd3;
	bound1 = bnd1;
	bound2 = bnd2;
	bound3 = bnd3;
	bound4 = 0;
	init(ndata);
}

ComplexArray::ComplexArray(int bnd1, int bnd2, int bnd3, int bnd4)
{
	ndata = bnd1 * bnd2 * bnd3 * bnd4;
	bound1 = bnd1;
	bound2 = bnd2;
	bound3 = bnd3;
	bound4 = bnd4;
	init(ndata);
}

ComplexArray::~ComplexArray()
{
//    cout << "\n Exit ComplexArray()";
	freemem();
}

void ComplexArray::init(int length)
{
	ndata = length;

//  set_new_handler(complexArrayxAlloc);

	if (ndata > 0)
	{
		ptr = new complex < double> [ndata];	//*sizeof(complex));
		assert(ptr != 0);	// terminate if memory not allocated
	}
	else
		cout << "\n Can't initiliaze ComplexArray.init() with length < 0 ";
}

void ComplexArray::freemem()
{
//	cout << "\n  ComplexArray::freemem() ";
	delete [] ptr;
	ptr = NULL;
	ndata  = 0;
	bound1 = 0;
	bound2 = 0;
	bound3 = 0;
	bound4 = 0;
}

void
ComplexArray::create(int bnd1, int bnd2, int bnd3)
{
//    cout << "\n ComplexArray::create() : " << endl;

	int size1;
	size1 = bnd1 * bnd2 * bnd3;

	if (size1 <= 0)
	{
		cout << "\n Can't create ComplexArray, "
		<< " length = 0 " << size1 << "\n";
		abort();
	}

//    set_new_handler(complexArrayAlloc);

	ndata = size1;

	bound1 = bnd1;

	bound2 = bnd2;

	bound3 = bnd3;

	delete[] ptr;

	ptr = new std::complex < double> [ndata];	//*sizeof(complex));

	assert(ptr != 0);	// terminate if memory not allocated

	for (int i = 0; i < ndata; i++)
		ptr[i] = ZERO;

//  cout << "  size = " << ndata
//	 << "     ptr = " << ptr << endl;
}

void
ComplexArray::create(int bnd1, int bnd2, int bnd3, int bnd4)
{
//    cout << "\n ComplexArray::create() : ";

	int size1;
	size1 = bnd1 * bnd2 * bnd3 * bnd4;

	if (size1 <= 0)
	{
		cout << "\n Can't create ComplexArray, "
		<< " length = 0 " << size1 << "\n";
		abort();
	}

//    set_new_handler(complexArrayAlloc);

	ndata = size1;

	bound1 = bnd1;

	bound2 = bnd2;

	bound3 = bnd3;

	bound4 = bnd4;

	delete[] ptr;

	ptr = new std::complex < double> [ndata];	//*sizeof(complex));

	assert(ptr != 0);	// terminate if memory not allocated

	for (int i = 0; i < ndata; i++)
		ptr[i] = ZERO;

//  cout << "  size = " << ndata
//	 << "     ptr = " << ptr << endl;
}

//////////////////////////////////////
//                                  //
// Operator functions               //
//                                  //
//////////////////////////////////////

/* Assignment:  nonstandard in that it returns void.  To make it standard,
 * replace void -> ComplexArray and uncomment the return *this; */
//#line 200
void ComplexArray::operator=(const ComplexArray & cd)
{
	if (ndata != cd.ndata)
	{
		cout << "\n errore in ComplexArray::operator=(), sizes don't agree.";
		abort();
	}

	for (int i = 0; i < ndata; i++)
		ptr[i] = cd.ptr[i];

	/* return *this; */
}

// Assignment of scalar:  all entries set to c
inline void
ComplexArray::operator=(std::complex < double> c)
{
	for (int i = 0; i < ndata; i++)
		ptr[i] = c;
}

/* Add two ComplexArray */
ComplexArray
ComplexArray::operator+(const ComplexArray &cd)
{

	if (ndata != cd.ndata)
	{
		cout << "\n Size mismatch in ComplexArray::operator+()";
		abort();
	}

	ComplexArray cd2(*this);

	for (int i = 0; i < ndata; i++)
		cd2.ptr[i] += cd.ptr[i];

	return cd2;
}

/* Accumulate sum of ComplexArray */
void
ComplexArray::operator+=(const ComplexArray & cd)
{

	if (ndata != cd.ndata)
	{
		cout << "\n Size mismatch in ComplexArray::operator+=()";
		abort();
	}

	for (int i = 0; i < ndata; i++)
		ptr[i] += cd.ptr[i];
}


/* Subtract two ComplexArray */
ComplexArray
ComplexArray::operator-(const ComplexArray &cd)
{

	if (ndata != cd.ndata)
	{
		cout << "\n Size mismatch in ComplexArray::operator-()";
		abort();
	}

	ComplexArray cd2(*this);

	for (int i = 0; i < ndata; i++)
		cd2.ptr[i] -= cd.ptr[i];

	return cd2;
}

/* Accumulate difference of arrays */
void
ComplexArray::operator-=(const ComplexArray & cd)
{
	if (ndata != cd.ndata)
	{
		cout << "\n Size mismatch in ComplexArray::operator-=()";
		abort();
	}

	for (int i = 0; i < ndata; i++)
		ptr[i] -= cd.ptr[i];

}

/* accumulate pointwise multiply */
void
ComplexArray::operator*=(const ComplexArray & cd)
{
	if (ndata != cd.ndata)
	{
		cout << "\n Size mismatch in ComplexArray::operator*=()";
		abort();
	}

	for (int i = 0; i < ndata; i++)
		ptr[i] *= cd.ptr[i];

}

/* Scale a ComplexArray cd by real r */
ComplexArray
operator*(double r, const ComplexArray &cd)
{

	ComplexArray cd2(cd);

	for (int i = 0; i < cd.ndata; i++)
		cd2.ptr[i] *= r;

	return cd2;
}

/* Scale a ComplexArray by real r */
ComplexArray
ComplexArray::operator*(double r)
{
	ComplexArray cd2(*this);

	for (int i = 0; i < ndata; i++)
		cd2.ptr[i] *= r;

	return cd2;
}

/* Scale a ComplexArray cd by complex number c */
ComplexArray
operator*(std::complex < double> c, const ComplexArray &cd)
{

	ComplexArray cd2(cd.ndata);

//	std::complex <real> *cd2;
//	cd2 = new std::complex <real> [cd.ndata];

	for (int i = 0; i < cd.ndata; i++)
		cd2.ptr[i] = c * cd.ptr[i];

	return cd2;
}

/* Scale a ComplexArray by a complex number c */
ComplexArray
ComplexArray::operator*(std::complex < double> c)
{
	ComplexArray cd(ndata);

	for (int i = 0; i < ndata; i++)
		cd.ptr[i] = ptr[i] * c;

	return cd;
}

/* Scale a ComplexArray by complex c in place */
void
ComplexArray::operator*=(std::complex < double> c)
{

	for (int i = 0; i < ndata; i++)
		ptr[i] *= c;
}

/* Scale a ComplexArray by complex c in place */
void
ComplexArray::operator*=(double r)
{

	for (int i = 0; i < ndata; i++)
		ptr[i] *= r;
}

/////////////////////////////////////////////
//                                         //
//  MEMBER FUNCTIONS:                      //
//                                         //
/////////////////////////////////////////////

// zeroes out the whole array
void ComplexArray::zero_out(void)
{
	for (int i = 0;i < ndata; i++)
		ptr[i] = std::complex < double> (0.0, 0.0);
}

/* Negates all the entries in the array */
void ComplexArray::negate(void)
{

	for (int i = 0;i < ndata; i++)
	{
		ptr[i] = -ptr[i];
//    d[i].real() = -d[i].real();
//    d[i].imag() = -d[i].imag();
	}
}

/* Randomizes w/ uniform distribution in [-.5, .5]*/
void ComplexArray::randomize(void)
{
	for (int i = 0;i < ndata; i++)
	{
		ptr[i] = std::complex < double> (rand() / (RAND_MAX + 1.) - .5,
		                                 rand() / (RAND_MAX + 1.) - .5);
//    d[i].x = rand()/(RAND_MAX+1.) -.5;
//    d[i].y = rand()/(RAND_MAX+1.) -.5;
	}
}

/*
//Write ComplexArray in binary form to the file fname
void ComplexArray::write(char *fname)
{
  FILE *fp;

  fp = dft_fopen(fname,"w");
  dft_fwrite(d,sizeof(scalar),ndata,fp);
  dft_fclose(fp);
}


void ComplexArray::write(FILE *fp)
{
  dft_fwrite(d, sizeof(scalar), ndata, fp);
}

void ComplexArray::writea(char *fname)
{
  FILE *fp;

  fp = dft_fopen(fname,"a");
  dft_fwrite(d,sizeof(scalar),ndata,fp);
  dft_fclose(fp);
}

void ComplexArray::read(char *fname)
{
  FILE *fp;
  fp = dft_fopen(fname,"r");
  dft_fread(d,sizeof(scalar),ndata,fp);
  dft_fclose(fp);


}

void ComplexArray::print()
{
	cout << "\n d[i] : \n ";
  for(int i = 0; i < ndata; i++)
    cout << ptr[i].real() << "+ i"<< ptr[i].imag() << ",";
}
*/
// Sum of absolute squares of all elements in cd
double abs2(const ComplexArray &cd)
{
	double cdcd(0);
	std::complex < double> c;

	for (int i = 0; i < cd.ndata; i++)
	{
		c =  cd.ptr[i];
		cdcd += c.real() * c.real() + c.imag() * c.imag();
	}

	return cdcd;
}

void
add_scale_abs2(const std::complex < double> &c, const ComplexArray & in,
               ComplexArray &out)
{
	if (in.ndata != out.ndata)
		cout << "\n different lengths %d != %d in ComplexArray:"
		<< "\n add_scale_abs2, in.ndata, out.ndata";

	for (int i = 0; i < in.ndata; i++)
	{
		//double z2 = in.ptr[i].real() * in.ptr[i].real() + in.ptr[i].imag() * in.ptr[i].imag();
		out.ptr[i] += std::complex < double> (c.real() * 22, c.imag() * 22);

	}

}

/* Take "dot-product" of two ComplexArray:  sum the diagonals of cd1^cd2 */
std::complex < double>
dot(const ComplexArray &cd1, const ComplexArray &cd2)
{
	//dft_log("ComplexArray::dot()\n");
	int i;
	std::complex < double> dot12;

	if (cd1.ndata != cd2.ndata)
		cout << "\n cd1.ndata != cd2.ndata in dot_ComplexArray()";


	dot12 = 0.0;

	for (i = 0; i < cd1.ndata; i++)
	{
		/* do dot12 += conj(cd1.d[i])*cd2.d[i]; */
		dot12 += std::complex < double>
		         (cd1.ptr[i].real() * cd2.ptr[i].real() +
		          cd1.ptr[i].imag() * cd2.ptr[i].imag(),
		          cd1.ptr[i].real() * cd2.ptr[i].imag() -
		          cd1.ptr[i].imag() * cd2.ptr[i].real());
//    dot12.x += cd1.d[i].x*cd2.d[i].x + cd1.d[i].y*cd2.d[i].y;
//    dot12.y += cd1.d[i].x*cd2.d[i].y - cd1.d[i].y*cd2.d[i].x;
	}

	return dot12;
}

/* Does cd2 += r * cd1 */
void
scale_accumulate(double r, const ComplexArray &cd1, ComplexArray &cd2)
{
	int i;

	for (i = 0; i < cd1.ndata; i++)
		cd2.ptr[i] += r * cd1.ptr[i];
}

/* Does cd2 += c * cd1 */
void
scale_accumulate(std::complex < double> c, const ComplexArray &cd1, ComplexArray &cd2)
{
	int i;

	for (i = 0; i < cd1.ndata; i++)
		cd2.ptr[i] += c * cd1.ptr[i];
}


/* Does cd3 = r1*cd1 + r2*cd2 */
void
scaled_sum(double r1, const ComplexArray &cd1,
           double r2, const ComplexArray &cd2,
           ComplexArray &cd3)
{
	if (cd1.ndata != cd2.ndata || cd1.ndata != cd3.ndata)
		cout << "\nIncompatible sizes in scaled_sum(r1,cd1,r2,cd2,cd3)";

	int i;

	for (i = 0; i < cd1.ndata; i++)
		cd3.ptr[i] = r1 * cd1.ptr[i] + r2 * cd2.ptr[i];
}

/* Does cd3 = c1*cd1 + c2*cd2 */
void
scaled_sum(std::complex < double> c1, const ComplexArray &cd1,
           std::complex < double> c2, const ComplexArray &cd2,
           ComplexArray &cd3)
{
	if (cd1.ndata != cd2.ndata || cd1.ndata != cd3.ndata)
		cout << "\n Incompatible sizes in scaled_sum(c1,cd1,c2,cd2,cd3)";

	int i;

	for (i = 0; i < cd1.ndata; i++)
		cd3.ptr[i] = c1 * cd1.ptr[i] + c2 * cd2.ptr[i];
}


void
point_mult(ComplexArray &in1, ComplexArray &in2, ComplexArray &out)
{
	if (in1.ndata != in2.ndata || in1.ndata != out.ndata)
		cout << "\n Incompatible sizes in point_multiply(in1, in2, out)";

	int i;

	for (i = 0; i < in1.ndata; i++)
	{
		out.ptr[i] = std::complex < double>
		             (in1.ptr[i].real() * in2.ptr[i].real() -
		              in1.ptr[i].imag() * in2.ptr[i].imag(),
		              in1.ptr[i].real() * in2.ptr[i].imag() +
		              in1.ptr[i].imag() * in2.ptr[i].real());
//    out.d[i].x = in1.d[i].x*in2.d[i].x - in1.d[i].y*in2.d[i].y;
//    out.d[i].y = in1.d[i].x*in2.d[i].y + in1.d[i].y*in2.d[i].x;
	}
}

// overloaded subscript operator for const complex Array
// const reference return creates an cvakue
const std::complex < double> &ComplexArray::operator()
(int ind1,int ind2, int ind3)const
{

	// check for subscript out of range error
	if (0 > ind1 || ind1 >= bound1)
	{
		cout << "\n subscript1 out of range error ";
		return ptr[0];
	}
	else if (0 > ind2 || ind2 >= bound2)
	{
		cout << "\n subscript2 out of range error ";
		return ptr[0];
	}
	else if (0 > ind3 || ind3 >= bound3)
	{
		cout << "\n subscript3 out of range error ";
		return ptr[0];
	}

	int ind = (ind1 * bound2 + ind2) * bound3 + ind3 ;

	return ptr[ind];

}

// overloaded subscript operator for const complex Array
// const reference return creates an cvakue
const std::complex < double> &ComplexArray::operator()
(int ind1,int ind2, int ind3, int ind4)const
{

	// check for subscript out of range error
	if (0 > ind1 || ind1 >= bound1)
	{
		cout << "\n subscript1 out of range error ";
		return ptr[0];
	}
	else if (0 > ind2 || ind2 >= bound2)
	{
		cout << "\n subscript2 out of range error ";
		return ptr[0];
	}
	else if (0 > ind3 || ind3 >= bound3)
	{
		cout << "\n subscript3 out of range error ";
		return ptr[0];
	}
	else if (0 > ind4 || ind4 >= bound4)
	{
		cout << "\n subscript4 out of range error ";
		return ptr[0];
	}

	int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;

	return ptr[ind];

}

// overloaded subscript operator for non-const complex Array
// const reference return creates an lvakue
complex < double>
&ComplexArray::operator()(int ind1,int ind2, int ind3)
{

	// check for subscript out of range error
	if (0 > ind1 || ind1 >= bound1)
	{
		cout << "\n subscript1 out of range  ";
		abort();
	}
	else if (0 > ind2 || ind2 >= bound2)
	{
		cout << "\n subscript2 out of range  ";
		abort();
	}
	else if (0 > ind3 || ind3 >= bound3)
	{
		cout << "\n subscript3 out of range  ";
		abort();
	}

	int ind = (ind1 * bound2 + ind2) * bound3 + ind3;

	return ptr[ind];

}

// overloaded subscript operator for non-const complex Array
// const reference return creates an lvakue
complex < double>
&ComplexArray::operator()(int ind1,int ind2, int ind3, int ind4)
{

	// check for subscript out of range error
	if (0 > ind1 || ind1 >= bound1)
	{
		cout << "\n subscript1 out of range  ";
		abort();
	}
	else if (0 > ind2 || ind2 >= bound2)
	{
		cout << "\n subscript2 out of range  ";
		abort();
	}
	else if (0 > ind3 || ind3 >= bound3)
	{
		cout << "\n subscript3 out of range  ";
		abort();
	}
	else if (0 > ind4 || ind4 >= bound4)
	{
		cout << "\n subscript4 out of range  ";
		abort();
	}

	int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;

	return ptr[ind];

}

