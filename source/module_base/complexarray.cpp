#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cstdlib>
#include <cassert>



#include "complexarray.h"
namespace ModuleBase
{

void complexArrayxAlloc()
{
	std::cout << "\n Allocation error for complexArray " << std::endl;
	abort();
}

ComplexArray::ComplexArray(const int bnd1, const int bnd2, const int bnd3, const int bnd4)
{
	bound1 = bnd1;
	bound2 = bnd2;
	bound3 = bnd3;
	bound4 = bnd4;
	init(this->getSize());
}

ComplexArray::~ComplexArray()
{
//    std::cout << "\n Exit ComplexArray()";
	freemem();
}

void ComplexArray::init(const int size)
{
	assert(size>=0);
//  set_new_handler(complexArrayxAlloc);
	if(size>0)
	{
		ptr = new std::complex<double> [size];	//*sizeof(std::complex));
		assert(ptr != 0);	// terminate if memory not allocated
	}
	else
	{
		ptr = nullptr;
	}
}

void ComplexArray::freemem()
{
//	std::cout << "\n  ComplexArray::freemem() ";
	delete [] ptr;
	ptr = nullptr;
	bound1 = 0;
	bound2 = 0;
	bound3 = 0;
	bound4 = 0;
}

void ComplexArray::create(const int bnd1, const int bnd2, const int bnd3, const int bnd4)
{
	delete [] ptr;

	bound1 = bnd1;
	bound2 = bnd2;
	bound3 = bnd3;
	bound4 = bnd4;
	const int size = this->getSize();
	this->init(size);
	this->zero_out();
}

ComplexArray::ComplexArray(const ComplexArray &cd)
{
	this->freemem();

	const int size = cd.getSize();
	this->init(size);
	for (int i = 0; i < size; i++)
		ptr[i] = cd.ptr[i];
		
	this->bound1 = cd.bound1;
	this->bound2 = cd.bound2;
	this->bound3 = cd.bound3;
	this->bound4 = cd.bound4;
}

ComplexArray::ComplexArray(ComplexArray &&cd)
{
	delete [] this->ptr;
	this->ptr   =cd.ptr;	cd.ptr   =nullptr;
	this->bound1=cd.bound1;	cd.bound1=0;
	this->bound2=cd.bound2;	cd.bound2=0;
	this->bound3=cd.bound3;	cd.bound3=0;
	this->bound4=cd.bound4;	cd.bound4=0;
}

ComplexArray& ComplexArray::operator=(ComplexArray &&cd)
{
	delete [] this->ptr;
	this->ptr   =cd.ptr;	cd.ptr   =nullptr;
	this->bound1=cd.bound1;	cd.bound1=0;
	this->bound2=cd.bound2;	cd.bound2=0;
	this->bound3=cd.bound3;	cd.bound3=0;
	this->bound4=cd.bound4;	cd.bound4=0;
	return *this;
}

//////////////////////////////////////
//                                  //
// Operator functions               //
//                                  //
//////////////////////////////////////

/* Assignment:  nonstandard in that it returns void.  To make it standard,
 * replace void -> ComplexArray and uncomment the return *this; */
//#line 200
ComplexArray &ComplexArray::operator=(const ComplexArray & cd)
{
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		ptr[i] = cd.ptr[i];
	return *this;
}

// Assignment of scalar:  all entries set to c
void ComplexArray::operator=(const std::complex < double> c)
{
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		ptr[i] = c;
}

/* Add two ComplexArray */
ComplexArray ComplexArray::operator+(const ComplexArray &cd)
{
	const int size = this->getSize();
	assert(size==cd.getSize());
	ComplexArray cd2(*this);
	for (int i = 0; i < size; i++)
		cd2.ptr[i] += cd.ptr[i];
	return cd2;
}

/* Accumulate sum of ComplexArray */
void ComplexArray::operator+=(const ComplexArray & cd)
{
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		ptr[i] += cd.ptr[i];
}


/* Subtract two ComplexArray */
ComplexArray ComplexArray::operator-(const ComplexArray &cd)
{
	const int size = this->getSize();
	assert(size==cd.getSize());
	ComplexArray cd2(*this);
	for (int i = 0; i < size; i++)
		cd2.ptr[i] -= cd.ptr[i];
	return cd2;
}

/* Accumulate difference of arrays */
void ComplexArray::operator-=(const ComplexArray & cd)
{
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		ptr[i] -= cd.ptr[i];
}

/* accumulate pointwise multiply */
void ComplexArray::operator*=(const ComplexArray & cd)
{
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		ptr[i] *= cd.ptr[i];
}

/* Scale a ComplexArray cd by real r */
ComplexArray operator*(const double r, const ComplexArray &cd)
{
	ComplexArray cd2(cd);
	const int size = cd.getSize();
	for (int i = 0; i < size; i++)
		cd2.ptr[i] *= r;
	return cd2;
}

/* Scale a ComplexArray by real r */
ComplexArray ComplexArray::operator*(const double r)
{
	ComplexArray cd2(*this);
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		cd2.ptr[i] *= r;
	return cd2;
}

/* Scale a ComplexArray cd by std::complex number c */
ComplexArray operator*(const std::complex < double> c, const ComplexArray &cd)
{
	const int size = cd.getSize();
	ComplexArray cd2(cd.getSize());
	for (int i = 0; i < size; i++)
		cd2.ptr[i] = c * cd.ptr[i];
	return cd2;
}

/* Scale a ComplexArray by a std::complex number c */
ComplexArray ComplexArray::operator*(std::complex < double> c)
{
	const int size = this->getSize();
	ComplexArray cd(size);
	for (int i = 0; i < size; i++)
		cd.ptr[i] = ptr[i] * c;
	return cd;
}

/* Scale a ComplexArray by std::complex c in place */
void ComplexArray::operator*=(std::complex < double> c)
{
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		ptr[i] *= c;
}

/* Scale a ComplexArray by std::complex c in place */
void ComplexArray::operator*=(double r)
{
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		ptr[i] *= r;
}

/* Judge if two ComplexArray is equal */
bool ComplexArray::operator==(const ComplexArray &cd2)const 
{
	const int size1 = this->getSize();
	const int size2 = cd2.getSize();
	const int b11 = this->getBound1();
	const int b12 = this->getBound2();
	const int b13 = this->getBound3();
	const int b14 = this->getBound4();
	const int b21 = cd2.getBound1();
	const int b22 = cd2.getBound2();
	const int b23 = cd2.getBound3();
	const int b24 = cd2.getBound4();
	if (size1 != size2) {return false;}
	if (b11 != b21) {return false;}
    if (b12 != b22) {return false;}
    if (b13 != b23) {return false;}
    if (b14 != b24) {return false;}
    for ( int i = 0;i <size1;++i) {if (this->ptr[i] != cd2.ptr[i]) {return false;} }
    return true;
}

/* Judge if two ComplexArray is not equal */
bool ComplexArray::operator!=(const ComplexArray &cd2)const 
{
	const int size1 = this->getSize();
	const int size2 = cd2.getSize();
	const int b11 = this->getBound1();
	const int b12 = this->getBound2();
	const int b13 = this->getBound3();
	const int b14 = this->getBound4();
	const int b21 = cd2.getBound1();
	const int b22 = cd2.getBound2();
	const int b23 = cd2.getBound3();
	const int b24 = cd2.getBound4();
	if (size1 != size2) {return true;}
	if (b11 != b21) {return true;}
    if (b12 != b22) {return true;}
    if (b13 != b23) {return true;}
    if (b14 != b24) {return true;}
    for ( int i = 0;i <size1;++i) {if (this->ptr[i] != cd2.ptr[i]) {return true;} }
    return false;
}

/////////////////////////////////////////////
//                                         //
//  MEMBER FUNCTIONS:                      //
//                                         //
/////////////////////////////////////////////

// zeroes out the whole array
void ComplexArray::zero_out(void)
{
	const int size = this->getSize();
	for (int i = 0;i < size; i++)
		ptr[i] = std::complex < double> (0.0, 0.0);
}

/* Negates all the entries in the array */
void ComplexArray::negate(void)
{
	const int size = this->getSize();
	for (int i = 0;i < size; i++)
	{
		ptr[i] = -ptr[i];
//    d[i].real() = -d[i].real();
//    d[i].imag() = -d[i].imag();
	}
}

/* Randomizes w/ uniform distribution in [-.5, .5]*/
void ComplexArray::randomize(void)
{
	const int size = this->getSize();
	for (int i = 0;i < size; i++)
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
  dft_fwrite(d,sizeof(scalar),size,fp);
  dft_fclose(fp);
}


void ComplexArray::write(FILE *fp)
{
  dft_fwrite(d, sizeof(scalar), size, fp);
}

void ComplexArray::writea(char *fname)
{
  FILE *fp;

  fp = dft_fopen(fname,"a");
  dft_fwrite(d,sizeof(scalar),size,fp);
  dft_fclose(fp);
}

void ComplexArray::read(char *fname)
{
  FILE *fp;
  fp = dft_fopen(fname,"r");
  dft_fread(d,sizeof(scalar),size,fp);
  dft_fclose(fp);


}

void ComplexArray::print()
{
	std::cout << "\n d[i] : \n ";
  for(int i = 0; i < size; i++)
    std::cout << ptr[i].real() << "+ i"<< ptr[i].imag() << ",";
}
*/
// Sum of absolute squares of all elements in cd
double abs2(const ComplexArray &cd)
{
	double cdcd= 0.0;
	const int size = cd.getSize();
	for (int i = 0; i < size; i++)
	{
		const std::complex < double> c = cd.ptr[i];
		cdcd += c.real() * c.real() + c.imag() * c.imag();
	}
	return cdcd;
}

void add_scale_abs2(const std::complex < double> &c, const ComplexArray & in,
               ComplexArray &out)
{
	assert(in.getSize() == out.getSize());
	const int size = in.getSize();
	for (int i = 0; i < size; i++)
	{
		//double z2 = in.ptr[i].real() * in.ptr[i].real() + in.ptr[i].imag() * in.ptr[i].imag();
		out.ptr[i] += std::complex < double> (c.real() * 22, c.imag() * 22);
	}
}

/* Take "dot-product" of two ComplexArray:  sum the diagonals of cd1^cd2 */
std::complex<double> dot(const ComplexArray &cd1, const ComplexArray &cd2)
{
	assert(cd1.getSize()==cd2.getSize());
	const int size = cd1.getSize();
	//dft_log("ComplexArray::dot()\n");
	std::complex < double> dot12(0.0,0.0);
	for (int i = 0; i < size; i++)
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
void scale_accumulate(double r, const ComplexArray &cd1, ComplexArray &cd2)
{
	assert(cd1.getSize()==cd2.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd2.ptr[i] += r * cd1.ptr[i];
}

/* Does cd2 += c * cd1 */
void scale_accumulate(const std::complex<double> c, const ComplexArray &cd1, ComplexArray &cd2)
{
	assert(cd1.getSize()==cd2.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd2.ptr[i] += c * cd1.ptr[i];
}


/* Does cd3 = r1*cd1 + r2*cd2 */
void scaled_sum(double r1, const ComplexArray &cd1,
           double r2, const ComplexArray &cd2,
           ComplexArray &cd3)
{
	assert(cd1.getSize()==cd2.getSize());
	assert(cd1.getSize()==cd3.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd3.ptr[i] = r1 * cd1.ptr[i] + r2 * cd2.ptr[i];
}

/* Does cd3 = c1*cd1 + c2*cd2 */
void scaled_sum(std::complex < double> c1, const ComplexArray &cd1,
           std::complex < double> c2, const ComplexArray &cd2,
           ComplexArray &cd3)
{
	assert(cd1.getSize()==cd2.getSize());
	assert(cd1.getSize()==cd3.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd3.ptr[i] = c1 * cd1.ptr[i] + c2 * cd2.ptr[i];
}


void point_mult(ComplexArray &in1, ComplexArray &in2, ComplexArray &out)
{
	assert(in1.getSize()==in2.getSize());
	assert(in1.getSize()==out.getSize());
	const int size = in1.getSize();
	for (int i = 0; i < size; i++)
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

// overloaded subscript operator for const std::complex Array
// const reference return creates an cvakue
const std::complex < double> &ComplexArray::operator()
	(const int ind1, const int ind2, const int ind3, const int ind4) const
{
	assert(ind1>=0);	assert(ind1<bound1);
	assert(ind2>=0);	assert(ind2<bound2);
	assert(ind3>=0);	assert(ind3<bound3);
	assert(ind4>=0);	assert(ind4<bound4);
	const int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;
	return ptr[ind];
}

// overloaded subscript operator for non-const std::complex Array
// const reference return creates an lvakue
std::complex<double>& ComplexArray::operator()
	(const int ind1, const int ind2, const int ind3, const int ind4)
{
	assert(ind1>=0);	assert(ind1<bound1);
	assert(ind2>=0);	assert(ind2<bound2);
	assert(ind3>=0);	assert(ind3<bound3);
	assert(ind4>=0);	assert(ind4<bound4);
	const int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;
	return ptr[ind];
}

}