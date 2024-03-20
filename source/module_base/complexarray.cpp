#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cstdlib>
#include <cassert>

#include "complexarray.h"
#include "global_function.h"
namespace ModuleBase
{
void complexArrayxAlloc(){ModuleBase::WARNING_QUIT("ComplexArray","Allocation error for complexArray");}

ComplexArray::ComplexArray(const int bnd1, const int bnd2, const int bnd3, const int bnd4){
	bound1 = bnd1;
	bound2 = bnd2;
	bound3 = bnd3;
	bound4 = bnd4;
	init(this->getSize());
}

ComplexArray::~ComplexArray(){
	freemem();
}
void ComplexArray::init(const int size){
	assert(size>=0);
	if(size>0){
		ptr = new std::complex<double> [size];
		assert(ptr != 0);}
	else
	{ptr = nullptr;}
}
void ComplexArray::freemem(){
	delete [] ptr;
	ptr = nullptr;
	bound1 = 0;
	bound2 = 0;
	bound3 = 0;
	bound4 = 0;
}
void ComplexArray::create(const int bnd1, const int bnd2, const int bnd3, const int bnd4){
	delete [] ptr;
	bound1 = bnd1;
	bound2 = bnd2;
	bound3 = bnd3;
	bound4 = bnd4;
	const int size = this->getSize();
	this->init(size);
	this->zero_out();
}
ComplexArray::ComplexArray(const ComplexArray &cd){
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
ComplexArray::ComplexArray(ComplexArray &&cd){
	delete [] this->ptr;
	this->ptr   =cd.ptr;	cd.ptr   =nullptr;
	this->bound1=cd.bound1;	cd.bound1=0;
	this->bound2=cd.bound2;	cd.bound2=0;
	this->bound3=cd.bound3;	cd.bound3=0;
	this->bound4=cd.bound4;	cd.bound4=0;
}
ComplexArray& ComplexArray::operator=(ComplexArray &&cd){
	delete [] this->ptr;
	this->ptr   =cd.ptr;	cd.ptr   =nullptr;
	this->bound1=cd.bound1;	cd.bound1=0;
	this->bound2=cd.bound2;	cd.bound2=0;
	this->bound3=cd.bound3;	cd.bound3=0;
	this->bound4=cd.bound4;	cd.bound4=0;
	return *this;}
ComplexArray &ComplexArray::operator=(const ComplexArray & cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		ptr[i] = cd.ptr[i];
	return *this;}
void ComplexArray::operator=(const std::complex < double> c){
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		ptr[i] = c;}
ComplexArray ComplexArray::operator+(const ComplexArray &cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	ComplexArray cd2(*this);
	for (int i = 0; i < size; i++)
		cd2.ptr[i] += cd.ptr[i];
	return cd2;}
void ComplexArray::operator+=(const ComplexArray & cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		ptr[i] += cd.ptr[i];
}
ComplexArray ComplexArray::operator-(const ComplexArray &cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	ComplexArray cd2(*this);
	for (int i = 0; i < size; i++)
		cd2.ptr[i] -= cd.ptr[i];
	return cd2;}
void ComplexArray::operator-=(const ComplexArray & cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		ptr[i] -= cd.ptr[i];
}
void ComplexArray::operator*=(const ComplexArray & cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		ptr[i] *= cd.ptr[i];
}
ComplexArray operator*(const double r, const ComplexArray &cd){
	ComplexArray cd2(cd);
	const int size = cd.getSize();
	for (int i = 0; i < size; i++)
		cd2.ptr[i] *= r;
	return cd2;}
ComplexArray ComplexArray::operator*(const double r){
	ComplexArray cd2(*this);
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		cd2.ptr[i] *= r;
	return cd2;}
ComplexArray operator*(const std::complex < double> c, const ComplexArray &cd){
	const int size = cd.getSize();
	ComplexArray cd2(cd.getSize());
	for (int i = 0; i < size; i++)
		cd2.ptr[i] = c * cd.ptr[i];
	return cd2;}
ComplexArray ComplexArray::operator*(const std::complex < double> c){
	const int size = this->getSize();
	ComplexArray cd(size);
	for (int i = 0; i < size; i++)
		cd.ptr[i] = ptr[i] * c;
	return cd;}
void ComplexArray::operator*=(const std::complex <double> c){
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		ptr[i] *= c;
}
void ComplexArray::operator*=(const double r){
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		ptr[i] *= r;
}
bool ComplexArray::operator==(const ComplexArray &cd2)const{
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
    	return true;}
bool ComplexArray::operator!=(const ComplexArray &cd2)const{
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
    	return false;}
void ComplexArray::zero_out(void){
	const int size = this->getSize();
	for (int i = 0;i < size; i++)
		ptr[i] = std::complex < double> (0.0, 0.0);
}
void ComplexArray::negate(void){
	const int size = this->getSize();
	for (int i = 0;i < size; i++){
		ptr[i] = -ptr[i];}
}
void ComplexArray::randomize(void){
	const int size = this->getSize();
	for (int i = 0;i < size; i++){
		ptr[i] = std::complex < double> (rand() / (RAND_MAX + 1.) - .5,
		                                 rand() / (RAND_MAX + 1.) - .5);}
}
double abs2(const ComplexArray &cd){
	double cdcd= 0.0;
	const int size = cd.getSize();
	for (int i = 0; i < size; i++){
		const std::complex < double> c = cd.ptr[i];
		cdcd += c.real() * c.real() + c.imag() * c.imag();}
	return cdcd;}
// void add_scale_abs2(const std::complex < double> &c, const ComplexArray & in, ComplexArray &out){
// 	assert(in.getSize() == out.getSize());
// 	const int size = in.getSize();
// 	for (int i = 0; i < size; i++)
// 		out.ptr[i] += std::complex < double> (c.real() * 22, c.imag() * 22);}
std::complex<double> dot(const ComplexArray &cd1, const ComplexArray &cd2){
	assert(cd1.getSize()==cd2.getSize());
	const int size = cd1.getSize();
	std::complex < double> dot12(0.0,0.0);
	for (int i = 0; i < size; i++){
		dot12 += std::complex < double>
		         (cd1.ptr[i].real() * cd2.ptr[i].real() +
		          cd1.ptr[i].imag() * cd2.ptr[i].imag(),
		          cd1.ptr[i].real() * cd2.ptr[i].imag() -
		          cd1.ptr[i].imag() * cd2.ptr[i].real());}
	return dot12;}
void scale_accumulate(double r, const ComplexArray &cd1, ComplexArray &cd2){
	assert(cd1.getSize()==cd2.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd2.ptr[i] += r * cd1.ptr[i];
}
void scale_accumulate(const std::complex<double> c, const ComplexArray &cd1, ComplexArray &cd2){
	assert(cd1.getSize()==cd2.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd2.ptr[i] += c * cd1.ptr[i];
}
void scaled_sum(double r1, const ComplexArray &cd1,double r2, const ComplexArray &cd2,ComplexArray &cd3){
	assert(cd1.getSize()==cd2.getSize());
	assert(cd1.getSize()==cd3.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd3.ptr[i] = r1 * cd1.ptr[i] + r2 * cd2.ptr[i];
}
void scaled_sum(std::complex < double> c1, const ComplexArray &cd1,std::complex < double> c2, const ComplexArray &cd2,ComplexArray &cd3){
	assert(cd1.getSize()==cd2.getSize());
	assert(cd1.getSize()==cd3.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd3.ptr[i] = c1 * cd1.ptr[i] + c2 * cd2.ptr[i];}
void point_mult(ComplexArray &in1, ComplexArray &in2, ComplexArray &out){
	assert(in1.getSize()==in2.getSize());
	assert(in1.getSize()==out.getSize());
	const int size = in1.getSize();
	for (int i = 0; i < size; i++){
		out.ptr[i] = std::complex < double>
		             (in1.ptr[i].real() * in2.ptr[i].real() -
		              in1.ptr[i].imag() * in2.ptr[i].imag(),
		              in1.ptr[i].real() * in2.ptr[i].imag() +
		              in1.ptr[i].imag() * in2.ptr[i].real());}
}
const std::complex <double> &ComplexArray::operator()(const int ind1, const int ind2, const int ind3, const int ind4) const{
	assert(ind1>=0);	assert(ind1<bound1);
	assert(ind2>=0);	assert(ind2<bound2);
	assert(ind3>=0);	assert(ind3<bound3);
	assert(ind4>=0);	assert(ind4<bound4);
	const int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;
	return ptr[ind];}
std::complex<double>& ComplexArray::operator()(const int ind1,const int ind2,const int ind3,const int ind4){
	assert(ind1>=0);	assert(ind1<bound1);
	assert(ind2>=0);	assert(ind2<bound2);
	assert(ind3>=0);	assert(ind3<bound3);
	assert(ind4>=0);	assert(ind4<bound4);
	const int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;
	return ptr[ind];}
}
