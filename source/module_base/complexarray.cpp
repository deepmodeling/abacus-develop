#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <vector>
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
	vec.resize(size, std::complex<double>(0.0, 0.0));
}
void ComplexArray::freemem(){
	vec.resize(0);
	bound1 = 0;
	bound2 = 0;
	bound3 = 0;
	bound4 = 0;
}
void ComplexArray::create(const int bnd1, const int bnd2, const int bnd3, const int bnd4){
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
		vec[i] = cd.vec[i];
	this->bound1 = cd.bound1;
	this->bound2 = cd.bound2;
	this->bound3 = cd.bound3;
	this->bound4 = cd.bound4;
}
ComplexArray::ComplexArray(ComplexArray &&cd){
	this->vec   =cd.vec;	cd.vec.resize(0);
	this->bound1=cd.bound1;	cd.bound1=0;
	this->bound2=cd.bound2;	cd.bound2=0;
	this->bound3=cd.bound3;	cd.bound3=0;
	this->bound4=cd.bound4;	cd.bound4=0;
}
ComplexArray& ComplexArray::operator=(ComplexArray &&cd){
	this->vec   =cd.vec;	cd.vec.resize(0);
	this->bound1=cd.bound1;	cd.bound1=0;
	this->bound2=cd.bound2;	cd.bound2=0;
	this->bound3=cd.bound3;	cd.bound3=0;
	this->bound4=cd.bound4;	cd.bound4=0;
	return *this;}
ComplexArray &ComplexArray::operator=(const ComplexArray & cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		vec[i] = cd.vec[i];
	return *this;}
void ComplexArray::operator=(const std::complex <double> c){
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		vec[i] = c;}
ComplexArray ComplexArray::operator+(const ComplexArray &cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	ComplexArray cd2(*this);
	for (int i = 0; i < size; i++)
		cd2.vec[i] += cd.vec[i];
	return cd2;}
void ComplexArray::operator+=(const ComplexArray & cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		vec[i] += cd.vec[i];
}
ComplexArray ComplexArray::operator-(const ComplexArray &cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	ComplexArray cd2(*this);
	for (int i = 0; i < size; i++)
		cd2.vec[i] -= cd.vec[i];
	return cd2;}
void ComplexArray::operator-=(const ComplexArray & cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		vec[i] -= cd.vec[i];
}
void ComplexArray::operator*=(const ComplexArray & cd){
	const int size = this->getSize();
	assert(size==cd.getSize());
	for (int i = 0; i < size; i++)
		vec[i] *= cd.vec[i];
}
ComplexArray operator*(const double r, const ComplexArray &cd){
	ComplexArray cd2(cd);
	const int size = cd.getSize();
	for (int i = 0; i < size; i++)
		cd2.vec[i] *= r;
	return cd2;}
ComplexArray ComplexArray::operator*(const double r){
	ComplexArray cd2(*this);
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		cd2.vec[i] *= r;
	return cd2;}
ComplexArray operator*(const std::complex < double> c, const ComplexArray &cd){
	const int size = cd.getSize();
	ComplexArray cd2(cd.getSize());
	for (int i = 0; i < size; i++)
		cd2.vec[i] = c * cd.vec[i];
	return cd2;}
ComplexArray ComplexArray::operator*(const std::complex < double> c){
	const int size = this->getSize();
	ComplexArray cd(size);
	for (int i = 0; i < size; i++)
		cd.vec[i] = vec[i] * c;
	return cd;}
void ComplexArray::operator*=(const std::complex <double> c){
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		vec[i] *= c;
}
void ComplexArray::operator*=(const double r){
	const int size = this->getSize();
	for (int i = 0; i < size; i++)
		vec[i] *= r;
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
    	for ( int i = 0;i <size1;++i) {if (this->vec[i] != cd2.vec[i]) {return false;} }
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
    	for ( int i = 0;i <size1;++i) {if (this->vec[i] != cd2.vec[i]) {return true;} }
    	return false;}
void ComplexArray::zero_out(void){
	const int size = this->getSize();
	for (int i = 0;i < size; i++)
		vec[i] = std::complex < double> (0.0, 0.0);
}
void ComplexArray::negate(void){
	const int size = this->getSize();
	for (int i = 0;i < size; i++){
		vec[i] = -vec[i];}
}
void ComplexArray::randomize(void){
	const int size = this->getSize();
	for (int i = 0;i < size; i++){
		vec[i] = std::complex < double> (rand() / (RAND_MAX + 1.) - .5,
		                                 rand() / (RAND_MAX + 1.) - .5);}
}
double abs2(const ComplexArray &cd){
	double cdcd= 0.0;
	const int size = cd.getSize();
	for (int i = 0; i < size; i++){
		const std::complex < double> c = cd.vec[i];
		cdcd += c.real() * c.real() + c.imag() * c.imag();}
	return cdcd;}
// void add_scale_abs2(const std::complex < double> &c, const ComplexArray & in, ComplexArray &out){
// 	assert(in.getSize() == out.getSize());
// 	const int size = in.getSize();
// 	for (int i = 0; i < size; i++)
// 		out.vec[i] += std::complex < double> (c.real() * 22, c.imag() * 22);}
std::complex<double> dot(const ComplexArray &cd1, const ComplexArray &cd2){
	assert(cd1.getSize()==cd2.getSize());
	const int size = cd1.getSize();
	std::complex < double> dot12(0.0,0.0);
	for (int i = 0; i < size; i++){
		dot12 += std::complex < double>
		         (cd1.vec[i].real() * cd2.vec[i].real() +
		          cd1.vec[i].imag() * cd2.vec[i].imag(),
		          cd1.vec[i].real() * cd2.vec[i].imag() -
		          cd1.vec[i].imag() * cd2.vec[i].real());}
	return dot12;}
void scale_accumulate(double r, const ComplexArray &cd1, ComplexArray &cd2){
	assert(cd1.getSize()==cd2.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd2.vec[i] += r * cd1.vec[i];
}
void scale_accumulate(const std::complex<double> c, const ComplexArray &cd1, ComplexArray &cd2){
	assert(cd1.getSize()==cd2.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd2.vec[i] += c * cd1.vec[i];
}
void scaled_sum(double r1, const ComplexArray &cd1,double r2, const ComplexArray &cd2,ComplexArray &cd3){
	assert(cd1.getSize()==cd2.getSize());
	assert(cd1.getSize()==cd3.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd3.vec[i] = r1 * cd1.vec[i] + r2 * cd2.vec[i];
}
void scaled_sum(std::complex < double> c1, const ComplexArray &cd1,std::complex < double> c2, const ComplexArray &cd2,ComplexArray &cd3){
	assert(cd1.getSize()==cd2.getSize());
	assert(cd1.getSize()==cd3.getSize());
	const int size = cd1.getSize();
	for (int i = 0; i < size; i++)
		cd3.vec[i] = c1 * cd1.vec[i] + c2 * cd2.vec[i];}
void point_mult(ComplexArray &in1, ComplexArray &in2, ComplexArray &out){
	assert(in1.getSize()==in2.getSize());
	assert(in1.getSize()==out.getSize());
	const int size = in1.getSize();
	for (int i = 0; i < size; i++){
		out.vec[i] = std::complex < double>
		             (in1.vec[i].real() * in2.vec[i].real() -
		              in1.vec[i].imag() * in2.vec[i].imag(),
		              in1.vec[i].real() * in2.vec[i].imag() +
		              in1.vec[i].imag() * in2.vec[i].real());}
}
const std::complex <double> &ComplexArray::operator()(const int ind1, const int ind2, const int ind3, const int ind4) const{
	assert(ind1>=0);	assert(ind1<bound1);
	assert(ind2>=0);	assert(ind2<bound2);
	assert(ind3>=0);	assert(ind3<bound3);
	assert(ind4>=0);	assert(ind4<bound4);
	const int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;
	return vec[ind];}
std::complex<double>& ComplexArray::operator()(const int ind1,const int ind2,const int ind3,const int ind4){
	assert(ind1>=0);	assert(ind1<bound1);
	assert(ind2>=0);	assert(ind2<bound2);
	assert(ind3>=0);	assert(ind3<bound3);
	assert(ind4>=0);	assert(ind4<bound4);
	const int ind = ((ind1 * bound2 + ind2) * bound3 + ind3) * bound4 + ind4;
	return vec[ind];}
}
