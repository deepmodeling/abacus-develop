#ifndef MYMATH_H
#define MYMATH_H
#if defined __FFTW2
#include "fftw.h"
#elif defined __FFTW3
#include "fftw3.h"
#endif
namespace ModuleBase
{
// in mymath1.cpp
#ifdef __FFTW3
void fftw_zeros(fftw_complex *data,int n);//not be used now!
#endif

double rndm();//not be used now!
void simpson(const int mesh,const double *func,const double *rab, double &asum);//not be used now!

double my_erf(double x);//not be used now!
double my_erfc(double x);//not be used now!
double gauss_freq(double x);//not be used now!

void heapsort(int n, double *r, int *ind);
void heapAjust(double r[], int ind[], int s, int m);//not be used now!
void hpsort(int n, double *ra, int *ind);

}

#endif // MYMATH_H


