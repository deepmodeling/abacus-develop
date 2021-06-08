#ifndef STO_CHEBYCHEV_H
#define STO_CHEBYCHEV_H

#include "tools.h"
#include "../src_global/mymath.h"
//----------------------------------------------
// Chebychev Filtering
//----------------------------------------------

class Stochastic_Chebychev
{

	public:

    // constructor and deconstructor
    Stochastic_Chebychev();
    ~Stochastic_Chebychev();

    void init(int &dim, int& chetype);
    
    void calcoef(double fun(double));

    complex<double> sumallterms();
    
    void calfinalvec(
		void fun(complex<double> *in, complex<double> *out, const int), 
		complex<double> *wavein, 
		complex<double> *waveout, 
		const int m = 1);

    bool checkconverge(
		void tfun(complex<double> *in, complex<double> *out, const int), 
 		complex<double> *wavein,
		double& tmax, 
		double& tmin, 
		double stept);

    void calpolyval(
		void fun(complex<double> *in, complex<double> *out, const int), 
		complex<double> *wavein, 
		const int m =1);

    int norder;
    int extend;
    int norder2;  // 2 * norder

    double* coef;  

	// expansion coefficient of each order
	// only first norder coefficients are usefull
    double* dcoef; //[norder2]

	// temporary complex expansion coefficient of each order
	// only first norder coefficients are usefull
    fftw_complex *ccoef;  //[norder2]

    double *polyvalue;

    complex<double> *vecn;

    fftw_plan plancoef;

    bool initplan;
	bool initcoef;
	bool getcoef;
	bool getpolyval;

    private:

    int ndim; //dim of vector

    template<class T>
    void recurs(
		T *arraynp1, 
		T* arrayn, 
		T *arrayn_1, 
		void fun(T *in,T *out, const int),
		const int m);

};


template<class T>
void Stochastic_Chebychev::recurs(
		T *arraynp1, 
		T* arrayn, 
		T *arrayn_1, 
		void fun(T *in,T *out, const int),
		const int m)
{
    fun(arrayn,arraynp1,m);
    for(int i = 0; i < ndim * m; ++i)
    {
        arraynp1[i]=2*arraynp1[i]-arrayn_1[i];
    }
}
#endif// Eelectrons_Chebychev
