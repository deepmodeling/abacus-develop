#include "./constants.h"
#include "./tool_quit.h"
#include "assert.h"
#include "./global_function.h"
namespace ModuleBase
{

template<typename REAL>
template<class T>
void Chebyshev<REAL>:: calcoef_real(T *ptr, REAL (T::*fun)(REAL))
{
    std::complex<REAL> *pcoef = (std::complex<REAL> *)this->fftw.ccoef;

    //three point = 2/3 M + 1/3 T;
    //-----------------------------------------------
    //(M)iddle point integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)((REAL) cos((i+0.5)*ModuleBase::TWO_PI/norder2));
    }

    //this->fftw.dcoef --FFT--> fftw.pcoef
    this->fftw.execute_fftw();

    for(int i = 0; i<norder; ++i)
    {
        REAL phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef_real[i] = (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3;
        }
        else
        {
            coef_real[i] = (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3;
        }
    }

    //-----------------------------------------------
    //(T)rapezoid integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)(cos(i*ModuleBase::TWO_PI/norder2));
    }

    //this->fftw.dcoef --FFT--> fftw.pcoef
    this->fftw.execute_fftw();

    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef_real[i] += real(pcoef[i]) / norder2 * 1 / 3;
        }
        else
        {
            coef_real[i] += real(pcoef[i]) / norder2 * 2 / 3;
        }
    } 

    getcoef_real = true;
	return;
}

template<typename REAL>
template<class T>
void Chebyshev<REAL>::calcoef_complex(T *ptr, std::complex<REAL> (T::*fun)(std::complex<REAL>))
{
     std::complex<REAL> *pcoef = ( std::complex<REAL> *)this->fftw.ccoef;

    //three point = 2/3 M + 1/3 T;
    //-----------------------------------------------
    //(M)iddle point integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)(cos((i+0.5)*ModuleBase::TWO_PI/norder2)).real();
    }
    this->fftw.execute_fftw();
    for(int i = 0; i<norder; ++i)
    {
        REAL phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef_complex[i].real( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3 );
        }
        else
        {
            coef_complex[i].real( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3 );
        }
    }

    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)(cos((i+0.5)*ModuleBase::TWO_PI/norder2)).imag();
    }
    this->fftw.execute_fftw();
    for(int i = 0; i<norder; ++i)
    {
        REAL phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef_complex[i].imag( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3 );
        }
        else
        {
            coef_complex[i].imag( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3 );
        }
    }

    //-----------------------------------------------
    //(T)rapezoid integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)(cos(i*ModuleBase::TWO_PI/norder2)).real();
    }
    this->fftw.execute_fftw();
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef_complex[i].real( real(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3 );
        }
        else
        {
            coef_complex[i].real( real(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3 );
        }
    } 

    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun)(cos(i*ModuleBase::TWO_PI/norder2)).imag();
    }
    this->fftw.execute_fftw();
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef_complex[i].imag( imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3 );
        }
        else
        {
            coef_complex[i].imag( imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3 );
        }
    } 

    getcoef_complex = true;
	return;
}

template<typename REAL>
template<class T>
void Chebyshev<REAL>::calcoef_pair(T *ptr, REAL (T::*fun1)(REAL), REAL (T::*fun2)(REAL))
{
     std::complex<REAL> *pcoef = ( std::complex<REAL> *)this->fftw.ccoef;

    //three point = 2/3 M + 1/3 T;
    //-----------------------------------------------
    //(M)iddle point integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun1)(cos((i+0.5)*ModuleBase::TWO_PI/norder2));
    }
    this->fftw.execute_fftw();
    for(int i = 0; i<norder; ++i)
    {
        REAL phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef_complex[i].real( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3 );
        }
        else
        {
            coef_complex[i].real( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3 );
        }
    }

    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun2)(cos((i+0.5)*ModuleBase::TWO_PI/norder2));
    }
    this->fftw.execute_fftw();
    for(int i = 0; i<norder; ++i)
    {
        REAL phi=i*ModuleBase::PI/norder2;
        if(i == 0)
        {
            coef_complex[i].imag( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3 );
        }
        else
        {
            coef_complex[i].imag( (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3 );
        }
    }

    //-----------------------------------------------
    //(T)rapezoid integral method part
    //-----------------------------------------------
    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun1)(cos(i*ModuleBase::TWO_PI/norder2));
    }
    this->fftw.execute_fftw();
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef_complex[i].real( real(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3 );
        }
        else
        {
            coef_complex[i].real( real(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3 );
        }
    } 

    for(int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i]=(ptr->*fun2)(cos(i*ModuleBase::TWO_PI/norder2));
    }
    this->fftw.execute_fftw();
    for(int i = 0; i<norder; ++i)
    {
        if(i == 0)
        {
            coef_complex[i].imag( imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3 );
        }
        else
        {
            coef_complex[i].imag( imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3 );
        }
    } 

    getcoef_complex = true;
	return;
}

template<typename REAL>
template<class T>
void Chebyshev<REAL>::calfinalvec_real(T *ptr, 
	void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
	std::complex<REAL> *wavein, 
	std::complex<REAL> *waveout, 
	const int N, const int LDA, const int m)
{
    if(!getcoef_real) ModuleBase::WARNING_QUIT("Chebyshev<REAL>", "Please calculate coef_real first!");

    std::complex<REAL> *arraynp1;
	std::complex<REAL> *arrayn;
	std::complex<REAL> *arrayn_1;
    assert(N>=0 && LDA >= N);
    int ndmxt;
    if(m == 1) ndmxt = N * m;
    else       ndmxt = LDA * m; 

    arraynp1 = new std::complex<REAL> [ndmxt];
    arrayn = new std::complex<REAL> [ndmxt];
    arrayn_1 = new std::complex<REAL> [ndmxt];
  
    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);
    
    (ptr->*funA)(arrayn_1, arrayn, m);
    
    //0- & 1-st order
    for(int i = 0; i < ndmxt; ++i)
    {
        waveout[i] = coef_real[0] * arrayn_1[i] + coef_real[1] * arrayn[i];
    }    

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs_complex(ptr, funA, arraynp1, arrayn, arrayn_1, N, LDA, m);
        for(int i = 0; i < ndmxt; ++i)
        {
            waveout[i] += coef_real[ior] * arraynp1[i];
        }
        std::complex<REAL> * tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return;
}

template<typename REAL>
template<class T>
void Chebyshev<REAL>::calfinalvec_complex(T *ptr, 
	void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
	std::complex<REAL> *wavein, 
	std::complex<REAL> *waveout, 
	const int N, const int LDA, const int m)
{
    if(!getcoef_complex) ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef_complex first!");

    std::complex<REAL> *arraynp1;
	std::complex<REAL> *arrayn;
	std::complex<REAL> *arrayn_1;
    assert(N>=0 && LDA >= N);
    int ndmxt;
    if(m == 1) ndmxt = N * m;
    else       ndmxt = LDA * m;

    arraynp1 = new std::complex<REAL> [ndmxt];
    arrayn = new std::complex<REAL> [ndmxt];
    arrayn_1 = new std::complex<REAL> [ndmxt];
  
    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);
    
    (ptr->*funA)(arrayn_1, arrayn,m);
    
    //0- & 1-st order
    for(int i = 0; i < ndmxt; ++i)
    {
        waveout[i] = coef_complex[0] * arrayn_1[i] + coef_complex[1] * arrayn[i];
    }    

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs_complex(ptr, funA, arraynp1, arrayn, arrayn_1, N, LDA, m);
        for(int i = 0; i < ndmxt; ++i)
        {
            waveout[i] += coef_complex[ior] * arraynp1[i];
        }
        std::complex<REAL> * tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return;
}

template<typename REAL>
template<class T>
void Chebyshev<REAL>::calpolyvec_complex(T *ptr, 
	void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
	std::complex<REAL> *wavein, 
	std::complex<REAL> *polywaveout, 
	const int N, const int LDA, const int m)
{

    assert(N>=0 && LDA >= N);
    const int ndmxt = LDA * m; 

    std::complex<REAL> *arraynp1 = polywaveout + 2 * ndmxt;
	std::complex<REAL> *arrayn = polywaveout + ndmxt;
	std::complex<REAL> *arrayn_1 = polywaveout;
    
    std::complex<REAL> *tmpin = wavein, *tmpout = arrayn_1;
    for(int i = 0 ; i < m ; ++i)
    {
        ModuleBase::GlobalFunc::DCOPY(tmpin, tmpout, N);
        tmpin += LDA;
        tmpout += LDA;
    }
    
    //1-st order
    (ptr->*funA)(arrayn_1, arrayn, m);

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs_complex(ptr, funA, arraynp1, arrayn, arrayn_1, N, LDA, m);
        arrayn_1 += ndmxt;
        arrayn += ndmxt;
        arraynp1 += ndmxt; 
    }
    return;
}

template<typename REAL>
template<class T>
void Chebyshev<REAL>::tracepolyA(
	T *ptr, void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
	std::complex<REAL> *wavein, 
	const int N, const int LDA, const int m)
{
    std::complex<REAL> *arraynp1;
	std::complex<REAL> *arrayn;
	std::complex<REAL> *arrayn_1;
    assert(N>=0 && LDA >= N);
    int ndmxt;
    if(m == 1) ndmxt = N * m;
    else       ndmxt = LDA * m;

    arraynp1 = new std::complex<REAL> [ndmxt];
    arrayn = new std::complex<REAL> [ndmxt];
    arrayn_1 = new std::complex<REAL> [ndmxt];

    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);

    (ptr->*funA)(arrayn_1, arrayn,m);

    polytrace[0] = this->ddot_real(wavein,wavein,N,LDA,m);
    polytrace[1] = this->ddot_real(wavein,arrayn,N,LDA,m);

    //more than 1-st orders
    for(int ior = 2; ior < norder; ++ior)
    {
        recurs_complex(ptr, funA, arraynp1, arrayn, arrayn_1, N, LDA, m);
        polytrace[ior] = this->ddot_real(wavein,arraynp1,N,LDA,m);   
        std::complex<REAL>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }

    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return;
}

template<typename REAL>
template<class T>
void Chebyshev<REAL>::recurs_complex(
    T *ptr, void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int),
	std::complex<REAL>* arraynp1,  //v_{n+1}
	std::complex<REAL>* arrayn,    //v_n
	std::complex<REAL>* arrayn_1,  //v_{n-1}
	const int N, const int LDA,  const int m)
{
    (ptr->*funA)(arrayn,arraynp1,m);
	for(int ib = 0 ; ib < m ; ++ib)
	{
    	for(int i = 0; i < N; ++i)
    	{
        	arraynp1[i+ib*LDA]=REAL(2.0)*arraynp1[i+ib*LDA]-arrayn_1[i+ib*LDA];
    	}
	}
}

template<typename REAL>
template<class T>
bool Chebyshev<REAL>::checkconverge(
	T *ptr, void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
 	std::complex<REAL> *wavein, const int N,
	REAL& tmax, 
	REAL& tmin, 
	REAL stept)
{
    bool converge = true;
    std::complex<REAL> *arraynp1;
	std::complex<REAL> *arrayn;
	std::complex<REAL> *arrayn_1;

    arraynp1 = new std::complex<REAL> [N];
    arrayn = new std::complex<REAL> [N];
    arrayn_1 = new std::complex<REAL> [N];

    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, N);

    if(tmin == tmax) 
	{
		tmax += stept;
	}

    (ptr->*funA)(arrayn_1, arrayn,1);
    REAL sum1,sum2;
    REAL t;
#ifdef __MPI
    sum1=ModuleBase::GlobalFunc::ddot_real(N,arrayn_1,arrayn_1);
    sum2=ModuleBase::GlobalFunc::ddot_real(N,arrayn_1,arrayn);
#else
    sum1=this->ddot_real(arrayn_1,arrayn_1,N);
    sum2=this->ddot_real(arrayn_1,arrayn,N);
#endif
    t = sum2 / sum1 * (tmax - tmin) / 2 + (tmax + tmin) / 2;
    if(t < tmin || tmin == 0)
    {
        converge = false;
        tmin = t-stept;
    }
    if(t > tmax)
    {
        converge = false;
        tmax = t+stept;
    }

    for(int ior = 2; ior < norder; ++ior)
    {
        (ptr->*funA)(arrayn,arraynp1,1);
#ifdef __MPI
        sum1=ModuleBase::GlobalFunc::ddot_real(N,arrayn,arrayn);
        sum2=ModuleBase::GlobalFunc::ddot_real(N,arrayn,arraynp1);
#else
    sum1=this->ddot_real(arrayn,arrayn,N);
    sum2=this->ddot_real(arrayn,arraynp1,N);
#endif
        t = sum2/sum1 * (tmax - tmin) / 2 + (tmax + tmin) / 2;
        if(t < tmin)
        {
            converge = false;
            tmin = t - stept;
        }
        else if(t > tmax)
        {
            converge = false;
            tmax = t + stept;
        }
        for(int i = 0; i < N; ++i)
        {
            arraynp1[i]=2.0*arraynp1[i]-arrayn_1[i];
        }
        std::complex<REAL>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem; 
    }
  
    delete [] arraynp1;
    delete [] arrayn;
    delete [] arrayn_1;
    return converge;
}

}