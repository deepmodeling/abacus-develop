#ifndef STO_CHEBYCHEV_H
#define STO_CHEBYCHEV_H
#include <complex>
#include "fftw3.h"

namespace ModuleBase
{
//template class for fftw
template<typename T>
class FFTW;

/**
 * @brief A class to treat the Chebyshev expansion.
 *
 * @author qianrui on 2022-05-18
 * @details
 * Math:
 * I. 
 * Chebyshev polynomial:
 * T_0(x) = 1;
 * T_1(x) = x;
 * T_2(x) = 2x^2 -1;
 * T_3(x) = 4x^3-3x;
 * T_{n+2}(x) = 2xT_{n+1}(x) - T_n(x)
 * II. 
 * Any analytical function f(x) can be expanded by Chebyshev polynomial:
 * f(x) = \sum_{n=0}^{norder-1} C_n[f]*T_n(x) (|x| < 1),
 * where C_n[f] = \frac{2-\delta_{n0}}{\pi} \int_0^\pi f(cos(\theta))cos(n\theta) d\theta
 * Here C_n can be calculate with FFT.
 * III. 
 * Any functions of linear Operator or matrix f(\hat{A}) or f(A) can also be expanded as well:
 * f(A) = \sum_{n=0}^{norder-1} C_n[f]*T_n(A) (|all eigenvalues of A| < 1).
 * f(A)v = \sum_{n=0}^{norder-1} C_n[f]*T_n(A)v, where v is column vector
 *       = \sum_{n=0}^{norder-1} C_n[f]*v_n, where v_n = T_n(A)v, v_0 = v
 * v_{n+2} = 2Av_{n+1} - v_n
 * IV.
 * v^+f(A)v = \sum_{n=0}^{norder-1} C_n[f]*v^+v_n = \sum_{n=0}^{norder-1} C_n[f] * w_n, 
 * where w_n = v^+ * v_n = v^+ * T_n(A) * v
 * 
 * USAGE: 
 * Chebyshev che(10); // constructe a chebyshev expansion of 10 orders (n=0,1,...,9)
 * 1. che.calcoef_real(&a, &A::cos) 						// calculate C_n[f], where f is a.cos
 *    for(int i=0;i<10;++i) cout<<che.coef_real[i]<<endl; 	//Then we print C_n[f]
 * 
 *    che.calcoef_complex(&b, &B::expi) 					// calculate C_n[g], where g is b.expi
 *    for(int i=0;i<10;++i) cout<<che.coef_complex[i]<<endl; //Then we print C_n[g]
 * 
 *    che.calcoef_pair(&c, &C::cos, &C::sin) 				// calculate C_n[g], where g is (c.cos, c.sin)
 *    for(int i=0;i<10;++i) cout<<che.coef_complex[i]<<endl; //Then we print C_n[g]
 * 
 * 2. che.calcoef_real(&occ, &Occupy::fd)
 * 	  che.calfinalvec_real(&hamilt, &Hamilt::hpsi, psi_in, psi_out, npw);
 *    //calculate f(H)|psi>, where f is occ.fd and H is hamilt.hpsi
 * 
 *    che.calcoef_complex(&b, &B::expi)
 * 	  che.calfinalvec_complex(&hamilt, &Hamilt::hpsi, psi_in, psi_out, npw, npwx, nbands);
 *    //calculate exp(iH)|psi_i>
 * 
 * 3. che.tracepolyA(&hamilt, &Hamilt::hpsi, psi_in, npw, npwx, nbands)
 * 	  //calculate \sum_i^{nbands} <psi_i|T_n(H)|psi_i>
 * 
 * 4. che.calcoef_complex(&fun, &toolfunc::expi);  //calculate C_n[exp(ix)]
 * 	  che.getpolyval(PI/4, T, norder);             //get T_n(pi/4)
 *    std::complex<double> sum(0,0);
 *    for(int i = 0; i < norder ; ++i)
 *    {
 * 		sum += che.coef_complex[i]*T[i];          //sum = exp(i*pi/4) = \sum_n C_n[exp(ix)]*T_n(pi/4)
 *    }
 * 
 * 5. che.recurs_complex(&hamilt, &Hamilt::hpsi, vp1, v, vm1, npw)
 *    //calculate vp1: |vp1> = 2 H|v> - |vm1>;
 * 
 */
template<typename REAL>
class Chebyshev
{

public:

    // constructor and deconstructor
    Chebyshev(const int norder);
    ~Chebyshev();

public:
	// I.
	// Calculate coefficients C_n[f], where f is a function of real number
	template<class T>
    void calcoef_real(T *ptr, REAL (T::*fun)(REAL));
	// Calculate coefficients C_n[g], where g is a function of complex number
	template<class T>
    void calcoef_complex(T *ptr, std::complex<REAL> (T::*fun)(std::complex<REAL>));
	// Calculate coefficients C_n[g], where g is a general complex function g(x)=(g1(x), g2(x)) e.g. exp(ix)=(cos(x), sin(x))
	template<class T>
	void calcoef_pair(T *ptr, REAL (T::*fun1)(REAL), REAL (T::*fun2)(REAL));
    
	// II.
	// Calculate the final vector f(A)v = \sum_{n=0}^{norder-1} C_n[f]*v_n
	// Here funA(in, out) means the map v -> Av : funA(v, Av)
	// Here m represents we treat m vectors at the same time: f(A)[v1,...,vm] and funA(in,out,m) means [v1,...,vm] -> A[v1,...,vm]
	// N is dimension of vector, and LDA is the distance between the first number of v_n and v_{n+1}.
	// LDA >= max(1, N). It is the same as the BLAS lib.
	// calfinalvec_real uses C_n[f], where f is a function of real number and A is a real Operator.
	template<class T>
    void calfinalvec(T *ptr, 
		void (T::*funA)(REAL *in, REAL *out, const int), 
		REAL *wavein, REAL *waveout, 
		const int N, const int LDA = 1,  const int m = 1); //do not define yet

	// calfinalvec_real uses C_n[f], where f is a function of real number and A is a complex Operator.
	template<class T>
    void calfinalvec_real(T *ptr, 
		void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
		std::complex<REAL> *wavein, std::complex<REAL> *waveout, 
		const int N, const int LDA = 1,  const int m = 1);

	// calfinalvec_complex uses C_n[g], where g is a function of complex number and A is a complex Operator.
	template<class T>
	void calfinalvec_complex(T *ptr,
		void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
		std::complex<REAL> *wavein, std::complex<REAL> *waveout, 
		const int N, const int LDA = 1,  const int m = 1);

	// III.
	// \sum_i v_i^+f(A)v_i = \sum_{i,n=0}^{norder-1} C_n[f]*v_i^+v_{i,n} = \sum_{n=0}^{norder-1} C_n[f] * w_n
	// calculate the sum of diagonal elements (Trace) of T_n(A) in v-represent: w_n = \sum_i v_i^+ * T_n(A) * v_i
	// i = 1,2,...m
	template<class T>
	void tracepolyA(
		T *ptr, void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
		std::complex<REAL> *wavein, 
		const int N, const int LDA = 1,  const int m = 1);
	
	// get T_n(x)
	void getpolyval(REAL x, REAL* polyval, const int N);

	// get each order of vector: {T_0(A)v, T_1(A)v, ..., T_n(A)v}
	// Note: use it carefully, it will cost a lot of memory!
	// calpolyvec_real: f(x) = \sum_n C_n*T_n(x), f is a real function
	template<class T>
	void calpolyvec_real(T *ptr, 
		void (T::*funA)(REAL *in, REAL *out, const int), 
		REAL *wavein, REAL *waveout, 
		const int N, const int LDA = 1,  const int m = 1);//do not define yet
	// calpolyvec_complex: f(x) = \sum_n C_n*T_n(x), f is a complex function
	template<class T>
	void calpolyvec_complex(T *ptr, 
		void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
		std::complex<REAL> *wavein, std::complex<REAL> *waveout, 
		const int N, const int LDA = 1,  const int m = 1);

	
	// IV.
	// recurs fomula: v_{n+1} = 2Av_n - v_{n-1}
	// get v_{n+1} from v_n and v_{n-1}
	// recurs_complex: A is a real operator
	template<class T>
    void recurs_real(
		T *ptr, void (T::*funA)(REAL *in, REAL *out, const int),
		REAL* arraynp1,  //v_{n+1}
		REAL* arrayn,    //v_n
		REAL* arrayn_1,  //v_{n-1}
		const int N, const int LDA = 1,  const int m = 1);
	// recurs_complex: A is a complex operator
	template<class T>
    void recurs_complex(
		T *ptr, void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int),
		std::complex<REAL>* arraynp1,  //v_{n+1}
		std::complex<REAL>* arrayn,    //v_n
		std::complex<REAL>* arrayn_1,  //v_{n-1}
		const int N, const int LDA = 1,  const int m = 1);
	//return 2xTn-Tn_1
	REAL recurs(const REAL x, const REAL Tn, const REAL Tn_1);
	
	// V.
	// auxiliary function
	// Abs of all eigenvalues of A should be less than 1.
	// Thus \hat(a) = \frac{(A - (tmax+tmin)/2)}{(tmax-tmin)/2}
	// tmax >= all eigenvalues; tmin <= all eigenvalues
	// Here we check if the trial number tmax(tmin) is the upper(lower) bound of eigenvalues and return it.
    template<class T>
	bool checkconverge(
		T *ptr, void (T::*funA)(std::complex<REAL> *in, std::complex<REAL> *out, const int), 
 		std::complex<REAL> *wavein, const int N,
		REAL& tmax, //trial number for upper bound
		REAL& tmin, //trial number for lower bound
		REAL stept); //tmax = max() + stept, tmin = min() - stept

public:
	//Members:
    int norder;   // order of Chebyshev expansion
    int norder2;  // 2 * norder * EXTEND

    REAL* coef_real; // expansion coefficient of each order
	std::complex<REAL>* coef_complex; // expansion coefficient of each order
	FFTW<REAL> fftw; //use for fftw
    REAL *polytrace; //w_n = \sum_i v^+ * T_n(A) * v

	bool getcoef_real;    //coef_real has been calculated
	bool getcoef_complex; //coef_complex has been calculated

private:
	//SI.
	//calculate dot product <psi_L|psi_R>
    REAL ddot_real(
    const std::complex<REAL>* psi_L,
    const std::complex<REAL>* psi_R,
	const int N, const int LDA = 1,  const int m = 1);

    
};

template<>
class FFTW<double>
{
public:
	FFTW(const int norder2_in);
	~FFTW();
	void execute_fftw();
    double* dcoef; //[norder2]
	fftw_complex *ccoef;
	fftw_plan coef_plan;
};

#ifdef __ENABLE_FLOAT_FFTW
template<>
class FFTW<float>
{
public:
	FFTW(const int norder2_in);
	~FFTW();
	void execute_fftw();
    float* dcoef; //[norder2]
	fftwf_complex *ccoef;
	fftwf_plan coef_plan;
};
#endif
}

#include "math_chebyshev_def.h"

#endif
